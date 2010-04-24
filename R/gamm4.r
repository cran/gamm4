## Version of gamm using lme4 as fit engine. (c) Simon N. Wood 2009
## Reparameterization trick as Wood (2004,2006). 
## fooling lmer using Fabian Scheipl's trick.


gamm4.setup<-function(formula,pterms,data=stop("No data supplied to gamm4.setup"),knots=NULL)
# set up the model matrix, penalty matrices and auxilliary information about the smoothing bases
# needed for a gamm fit.
# There is an implicit assumption that any rank deficient penalty does not penalize 
# the constant term in a basis. 
{ 
  ## first simply call `gam.setup'....

  G <- mgcv:::gam.setup(formula,pterms,data=data,knots=knots,sp=NULL,
                    min.sp=NULL,H=NULL,absorb.cons=TRUE)
 
  if (!is.null(G$L)) stop("gamm4 can not handle linked smoothing parameters (probably from use of `id')")
  # now perform re-parameterization...

  first.f.para<-G$nsdf+1
  first.r.para<-1
 
  G$Xf <- G$X # full GAM model matrix, treating smooths as fixed effects
  random<-list()
  random.i<-0

  X <- G$X[,1:G$nsdf,drop=FALSE] # accumulate fixed effects into here

  xlab <- rep("",0)
  if (G$m)
  for (i in 1:G$m) 
  { sm <- G$smooth[[i]]
    sm$X <- G$X[,sm$first.para:sm$last.para]
 
    if (!sm$fixed) random.i <- random.i+1
   
    ZSZ <- list()
    if (!sm$fixed) { 
      if (length(sm$S)>1) stop("gamm4 can only handle single penalty smooths")
      ZSZ<-sm$S[[1]]
    }
    XZ<-sm$X
    k <- ncol(sm$X);j<-0
   
    if (!sm$fixed) 
    { sm$ZSZ <- ZSZ            # store these too - for construction of Vp matrix
      ev<-eigen(ZSZ,symmetric=TRUE)
      null.rank <- sm$df - sm$rank
      mult.pen <- FALSE

      p.rank <- ncol(XZ) - null.rank
      if (p.rank>ncol(XZ)) p.rank <- ncol(XZ)
      U<-ev$vectors
      D<-ev$values[1:p.rank]
      D<-1/sqrt(D)
      XZU<-XZ%*%U
      if (p.rank<k-j) Xf<-XZU[,(p.rank+1):(k-j),drop=FALSE]
      else Xf<-matrix(0,nrow(sm$X),0) # no fixed terms left
      
      Xr<-t(t(XZU[,1:p.rank])*D)
      n.para<-k-j-p.rank # indices for fixed parameters
      sm$first.f.para<-first.f.para
      first.f.para<-first.f.para+n.para
      sm$last.f.para<-first.f.para-1
      n.para<-ncol(Xr) # indices for random parameters
      sm$first.r.para<-first.r.para
      first.r.para<-first.r.para+n.para
      sm$last.r.para<-first.r.para-1
    
      sm$D<-D;sm$U<-U # information (with qrc) for backtransforming to original space 

      term.name <- paste("Xr.",random.i,sep="")
      term.name <- new.name(term.name,names(data))
      form <- as.formula(paste("~",term.name,"-1",sep=""))
     
      random[[random.i]] <- Xr
      names(random)[random.i] <- term.name
      attr(random[[random.i]],"s.label") <- sm$label
      sm$lmer.name <- term.name ## store the name by which this is identified in lmer call
      ## eval(parse(text=paste("G$",term.name,"<-Xr",sep="")))
    } else # term is fixed, so model matrix appended to fixed matrix
    { Xf <- XZ # whole term goes to fixed 
      n.para <- ncol(Xf)       # now define where the parameters of this term live 
      sm$first.f.para <- first.f.para
      first.f.para <- first.f.para+n.para
      sm$last.f.para <- first.f.para-1
    }
    ## now add appropriate column names to Xf.
    ## without these, summary.lme will fail
    
    if (ncol(Xf)) {
      Xfnames<-rep("",ncol(Xf)) 
      k<-length(xlab)+1
      for (j in 1:ncol(Xf)) {
        xlab[k] <- Xfnames[j] <-
        new.name(paste(sm$label,"Fx",j,sep=""),xlab)
        k <- k + 1
      } 
      colnames(Xf) <- Xfnames
    }

    X<-cbind(X,Xf) # add fixed model matrix to overall X
  
    sm$X <- NULL
  
    G$smooth[[i]] <- sm  ## replace smooth object with transformed version 
  }
 
  G$random<-random ## named list of Random effect matrices
  G$X<-X  ## fixed effects model matrix

  G
} ## end of gamm4 setup


gamm4 <- function(formula,random=NULL,family=gaussian(),data=list(),weights=NULL,
      subset=NULL,na.action,knots=NULL,...)
# Routine to fit a GAMM to some data. Fixed and smooth terms are defined in the formula, but the wiggly 
# parts of the smooth terms are treated as random effects. The onesided formula random defines additional 
# random terms. 

{ if (!require("lme4")) stop("gamm4() requires package lme4 to be installed")
  if (!require("mgcv")) stop("gamm4() requires package mgcv to be installed")
  if (!is.null(random)) {
    if (!inherits(random,"formula")) stop("gamm4 requires `random' to be a formula")
    random.vars <- all.vars(random)
  } else random.vars <- NULL

  # create model frame.....
  gp<-interpret.gam(formula) # interpret the formula 
  
  mf<-match.call(expand.dots=FALSE)
 
  mf$formula<-gp$fake.formula
  mf$family<-mf$scale<-mf$knots<-mf$random <- mf$...<-NULL ## mf$weights?
  mf$drop.unused.levels<-TRUE
  mf[[1]]<-as.name("model.frame")
  pmf <- mf
  gmf <- eval(mf, parent.frame()) # the model frame now contains all the data, for the gam part only 
  gam.terms <- attr(gmf,"terms") # terms object for `gam' part of fit -- need this for prediction to work properly

  if (length(random.vars)) {
    mf$formula <- as.formula(paste(paste(deparse(gp$fake.formula,
            backtick = TRUE), collapse = ""), "+", paste(random.vars,
            collapse = "+")))
    mf <- eval(mf, parent.frame())
  } else mf <- gmf
  rm(gmf)

  if (nrow(mf)<2) stop("Not enough (non-NA) data to do anything meaningful")
  Terms <- attr(mf,"terms")    
  
  ## summarize the *raw* input variables
  ## note can't use get_all_vars here -- buggy with matrices
  vars <- all.vars(gp$fake.formula[-2]) ## drop response here
  inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))
  dl <- eval(inp, data, parent.frame())
  names(dl) <- vars ## list of all variables needed
  var.summary <- mgcv:::variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data
  rm(dl) ## save space 

  pmf$formula <- gp$pf
  pmf <- eval(pmf, parent.frame()) # pmf contains all data for non-smooth part 
  
  pTerms <- attr(pmf,"terms")

  if (is.character(family)) family<-eval(parse(text=family))
  if (is.function(family)) family <- family()
  if (is.null(family$family)) stop("family not recognized")
  if (family$family == "gaussian" && family$link == "identity") linear <- TRUE else linear <- FALSE
  # now call gamm4.setup 

  G<-gamm4.setup(gp,pterms=pTerms,data=mf,knots=knots)
  
  G$var.summary <- var.summary    

  n.sr <- length(G$random) # number of random smooths (i.e. s(...,fx=FALSE,...) terms)

  if (is.null(random)&&n.sr==0) 
  stop("gamm4 models must have at least 1 smooth with unknown smoothing parameter or at least one other random effect")

  g<-as.factor(G$y*0+1) ## needed, whatever codetools says

  offset.name <- attr(mf,"names")[attr(attr(mf,"terms"),"offset")]

  yname <- new.name("y",names(mf))
  eval(parse(text=paste("mf$",yname,"<-G$y",sep="")))
  Xname <- new.name("X",names(mf))
  eval(parse(text=paste("mf$",Xname,"<-G$X",sep="")))
    
  lme4.formula <- paste(yname,"~",Xname,"-1")
  if (length(offset.name)) 
  { lme4.formula <- paste(lme4.formula,"+",offset.name) 
  }
  ## next add the random effect dummy variables for the smooth
  r.name <- names(G$random) 
  if (n.sr) for (i in 1:n.sr) # adding the constructed variables to the model frame avoiding name duplication
  { mf[[r.name[i]]] <- factor(rep(1:ncol(G$random[[i]]),length=nrow(G$random[[i]])))
    lme4.formula <- paste(lme4.formula,"+ (1|",r.name[i],")")
  }
  
  if (!is.null(random)) { ## append the regular random effects
    lme4.formula <- paste(lme4.formula,"+",substring(deparse(random),first=2))
  }
  
  lme4.formula <- as.formula(lme4.formula)
    
  if (linear) b <- lmer(lme4.formula,data=mf,family=family,weights=G$w,doFit=FALSE)
  else  b <- glmer(lme4.formula,data=mf,family=family,weights=G$w,doFit=FALSE)

  if (n.sr) { ## use Fabian Scheipl's trick of overwriting dummy slots
     tn <- names(b$FL$fl) 
     ## some names go with more than one element of b$FL$trms, so...
     tn <- tn[attr(b$FL$fl,"assign")] ## group name associated with each element of b$FL$trms 
     ind <- 1:length(tn)
     sn <- names(G$random)
     for (i in 1:n.sr) {
       k <- ind[sn[i]==tn] ## which trm should contain G$random[[i]] 
       b$FL$trms[[k]]$A <- b$FL$trms[[k]]$Zt <- as(t(G$random[[i]]),"dgCMatrix")
       attr(G$random[[i]],"s.label") -> sl
       attr(b$FL$trms[[k]]$ST,"dimnames") <- list(sl,sl)
     }
  }

  ret <- list()

  if (linear) ret$mer <- do.call(lme4:::lmer_finalize,b)
  else ret$mer <- do.call(lme4:::glmer_finalize,b)

  rm(b)

  

  ### .... fitting finished

  ## now fake a gam object 
    
  object<-list(model=mf,formula=formula,smooth=G$smooth,nsdf=G$nsdf,family=family,
                 df.null=nrow(G$X),y=ret$mer@y,terms=gam.terms,pterms=pTerms,xlevels=G$xlevels,
                 contrasts=G$contrasts,assign=G$assign,na.action=attr(mf,"na.action"),
                 cmX=G$cmX,var.summary=G$var.summary)
  
  ## to unpack coefficients look at names(ret$lme$flist), ret$lme@Zt, ranef(), fixef()


  # Transform  parameters back to the original space....
    bf <- as.numeric(lme4::fixef(ret$mer)) ## the fixed effects
    br <- lme4::ranef(ret$mer) ## a named list
    if (G$nsdf) p<-bf[1:G$nsdf] else p<-array(0,0)
    if (G$m>0) for (i in 1:G$m)
    { fx <- G$smooth[[i]]$fixed 
      first<-G$smooth[[i]]$first.f.para;last<-G$smooth[[i]]$last.f.para
      if (first <=last) beta<-bf[first:last] else beta<-array(0,0)
      if (fx) b <- beta 
      else # not fixed so need to undo transform of random effects etc. 
      { 
        b <- as.numeric(br[[G$smooth[[i]]$lmer.name]][[1]])     
       
        b <- c(G$smooth[[i]]$D*b,beta) # single penalty case
        b<-G$smooth[[i]]$U%*%b 
      }
      ## if (is.null(G$smooth[[i]]$C)) nc <- 0 else nc <- nrow(G$smooth[[i]]$C) 
      ## if (nc) b <- qr.qy(G$smooth[[i]]$qrc,c(rep(0,nc),b))
      object$smooth[[i]]$first.para<-length(p)+1
      p<-c(p,b)
      object$smooth[[i]]$last.para<-length(p)
    }
 
    object$coefficients<-p

    
    ## need to drop smooths from Zt and then
    ## form Z'phiZ + I \sigma^2
    vr <- lme4::VarCorr(ret$mer) ## list of ranef cov matrices in the same order as Gp

    scale <- as.numeric(attr(vr,"sc"))^2 ## get the scale parameter
    if (!is.finite(scale)) { 
      scale <- 1
      object$scale.estimated <- FALSE
    } else object$scale.estimated <- TRUE
    
    sp <- rep(-1,n.sr)

    Zt <- Matrix(0,0,ncol(ret$mer@Zt))
    if (n.sr==0) sn <- NULL ## names by which smooths are known in mer
    rn <- names(vr)
    for (i in 1:length(vr)) {
      if (is.null(sn)||!rn[i]%in%sn) { ## append non smooth r.e.s to Zt
        ind <- (ret$mer@Gp[i]+1):ret$mer@Gp[i+1]
        Zt <- rBind(Zt,ret$mer@Zt[ind,])
      } else if (!is.null(sn)) { ## extract smoothing parameters for smooth r.e.s
        k <- (1:n.sr)[rn[i]==sn] ## where in original smooth ordering is current smooth
        if (as.numeric(vr[[i]]>0)) sp[k] <- scale/as.numeric(vr[[i]]) else 
        sp[k] <- 1e10
      }
    }
    phi <- Matrix(0,nrow(Zt),nrow(Zt))
    for (i in 1:length(vr)) {
      k <- 0
      if (is.null(sn)||!rn[i]%in%sn) {
        nc <- ncol(vr[[i]])
        ind <- 1:nc + k;k <- k + nc
        phi[ind,ind] <- vr[[i]]
      } 
    }

   # weights <- NULL
   # if (is.null(weights)) object$prior.weights <- object$y*0+1 else object$prior.weights <- weights 
    
    object$prior.weights <- ret$mer@pWt

    if (length(ret$mer@var)==0) { 
      V <- Diagonal(ncol(Zt))*scale
      object$weights <- object$prior.weights
    } else 
    { V <- Diagonal(x=1/ret$mer@var)*scale ## the response variance conditional on the r.e.s
       object$weights <- 1/ret$mer@var
    }
    if (nrow(Zt)>0) V <- V + crossprod(Zt,phi%*%Zt) ## data or pseudodata cov matrix, treating smooths as fixed now
    G$Xf <- as(G$Xf,"dgCMatrix")
    XVX <- t(G$Xf)%*%solve(V,G$Xf) ## X'V^{-1}X
    
    object$sp <- sp

    S<-matrix(0,ncol(G$Xf),ncol(G$Xf)) # penalty matrix
    first <- G$nsdf+1
    k <- 1
    if (G$m>0) for (i in 1:G$m) # Accumulate the total penalty matrix
    { n.para <- object$smooth[[i]]$last.para - object$smooth[[i]]$first.para + 1
      last <- first + n.para - 1 
      if (!object$smooth[[i]]$fixed)
      { S[first:last,first:last] <- S[first:last,first:last] + 
                  object$smooth[[i]]$ZSZ*object$sp[k]
          k <- k+1
														      
      }
      first <- last + 1 
    }
   
    ## Vb <- solve(XVX+S/scale) # covariance matrix - in constraint space
    ev <- eigen(XVX+S/scale,symmetric=TRUE)
    ind <- ev$values != 0
    iv <- ev$values;iv[ind] <- 1/ev$values[ind]
    Vb <- ev$vectors%*%(iv*t(ev$vectors))

    object$edf<-rowSums(Vb*t(XVX))
   
    object$sig2 <- scale
    if (linear) { object$method <- "lmer.REML"
    } else { object$method <- "glmer.ML"}

    Vb <- Vb
    object$Vp <- as(Vb,"matrix")
  
    object$Ve <- as(Vb%*%XVX%*%Vb,"matrix")
 
    class(object) <- "gam"
    object$fitted.values <- predict.gam(object,type="response")
    object$residuals <- lme4:::residuals(ret$mer) 

    if (G$nsdf>0) term.names<-colnames(G$X)[1:G$nsdf] else term.names<-array("",0)
    n.smooth<-length(G$smooth) 
    if (n.smooth)
    for (i in 1:n.smooth)
    { k<-1
      for (j in object$smooth[[i]]$first.para:object$smooth[[i]]$last.para)
      { term.names[j]<-paste(object$smooth[[i]]$label,".",as.character(k),sep="")
        k<-k+1
      }
    }
    names(object$coefficients)<-term.names  # note - won't work on matrices!!
   

    object$gcv.ubre <- deviance(ret$mer)

    ret$gam<-object
    ret

} ## end of gamm4



print.gamm4.version <- function()
{ library(help=gamm4)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  cat(paste("This is gamm4 ",version,"\n",sep=""))
}

.onAttach <- function(...) { 
  print.gamm4.version()
 
}

.onUnload <- function(libpath) {}

.First.lib <- function(lib, pkg) {
  print.gamm4.version()
 
}
