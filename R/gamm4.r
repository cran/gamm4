## Version of gamm using lme4 as fit engine. (c) Simon N. Wood 2009-20
## Reparameterization trick as Wood (2004,2006). 
## fooling lmer using Fabian Scheipl's trick (now adapted for lme4 >1.0).

gamm4.setup<-function(formula,pterms,
                      data=stop("No data supplied to gamm.setup"),knots=NULL)
## set up the model matrix, penalty matrices and auxilliary information about the smoothing bases
## needed for a gamm4 fit.
## There is an implicit assumption that any rank deficient penalty does not penalize 
## the constant term in a basis. 
## 1. Calls gam.setup, as for a gam to produce object G suitable for estimating a gam.
## 2. Works through smooth list, G$smooth, modifying so that... 
##    i) Smooths are reparameterized to have a sequence of (portion of) identity matrix
##       penalties. 
##    ii) 'random' list is accumulated containing random effect model matrices for terms.     
##    iii) Sparse version of full model matrix in original parameterization is also accumulated
##    iv) Various indices are created for moving between the parameterizations.
{ 
  ## first simply call `gam.setup'....

  G <- mgcv:::gam.setup(formula,pterms,
                 data=data,knots=knots,sp=NULL,
                 min.sp=NULL,H=NULL,absorb.cons=TRUE,sparse.cons=0,gamm.call=TRUE)
 
  if (!is.null(G$L)) stop("gamm can not handle linked smoothing parameters (probably from use of `id' or adaptive smooths)")
  # now perform re-parameterization...

  first.f.para <- G$nsdf+1 
  
  random <- list()
  
  if (G$nsdf>0) ind <- 1:G$nsdf else ind <- rep(0,0)  
  X <- G$X[,ind,drop=FALSE] # accumulate fixed effects into here

  xlab <- rep("",0)
  
  G$Xf <- as(X,"dgCMatrix") ## sparse version of full matrix, treating smooths as fixed

  first.para <- G$nsdf+1

  used.names <- names(data) ## keep track of all variable names already used

  if (G$m) for (i in 1:G$m) { ## work through the smooths
    
    sm <- G$smooth[[i]]
    sm$X <- G$X[,sm$first.para:sm$last.para,drop=FALSE]
    rasm <- mgcv::smooth2random(sm,used.names,type=2) ## convert smooth to random effect and fixed effects
    used.names <- c(used.names,names(rasm$rand))    

    sm$fixed <- rasm$fixed

    ## deal with creation of sparse full model matrix  
    if (!is.null(sm$fac)) { 
      flev <- levels(sm$fac) ## grouping factor for smooth
      n.lev <- length(flev)
      for (k in 1:n.lev) {
        G$Xf <- cbind2(G$Xf,as(sm$X*as.numeric(sm$fac==flev[k]),"dgCMatrix"))
      }
    } else { 
      n.lev <- 1
      G$Xf <- cbind2(G$Xf,as(sm$X,"dgCMatrix"))
    }

    ## now append random effects to main list
    n.para <- 0 ## count random coefficients
    #rinc <- rind <- rep(0,0)
    if (!sm$fixed) {
      for (k in 1:length(rasm$rand)) n.para <- n.para + ncol(rasm$rand[[k]])
      sm$lmer.name <- names(rasm$rand)
      random <- c(random,rasm$rand)
      sm$trans.D <- rasm$trans.D
      sm$trans.U <- rasm$trans.U ## matrix mapping fit coefs back to original
    }

    ## ensure stored first and last para relate to G$Xf in expanded version

    sm$last.para <- first.para + ncol(rasm$Xf) + n.para - 1
    sm$first.para <- first.para
    first.para <- sm$last.para + 1    

    if (ncol(rasm$Xf)) {
      Xfnames <- rep("",ncol(rasm$Xf)) 
      k <- length(xlab)+1
      for (j in 1:ncol(rasm$Xf)) {
        xlab[k] <- Xfnames[j] <-
        new.name(paste(sm$label,"Fx",j,sep=""),xlab)
        k <- k + 1
      } 
      colnames(rasm$Xf) <- Xfnames
    }

    X <- cbind(X,rasm$Xf) # add fixed model matrix to overall fixed X
   
    sm$first.f.para <- first.f.para
    first.f.para <- first.f.para + ncol(rasm$Xf)
    sm$last.f.para <- first.f.para - 1 ## note less than sm$first.f.para => no fixed

    ## store indices of random parameters in smooth specific array
    sm$rind <- rasm$rind 
    sm$rinc <- rasm$rinc 

    sm$pen.ind <- rasm$pen.ind ## pen.ind==i TRUE for coef penalized by ith penalty

    sm$n.para <- n.para
 
    sm$X <- NULL ## delete model matrix
  
    G$smooth[[i]] <- sm  ## replace smooth object with extended version 
  }
 
  G$random <- random ## named list of random effect matrices
  G$X <- X  ## fixed effects model matrix

  G
} ## end of gamm4.setup


gamm4 <- function(formula,random=NULL,family=gaussian(),data=list(),weights=NULL,
      subset=NULL,na.action,knots=NULL,drop.unused.levels=TRUE,REML=TRUE,
      control=NULL,start=NULL,verbose=0L,nAGQ=1L,python_cholmod=FALSE,...) {
# Routine to fit a GAMM to some data. Fixed and smooth terms are defined in the formula, but the wiggly 
# parts of the smooth terms are treated as random effects. The onesided formula random defines additional 
# random terms. 


  if (!is.null(random)) {
    if (!inherits(random,"formula")) stop("gamm4 requires `random' to be a formula")
    random.vars <- all.vars(random)
  } else random.vars <- NULL

  # create model frame.....
  gp <- interpret.gam(formula) # interpret the formula 
  
  mf <- match.call(expand.dots=FALSE)
 
  mf$formula <- gp$fake.formula
  mf$REML <- mf$verbose <- mf$control <- mf$start <- mf$family <- mf$scale <-
             mf$knots <- mf$random <- mf$nAGQ <- mf$python_cholmod <- mf$... <-NULL ## mf$weights?
  mf$drop.unused.levels <- drop.unused.levels
  mf[[1]] <- as.name("model.frame")
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
  
  ## summarize the *raw* input variables
  ## note can't use get_all_vars here -- buggy with matrices
  vars <- all.vars(gp$fake.formula[-2]) ## drop response here
  inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))
  dl <- eval(inp, data, parent.frame())
  names(dl) <- vars ## list of all variables needed
  var.summary <- mgcv:::variable.summary(gp$pf,dl,nrow(mf)) ## summarize the input data

  ## lmer offset handling work around...
  mvars <- vars[!vars%in%names(mf)] ## variables not in mf raw -- can cause lmer problem
  if (length(mvars)>0) for (i in 1:length(mvars)) mf[[mvars[i]]] <- dl[[mvars[i]]] ## append raw versions to mf

  rm(dl) ## save space 

  pmf$formula <- gp$pf
  pmf <- eval(pmf, parent.frame()) # pmf contains all data for non-smooth part 
  pTerms <- attr(pmf,"terms")

  if (is.character(family)) family<-eval(parse(text=family))
  if (is.function(family)) family <- family()
  if (is.null(family$family)) stop("family not recognized")
  if (family$family == "gaussian" && family$link == "identity") linear <- TRUE else linear <- FALSE
  # now call gamm4.setup 

  G <- gamm4.setup(gp,pterms=pTerms,data=mf,knots=knots)
  
  G$var.summary <- var.summary    

  n.sr <- length(G$random) # number of random smooths (i.e. s(...,fx=FALSE,...) terms)

  if (is.null(random)&&n.sr==0) 
  stop("gamm4 models must have at least 1 smooth with unknown smoothing parameter or at least one other random effect")

  offset.name <- attr(mf,"names")[attr(attr(mf,"terms"),"offset")]

  yname <- new.name("y",names(mf))
  eval(parse(text=paste("mf$",yname,"<-G$y",sep="")))
  Xname <- new.name("X",names(mf))
  eval(parse(text=paste("mf$",Xname,"<-G$X",sep="")))
    
  lme4.formula <- paste(yname,"~",Xname,"-1")
  if (length(offset.name)) 
  { lme4.formula <- paste(lme4.formula,"+",offset.name) 
  }

  ## Basic trick is to call (g)lFormula to set up model, with simple i.i.d. dummy random effects for the 
  ## penalized component of each smooth. This results in columns of Z being produced for these dummy's,
  ## which can be over-written with the right thing. NOTE: that lambdat could also be modified, I think!!

  ## Add the random effect dummy variables for the smooth
  r.name <- names(G$random) 
  if (n.sr) for (i in 1:n.sr) # adding the constructed variables to the model frame avoiding name duplication
  { mf[[r.name[i]]] <- factor(rep(1:ncol(G$random[[i]]),length=nrow(G$random[[i]])))
    lme4.formula <- paste(lme4.formula,"+ (1|",r.name[i],")")
  }
  
  if (!is.null(random)) { ## append the regular random effects
    lme4.formula <- paste(lme4.formula,"+",
                    substring(paste(deparse(random,backtick=TRUE),collapse=""),first=2))
  }
  
  lme4.formula <- as.formula(lme4.formula)

  if (is.null(control)) control <- if (linear) lmerControl() else glmerControl() 

  ## NOTE: further arguments should be passed here... 
  b <- if (linear) lFormula(lme4.formula,data=mf,weights=G$w,REML=REML,control=control,...) else 
                   glFormula(lme4.formula,data=mf,family=family,weights=G$w,control=control,...)

    ## loop through random effect smooths and ingest them into Z
  if (n.sr) {
    tn <- names(b$reTrms$cnms) ## names associated with columns of Z (same order as Gp)
    ind <- 1:length(tn)
    sn <- names(G$random)

    sparse_summary <- as.data.frame(summary(b$reTrms$Zt))
    sparse_dims <- dim(b$reTrms$Zt)

    for (i in 1:n.sr) {
      k <- ind[sn[i]==tn] ## which term (variable) name represents random smooth i

      # Step 1: Extract indices and values from the transposed matrix
      indices <- which(t(G$random[[i]]) != 0, arr.ind = TRUE)
      values <- t(G$random[[i]])[indices]
      
      # Step 2: Create dense summary DataFrame for the smooth submatrix
      dense_summary <- data.frame(i = indices[, 1], j = indices[, 2], x = values)
      
      # Step 3: Adjust row indices according to b$reTrms$Gp[k]
      dense_summary$i <- dense_summary$i + b$reTrms$Gp[k]
      #ii <- (b$reTrms$Gp[k]+1):b$reTrms$Gp[k+1]
      #message("Debug: ii ", ii)
      #message("Debug: summary ", dense_summary$i)

      
      # Step 4: Remove rows from sparse_summary that have the same i as in dense_summary
      filtered_sparse_summary <- sparse_summary[!(sparse_summary$i %in% dense_summary$i), ]
      
      # Step 5: Combine the filtered sparse summary with the dense summary
      combined_summary <- rbind(filtered_sparse_summary, dense_summary)

      # Step 6: Update sparse_summary
      sparse_summary <- combined_summary
      
      # Step 7: Update labels
      b$reTrms$cnms[[k]] <- attr(G$random[[i]],"s.label") 
      }

    # Create the new sparse matrix
    new_sparse_mat <- sparseMatrix(
          i = sparse_summary$i,
          j = sparse_summary$j,
          x = sparse_summary$x,
          dims = c(sparse_dims[1], sparse_dims[2])
        )
        
    # Assign row and column names from the original matrix
    rownames(new_sparse_mat) <- rownames(b$reTrms$Zt)
    colnames(new_sparse_mat) <- colnames(b$reTrms$Zt)

    # Replace old Zt with the updated one
    b$reTrms$Zt <- new_sparse_mat
  }

  ## now do the actual fitting...
  ret <- list()
  #arg <- list(...)
  #arg <- arg[!(names(arg) %in% names(b))]
  #b <- c(b,arg) ## add '...' arguments for use with do.call
  b$control <- control; b$verbose=verbose; b$start=start
  if (linear) {
    ## Create the deviance function to be optimized:
    devfun <- do.call(mkLmerDevfun, b)
    ## Optimize the deviance function:
    opt <- optimizeLmer(devfun,start=start,verbose=verbose,control=control$optCtrl,optimizer=control$optimizer[[1]],calc.derivs=control$calc.derivs,boundary.tol=control$boundary.tol)
    ret$mer <- mkMerMod(environment(devfun), opt, b$reTrms, fr = b$fr)
  } else { ## generalized case...
    ## Create the deviance function for optimizing over theta:
    devfun <- do.call(mkGlmerDevfun, b)
    ## Optimize over theta using a rough approximation (i.e. nAGQ = 0):
    opt <- optimizeGlmer(devfun,start=start,verbose=verbose,control=control$optCtrl,optimizer=control$optimizer[[1]],calc.derivs=control$calc.derivs,boundary.tol=control$boundary.tol,nAGQ=0)
    ## Update the deviance function for optimizing over theta and beta:
    devfun <- updateGlmerDevfun(devfun, b$reTrms)
    ## Optimize over theta and beta:
    opt <- optimizeGlmer(devfun, stage=2,start=start,verbose=verbose,control=control$optCtrl,optimizer=control$optimizer[[1]],calc.derivs=control$calc.derivs,boundary.tol=control$boundary.tol,nAGQ=nAGQ)
    ## Package up the results:
    ret$mer <- mkMerMod(environment(devfun), opt, b$reTrms, fr = b$fr)
  }

  rm(b)

  ### .... fitting finished

  ## now fake a gam object 
    
  object<-list(model=mf,formula=formula,smooth=G$smooth,nsdf=G$nsdf,family=family,
                 df.null=nrow(G$X),y=getME(ret$mer,"y"),
                 terms=gam.terms,pterms=G$pterms,xlevels=G$xlevels,
                 contrasts=G$contrasts,assign=G$assign,na.action=attr(mf,"na.action"),
                 cmX=G$cmX,var.summary=G$var.summary)
  pvars <- all.vars(delete.response(object$terms))
  object$pred.formula <- if (length(pvars)>0) reformulate(pvars) else NULL

  ## to unpack coefficients look at names(ret$lme$flist), ret$lme@Zt, ranef(), fixef()
 
    ## let the GAM coefficients in the original parameterization be beta,
    ## and let them be beta' in the fitting parameterization. 
    ## Then beta = B beta'. B and B^{-1} can be efficiently accumulated
    ## and are useful for stable computation of the covariance matrix
    ## etc... 
  
    B <- Matrix(0,ncol(G$Xf),ncol(G$Xf))
    diag(B) <- 1
    Xfp <- G$Xf
    ## Transform  parameters back to the original space....
    bf <- as.numeric(lme4::fixef(ret$mer)) ## the fixed effects
    br <- lme4::ranef(ret$mer, condVar=FALSE) ## a named list
    if (G$nsdf) p <- bf[1:G$nsdf] else p <- array(0,0) ## fixed parametric componet
    if (G$m>0) for (i in 1:G$m) {
      fx <- G$smooth[[i]]$fixed 
      first <- G$smooth[[i]]$first.f.para; last <- G$smooth[[i]]$last.f.para
      if (first <=last) beta <- bf[first:last] else beta <- array(0,0)
      if (fx) b <- beta else { ## not fixed so need to undo transform of random effects etc. 
        b <- rep(0,0)
        for (k in 1:length(G$smooth[[i]]$lmer.name)) ## collect all coefs associated with this smooth
          b <- c(b,as.numeric(br[[G$smooth[[i]]$lmer.name[k]]][[1]]))     
        b <- b[G$smooth[[i]]$rind] ## make sure coefs are in order expected by smooth
        b <- c(b,beta) 
        b <- G$smooth[[i]]$trans.D*b
        if (!is.null(G$smooth[[i]]$trans.U)) b <- G$smooth[[i]]$trans.U%*%b ## transform back to original 
      }
      p <- c(p,b)
     
      ## now fill in B...
      ind <- G$smooth[[i]]$first.para:G$smooth[[i]]$last.para
      if (!fx) { 
         D <- G$smooth[[i]]$trans.D
         if (is.null(G$smooth[[i]]$trans.U)) B[ind,ind] <- Diagonal(length(D),D) else
         B[ind,ind] <- t(D*t(G$smooth[[i]]$trans.U))
      }
      ## and finally transform G$Xf into fitting parameterization...
      Xfp[,ind] <- as.matrix(G$Xf[,ind,drop=FALSE]%*%B[ind,ind,drop=FALSE])
      #rep <- G$Xf[,ind,drop=FALSE]%*%B[ind,ind,drop=FALSE]

      # sparse summaries
      #sum_rep <- summary(rep)
      #sum_Xfp <- summary(Xfp)

      # replace column indexes in rep
      #sum_rep$j <- ind[sum_rep$j]

      # remove columns from orig
      #sum_Xfp <- subset(sum_Xfp, !(j %in% ind))

      # merge
      #sum_comb <- rbind(sum_Xfp, sum_rep)

      # recreate updated sparse matrix
      #Xfp <- sparseMatrix(
      #  i = sum_comb$i,
      #  j = sum_comb$j,
      #  x = sum_comb$x,
      #  dims = dim(Xfp)
      #)
    }
 
    object$coefficients <- p

    ## need to drop smooths from Zt and then
    ## form Z'phiZ + I \sigma^2

    vr <- lme4::VarCorr(ret$mer) ## list of ranef variance components in the same order as Gp
    
    scale <- as.numeric(attr(vr,"sc"))^2 ## get the scale parameter
    if (!is.finite(scale) || scale==1) { ## NOTE: better test???
      scale <- 1
      object$scale.estimated <- FALSE
    } else object$scale.estimated <- TRUE
    
    sp <- rep(-1,n.sr)

    Zt <- Matrix(0,0,ncol(getME(ret$mer,"Zt")))
    if (n.sr==0) sn <- NULL ## names by which smooths are known in mer
    rn <- names(vr)
    ind <- rep(0,0) ## index the non-smooth random effects among the random effects
    for (i in 1:length(vr)) {
      if (is.null(sn)||!rn[i]%in%sn) { ## append non smooth r.e.s to Zt
        Gp <- getME(ret$mer,"Gp") ## group index ends
        ind <- c(ind,(Gp[i]+1):Gp[i+1])
      } else if (!is.null(sn)) { ## extract smoothing parameters for smooth r.e.s
        k <- (1:n.sr)[rn[i]==sn] ## where in original smooth ordering is current smoothing param
        if (as.numeric(vr[[i]]>0)) sp[k] <- scale/as.numeric(vr[[i]]) else 
        sp[k] <- 1e10
      }
    }

    if (length(ind)) { ## extract columns corresponding to non-smooth r.e.s 
      Zt <- getME(ret$mer,"Zt")[ind,] ## extracting random effects model matrix
      root.phi <- getME(ret$mer,"Lambdat")[ind,ind] ## and corresponding sqrt of cov matrix (phi)
    }

    object$prior.weights <- G$w
                          
    if (linear) {
      object$weights <- object$prior.weights 
      V <- Diagonal(n=length(object$weights),x=scale/object$weights) 
    } else { 
     # mu <- getME(ret$mer,"mu")
     # eta <- family$linkfun(mu)
      object$weights <- ret$mer@resp$sqrtWrkWt()^2
      ## object$prior.weights*family$mu.eta(eta)^2/family$variance(mu)
      V <- Diagonal(x=1/object$weights)*scale
      #V <- Diagonal(x=scale*family$variance(mu)/object$prior.weights)
    }

  
    if (nrow(Zt)>0) V <- V + crossprod(root.phi%*%Zt)*scale ## data or pseudodata cov matrix, treating smooths as fixed now

    ## NOTE: Cholesky probably better in the following - then pivoting 
    ##       automatic when solving....

    if (python_cholmod) {
      library(reticulate)
      py_require(packages="scikit-sparse>=0.4.16")
      cholmod <- import("sksparse.cholmod")
      R <- cholmod$cholesky(V)
    } else {
      R <- mgcv::mchol(V);piv <- attr(R,"pivot")
    }

    G$Xf <- as(G$Xf,"dgCMatrix")
    Xfp <- as(Xfp,"dgCMatrix")
    
    if (python_cholmod) {
      WX <- as.matrix(R$solve_Lt(Xfp,use_LDLt_decomposition=FALSE))
      XVX <- as.matrix(R$solve_Lt(G$Xf,use_LDLt_decomposition=FALSE))
    } else {

      if (is.null(piv)) {
      WX <- as(solve(t(R),Xfp),"matrix")    ## V^{-.5}Xp -- fit parameterization
      XVX <- as(solve(t(R),G$Xf),"matrix")  ## same in original parameterization

      } else {
      WX <- as(solve(t(R),Xfp[piv,]),"matrix")    ## V^{-.5}Xp -- fit parameterization
      XVX <- as(solve(t(R),G$Xf[piv,]),"matrix")  ## same in original parameterization
      }
    }

    qrz <- qr(XVX,LAPACK=TRUE)
    object$R <- qr.R(qrz);object$R[,qrz$pivot] <- object$R

    XVX <- crossprod(object$R) ## X'V^{-1}X original parameterization

    object$sp <- sp
    
    colx <- ncol(G$Xf)
    Sp <- matrix(0,colx,colx) # penalty matrix - fit param
    first <- G$nsdf+1
    k <- 1
    if (G$m>0) for (i in 1:G$m) { # Accumulate the total penalty matrix
      if (!object$smooth[[i]]$fixed) {
        ii <- object$smooth[[i]]$first.para:object$smooth[[i]]$last.para ## index this smooth's params
        for (j in 1:length(object$smooth[[i]]$S)) { ## work through penalty list
          ind <- ii[object$smooth[[i]]$pen.ind == j] ## index of currently penalized
          diag(Sp)[ind] <-  sqrt(object$sp[k]) ## diagonal penalty
          k <- k+1
        }														              }
      first <- last + 1 
    }
   
    ## Alternative cov matrix calculation. Basic
    ## idea is that cov matrix is computed stably in
    ## fitting parameterization, and then transformed to
    ## original parameterization. 
    qrx <- qr(rbind(WX,Sp/sqrt(scale)),LAPACK=TRUE)
    Ri <- backsolve(qr.R(qrx),diag(ncol(WX)))
    ind <- qrx$pivot;ind[ind] <- 1:length(ind)## qrx$pivot
    Ri <- Ri[ind,] ## unpivoted square root of cov matrix in fitting parameterization Ri Ri' = cov
    Vb <- B%*%Ri; Vb <- Vb%*%t(Vb)

    object$edf<-rowSums(Vb*t(XVX))
   
    object$df.residual <- length(object$y) - sum(object$edf)

    object$sig2 <- scale
    if (linear) { object$method <- "lmer.REML"
    } else { object$method <- "glmer.ML"}

    object$Vp <- as(Vb,"matrix")
  
    object$Ve <- as(Vb%*%XVX%*%Vb,"matrix")
   
    class(object) <- "gam"
   
    ## Restore original smooth list, if it was split to deal with t2 terms...
    if (!is.null(G$original.smooth)) {
      object$smooth <- G$smooth <- G$original.smooth
    }

    ## If prediction parameterization differs from fit parameterization, transform now...
    ## (important for t2 smooths, where fit constraint is not good for component wise 
    ##  prediction s.e.s)

    if (!is.null(G$P)) {
      object$coefficients <- G$P %*% object$coefficients
      object$Vp <- G$P %*% object$Vp %*% t(G$P) 
      object$Ve <- G$P %*% object$Ve %*% t(G$P) 
    }

    object$linear.predictors <- predict.gam(object,type="link")
    object$fitted.values <- object$family$linkinv(object$linear.predictors)
    
    object$residuals <- residuals(ret$mer) 

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
    names(object$coefficients) <- term.names  # note - won't work on matrices!!
    names(object$edf) <- term.names
    names(object$sp) <- names(G$sp)

    object$gcv.ubre <- if (isREML(ret$mer)) REMLcrit(ret$mer) else deviance(ret$mer)

    if (!is.null(G$Xcentre)) object$Xcentre <- G$Xcentre ## any column centering applied to smooths

    ret$gam<-object
    class(ret) <- c("gamm4","list")
    ret

} ## end of gamm4


print.gamm4.version <- function() {
  if (requireNamespace("gamm4", quietly = TRUE)) {
    version_info <- packageDescription("gamm4")$Version
  } else {
    # fallback: try to read DESCRIPTION manually (for devtools::load_all)
    desc_path <- file.path(dirname(sys.frame(1)$ofile), "..", "DESCRIPTION")
    if (file.exists(desc_path)) {
      dcf <- read.dcf(desc_path)
      version_info <- dcf[1, "Version"]
    } else {
      version_info <- "unknown"
    }
  }
  hello <- paste("This is gamm4 ", version_info, "\n", sep = "")
  packageStartupMessage(hello)
}

.onAttach <- function(...) { 
  print.gamm4.version()
 
}

.onUnload <- function(libpath) {}