# Fit tsBART.  Wrapper function for calling tsbartFit.cpp and tsbartProbit.cpp.


### For argument validation.
.ident <- function(...){
   # courtesy https://stackoverflow.com/questions/19966515/how-do-i-test-if-three-variables-are-equal-r
   args <- c(...)
   if( length( args ) > 2L ){
      #  recursively call ident()
      out <- c( identical( args[1] , args[2] ) , .ident(args[-1]))
   }else{
      out <- identical( args[1] , args[2] )
   }
   return( all( out ) )
}

### Tau solver for estimating prior probit scales (sd_control, sd_moderate) from data if indicated.
tau_calc = function(phat, rr){

   mu = qnorm(phat)
   f = function(x){rr-pnorm(mu+x)/pnorm(mu)}

   # Error handling for viable combination of phat and rr.
   mytry = try(uniroot(f, lower=-10, upper=10, extendInt = "yes"), silent=T)

   # If valid combination, return calculated tau.  Else print warning and calculate
   # tau for max rr possible given the input phat.
   if(is(mytry, "try-error")){

      # Print warning message.
      print('Warning: baseline risk (phat) and relative risk (rr) input combination invalid. Calculating tau based on max possible rr given phat.')

      # De-crement rr until it is possible to achieve rr with any tau, given phat.
      # Output this tau.
      err_flag = 1
      while(err_flag==1){
         rr = rr - .01
         mytry = try(uniroot(f, lower=-10, upper=10, extendInt = "yes"), silent=T)

         if(!is(mytry, "try-error")){
            tau = uniroot(f, lower=-10, upper=10, extendInt = "yes")$root
            err_flag = 0
            print(paste0('Max possible relative risk given phat: ', rr))
         }
      }

   } else{
      tau = uniroot(f, lower=-10, upper=10, extendInt = "yes")$root
   }

   return(tau)
}

### Main pbcf wrapper function.
pbcf <- function(y, pihat, z, x_control, x_moderate,
                  pihatpred=NULL, zpred=NULL, xpred_control=NULL, xpred_moderate=NULL,
                  nburn=100, nsim=1000, ntree_control=200, ntree_moderate=50,
                  sigq=.9, nu=3,
                  base_control=.95, power_control=2,
                  base_moderate=.25, power_moderate=3,
                  sd_control=2*sd(y), sd_moderate=sd(y),
                  use_muscale=T, use_tauscale=T,
                  pihat_in_trt=F,
                  set_probit_scales=F,
                  verbose=T, mh=F, save_inputs=T){

   ################################################################
   # Capture key arguments.
   ################################################################

   lambda = NULL  # Not needed in probit - to be deleted from package later on
   sighat = NULL  # Not needed in probit - to be deleted from package later on
   
   
   options(expressions=10000)

   inputs = cbind.data.frame(
      'arg'= c('nburn','nsim','ntree_control','ntree_moderate','lambda','sigq','sighat','nu','base_control','power_control',
              'base_moderate','power_moderate','sd_control','sd_moderate','use_muscale','use_tauscale',
              'pihat_in_trt','verbose','mh'),
      'value'= c(nburn,nsim,ntree_control,ntree_moderate,ifelse(is.null(lambda),"NULL",lambda),sigq,ifelse(is.null(sighat),"NULL",sighat),nu,base_control,power_control,
                base_moderate,power_moderate,sd_control,sd_moderate, use_muscale, use_tauscale,
                pihat_in_trt,verbose,mh)
   )

   ################################################################
   # Validate inputs.
   ################################################################

   #---------------------------------------------------------------
   # If not predicting, set pred objects to first three obs.
   # These are not output.
   #---------------------------------------------------------------
   predict = 1

   if(sum(is.null(pihatpred), is.null(zpred), is.null(xpred_control), is.null(xpred_moderate))>0){
      predict = 0
      pihatpred = pihat[1:min(3,length(pihat))]
      zpred = z[1:min(3,length(z))]
      xpred_control = x_control[1:min(3,nrow(x_control)),,drop=F]
      xpred_moderate = x_moderate[1:min(3,nrow(x_moderate)),,drop=F]
   }

   #---------------------------------------------------------------
   # Data size.
   #---------------------------------------------------------------

   # Check data size match.
   if( !.ident(length(y), length(pihat), nrow(x_control), nrow(x_moderate))){

      stop("Data size mismatch. The following should all be equal:
           length(y): ", length(y), "\n",
           "length(pihat): ", length(pihat), "\n",
           "nrow(x_control): ", nrow(x_control), "\n",
           "nrow(x_moderate): ", nrow(x_moderate), "\n")
   }

   # Check out-of-sample data size match.
   if( !.ident(nrow(xpred_control), nrow(xpred_moderate), length(pihatpred))){

      stop("Data size mismatch. The following should all be equal:
           length(pihatpred): ", nrow(pihatpred), "\n",
           "nrow(xpred_control): ", nrow(xpred_control), "\n",
           "nrow(xpred_control): ", nrow(xpred_control), "\n")
   }



   #---------------------------------------------------------------
   # Probit checks.
   #---------------------------------------------------------------

   # Data size match including yobs.
   if( !.ident(length(y), nrow(x_moderate), nrow(x_control))){
      stop("Data size mismatch. The following should all be equal:
           length(y): ", length(y), "\n",
           "nrow(x_control): ", nrow(x_control), "\n",
           "nrow(x_moderate): ", nrow(x_moderate), "\n")
   }

   # Yobs is only 0/1.  Y must not be only 0/1.
   yobs = y  # Observed 0/1
   y = ifelse(yobs == 1, 1.96, -1.96)
   
   if(length(unique(yobs))>2) warning("y should be binary 0/1")



   #---------------------------------------------------------------
   # Other inputs.
   #---------------------------------------------------------------

   if(any(is.na(yobs))) stop("Missing values in y")
   if(any(is.na(x_control))) stop("Missing values in x_control")
   if(any(is.na(x_moderate))) stop("Missing values in x_moderate")
   if(any(is.na(xpred_control))) stop("Missing values in xpred_control")
   if(any(is.na(xpred_moderate))) stop("Missing values in xpred_moderate")

   if(any(pihat>1 || pihat<0)) stop("pihat must be in 0-1 range.")
   if(any(pihatpred>1 || pihatpred<0)) stop("pihatpred must be in 0-1 range.")

   if(nburn<0) stop("nburn must be positive")
   if(nsim<0) stop("nsim must be positive")

   
   
   ################################################################
   # Scale/center response y. (For non-probit case.)
   ################################################################
   y_scale = y 

   ################################################################
   # Create model matrix and set up hyperparameters.
   ################################################################

   # Add propensity score to x_control and xpred_control matrices.
   x_control = cbind(x_control, pihat)
   xpred_control = cbind(xpred_control, pihatpred)

   # If indicated by user, add propensity score to x_moderate and xpred_moderate matrices.
   if(pihat_in_trt==TRUE){
      x_moderate[,ncol(x_moderate)+1] = pihat
      xpred_moderate[,ncol(xpred_moderate)+1] = pihatpred
   }

   
   # Model matrices.
   xx_control = makeModelMatrix(x_control)
   xxpred_control = makeModelMatrix(xpred_control)
   cutpoints_control = makeCutpoints(xx_control)

   xx_moderate = makeModelMatrix(x_moderate)
   xxpred_moderate = makeModelMatrix(xpred_moderate)
   cutpoints_moderate = makeCutpoints(xx_moderate)

   
   # Sighat and lambda calibration.
   df = cbind.data.frame(y_scale, xx_control)
   lmf = lm(y_scale ~ ., data=df)
   sighat = sigma(lmf)

   qchi = qchisq(1-sigq, nu)
   lambda = (sighat * sighat * qchi) / nu


   
   
   ################################################################
   # Set up probit parameters
   ################################################################
   phat = mean(unlist(yobs))
   offset = qnorm(phat)
   

   ################################################################
   # Probit-scaled default control_sd and moderate_sd, based on
   # estimates of baseline risk and relative risk in data, if indicated.
   ################################################################
   if (set_probit_scales==T) {

      phat = sum(yobs==1)/length(yobs)
      rrhat = (sum(yobs==1 & z==1) / sum(z==1)) /
         (sum(yobs==1 & z==0) / sum(z==0))

      sd_control = abs(qnorm(phat))
      sd_moderate = abs(tau_calc(phat, rrhat))

      print(paste0('Setting probit scales:'))
      print(paste0('sd_control: ', sd_control, ', sd_moderate: ', sd_moderate))

      # Update inputs.
      levels(inputs$value) <- c(levels(inputs$value), sd_control, sd_moderate)
      inputs$value[which(inputs$arg=="sd_control")] = sd_control
      inputs$value[which(inputs$arg=="sd_moderate")] = sd_moderate

   }

   ################################################################
   # Set up ordering, so that z=1 is first in in-samp and out-of-samp datasets.
   ################################################################

   perm = order(z, decreasing=TRUE)
   perm_oos = order(zpred, decreasing=TRUE)

 
   

   ################################################################
   # Call probitbcf.cpp
   ################################################################
   out = NULL

   out =  probitbcf(y = y_scale[perm], yobs=yobs[perm], z = z[perm], zpred = zpred[perm_oos],
                      x_con = t(xx_control[perm,]), x_mod = t(xx_moderate[perm,]),
                      xpred_con = t(xxpred_control[perm_oos,]), xpred_mod = t(xxpred_moderate[perm_oos,]),
                      xinfo_list_con = cutpoints_control, xinfo_list_mod = cutpoints_moderate,
                      nburn = nburn, nsim = nsim, ntree_con = ntree_control, ntree_mod = ntree_moderate,
                      offset=offset,
                      lambda=lambda, sigq=sigq, nu=nu,
                      base_con=base_control, power_con=power_control,
                      base_mod=base_moderate, power_mod=power_moderate,
                      con_sd=sd_control, mod_sd=sd_moderate,
                      use_muscale=use_muscale, use_tauscale=use_tauscale,
                      treef_name_="tsbtrees.txt", save_trees=FALSE, silent_mode=!verbose, trt_init = 1)

   
   
   ################################################################
   # Rescale output and restore to correct order.
   ################################################################

   # In-sample
   yhat = out$yhat[,order(perm)]
   mu = out$mu[,order(perm)]
   tau = out$tau[,order(perm)]
   sig_rescaled = out$sigma

   # Out-of-sample
   yhat_oos = out$yhat_oos[,order(perm_oos)]
   mu_oos = out$mu_oos[,order(perm_oos)]
   tau_oos = out$tau_oos[,order(perm_oos)]

   # Set up mu_sd and tau_sd draws.
   mu_sd = abs(sd_control * out$eta)
   tau_sd = abs(sd_moderate * (out$bscale1-out$bscale0))

   ################################################################
   # Adjust alpha's and add accept/reject indicator.
   # Note: bd function returns:
   #     alpha in (0,1) for births which are accepted
   #     -alpha in (-1,0) for deaths which are accepted
   #     10+alpha for rejected births
   #     -alpha-10 for rejected deaths.
   ################################################################

   if(mh){
      # Con metropolis info.
      bd_con = ifelse(out$alpha_con<=0,0,1)             # 1 = birth, 0 = death
      accept_con = ifelse(abs(out$alpha_con)<10,1,0)   # 1 = accepted, 0 = rejected
      alpha_con = ifelse(accept_con==1, abs(out$alpha_con), abs(out$alpha_con)-10)

      # Mod metropolis info.
      bd_mod = ifelse(out$alpha_mod<=0,0,1)             # 1 = birth, 0 = death
      accept_mod = ifelse(abs(out$alpha_mod)<10,1,0)   # 1 = accepted, 0 = rejected
      alpha_mod = ifelse(accept_mod==1, abs(out$alpha_mod), abs(out$alpha_mod)-10)

      # Assemble dataframes and convert bd to character.
      # Note: alpha is (nburn+nsim x ntree).
      metrop_con = cbind.data.frame(
         'iter' = rep(1:(nburn+nsim), times=ntree_control),
         'tree' = rep(1:ntree_control, each=nburn+nsim),
         'accept' = as.numeric(accept_con),
         'alpha' = as.numeric(alpha_con),
         'bd' = as.numeric(bd_con)
      )

      metrop_mod = cbind.data.frame(
         'iter' = rep(1:(nburn+nsim), times=ntree_moderate),
         'tree' = rep(1:ntree_moderate, each=nburn+nsim),
         'accept' = as.numeric(accept_mod),
         'alpha' = as.numeric(alpha_mod),
         'bd' = as.numeric(bd_mod)
      )

      metrop_con$bd = ifelse(metrop_con$bd==1,'birth','death')
      metrop_mod$bd = ifelse(metrop_mod$bd==1,'birth','death')

      # Combine into one metrop dataframe.
      metrop = rbind.data.frame(metrop_con, metrop_mod)
      metrop$tree = c(rep("control",nrow(metrop_con)), rep("moderate",nrow(metrop_mod)))
      rm(metrop_con); rm(metrop_mod)
   }

   ################################################################
   # Return output.
   ################################################################

   # Only include out-of-sample info if indicated by user.
   if(predict){
      out = list('yhat'=yhat,
                 'mu'=mu,
                 'tau'=tau,
                 'yhat_oos'=yhat_oos,
                 'mu_oos'=mu_oos,
                 'tau_oos'=tau_oos,
                 'sigma'=sig_rescaled,
                 'mu_sd' = mu_sd,
                 'tau_sd' = tau_sd,
                 'bscale1'=out$bscale1,
                 'bscale0' = out$bscale0,
                 'eta' = out$eta,
                 'gamma' = out$gamma_mu)
   } else{
      out = list('yhat'=yhat,
                 'mu'=mu,
                 'tau'=tau,
                 'sigma'=sig_rescaled,
                 'mu_sd' = mu_sd,
                 'tau_sd' = tau_sd,
                 'bscale1'=out$bscale1,
                 'bscale0' = out$bscale0,
                 'eta' = out$eta,
                 'gamma' = out$gamma_mu)
   }

   # Include metropolis info if indicated.
   if(mh){
      out$metrop = metrop
   }

   # Include inputs if indicated.
   if(save_inputs){
      out$inputs = inputs
   }

   # Return output.
   return(out)
   }
