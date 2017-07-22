# -----------------------------------------------------
# Simulator Function: var.sim
# -----------------------------------------------------

      pacman::p_load(tsDyn, truncnorm)

      var_sim <- function(sim.lags = 4, sim.vars = 4, len = 100, type = "const") {
        
        Bmat <- matrix(rtruncnorm(((sim.vars)^2)*sim.lags + if(type=="const") {sim.vars} else {0},-.9,.9, 0, 0.1),nrow = sim.vars)
        varcov<-matrix(rtruncnorm(sim.vars*sim.vars,0,.9, 0.2, 0.1),ncol = sim.vars)
        diag(varcov) <- 1
        varcov[lower.tri(varcov)] = t(varcov)[lower.tri(varcov)]
        varSim <-TVAR.sim(B=Bmat,nthresh=0,n=len, lag = sim.lags, include=type, varcov=varcov)
        return(varSim)
        
      }