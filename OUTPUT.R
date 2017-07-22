# --------------------------------------------------------
# OUTPUT
# --------------------------------------------------------

      pacman::p_load(vars)
      
      source("HVAR.R")

      # Load Simulators / Datasets
      source("sim_var.R")
      
      # Settings
      sim.runs <- 100
      sim.lags <- 4
      lags <- 6
      vars <- 10
      type <- "none"
      smpl <- 500
      fcst <- 0.05
      
      
      result <- matrix(NA, ncol = 5, nrow = sim.runs)
      for (i in 1:nrow(result)) {
        
        # Dataset
          # var.sim
          #dta <- var_sim(sim.lags, vars, len = smpl, type = type)
          # Randomly select "vars" different shares from a dataset of 20 (daily returns)
          x <- sample(c(1:20),vars)
          rows <- sample(1:(2615-smpl), 1)
          dta <- data.frame(read.csv("rets.csv"))[c(rows:(rows+smpl)),x+1] 
        
        dta.sample <- dta[1:ceiling(nrow(dta)*(1-fcst)),]
        dta.out <- dta[-c(1:ceiling(nrow(dta)*(1-fcst))),]
        
        suppressAll(var.conv <- VAR(dta.sample, p = lags, type = type))
        var.hrp2 <- HVAR(var.conv, max.lvl = 2)
        var.hrp3 <- HVAR(var.conv, max.lvl = 3)
        var.hrp4 <- HVAR(var.conv, max.lvl = 4)
        var.hrp <- HVAR(var.conv, max.lvl = NULL)
        
        fcst.conv <- sapply(predict(var.conv, n.ahead = nrow(dta.out))$fcst, function(x) x[,"fcst"])
        fcst.hrp2 <- sapply(predict(var.hrp2, n.ahead = nrow(dta.out))$fcst, function(x) x[,"fcst"])
        fcst.hrp3 <- sapply(predict(var.hrp3, n.ahead = nrow(dta.out))$fcst, function(x) x[,"fcst"])
        fcst.hrp4 <- sapply(predict(var.hrp4, n.ahead = nrow(dta.out))$fcst, function(x) x[,"fcst"])
        fcst.hrp <- sapply(predict(var.hrp, n.ahead = nrow(dta.out))$fcst, function(x) x[,"fcst"])
        
        result[i,] <- c(sum(sqrt(colMeans((fcst.conv - dta.out)^2))),
                        sum(sqrt(colMeans((fcst.hrp2 - dta.out)^2))),
                        sum(sqrt(colMeans((fcst.hrp3 - dta.out)^2))),
                        sum(sqrt(colMeans((fcst.hrp4 - dta.out)^2))),
                        sum(sqrt(colMeans((fcst.hrp - dta.out)^2))))
      }
      
      ggplot() + geom_density(aes(result[,1], color = "VAR")) + 
        geom_density(aes(result[,2], color = "HVAR2")) + 
        geom_density(aes(result[,3], color = "HVAR3")) + 
        geom_density(aes(result[,4], color = "HVAR4")) + 
        geom_density(aes(result[,5], color = "HVAR")) + 
        labs(x = "RMFE") + xlim(0,.5)
      