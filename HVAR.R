# -------------------------------------------------------
# HVAR wrapper
# -------------------------------------------------------

      HVAR <- function(VAR, max.lvl = NULL) {
        
        if (class(VAR) != "varest") stop("Please provide object of type 'varest'")
        type <- VAR$type
        y <- as.matrix(VAR$datamat[,c(1:VAR$K)])
        if (type == "none") x <- as.matrix(VAR$datamat[,-c(1:VAR$K)])
        if (type == "const") x <- as.matrix(VAR$datamat[,-c(1:VAR$K,ncol(VAR$datamat))])
        corr <- cor(x)
        cov <- cov(x)
        distmat <- ((1 - corr) / 2)^0.5
        distmat <- corr
        clustOrder <- hclust(dist(distmat), method = 'single')$order
        
        getRecBipart <- function(cov, sortIx, max.lvl = max.lvl) {
          
          w <- array(NA, c(ncol(cov) + if(type == "const") {1}else{0}, ncol(y), ceiling(log(ncol(cov))/log(2))))
          if (is.null(max.lvl)) max.lvl <- ceiling(log(ncol(cov))/log(2))
          
          # create recursion function within parent function to avoid use of globalenv
          recurFun <- function(cov, sortIx, res, max.lvl = max.lvl) {
            # get first half of sortIx which is a cluster order
            subIdx <- 1:trunc(length(sortIx)/2)
            
            # subdivide ordering into first half and second half
            cItems0 <- sortIx[subIdx]
            cItems1 <- sortIx[-subIdx]
            
            lvl <- min(which(is.na(w[c(cItems0[1]),1,])))
            
            if (lvl < max.lvl) {
              if (type == "const") {
                model <- lm(res ~ 1 + cbind(if (length(cItems0)==1) {x[,cItems0]} else {rowSums(x[,cItems0])}, 
                                            if (length(cItems1)==1) {x[,cItems1]} else {rowSums(x[,cItems1])}))
              } else {
                model <- lm(res ~ 0 + cbind(if (length(cItems0)==1) {x[,cItems0]} else {rowSums(x[,cItems0])}, 
                                            if (length(cItems1)==1) {x[,cItems1]} else {rowSums(x[,cItems1])}))
              }
            }
            if (lvl == max.lvl) {
              if (type == "const") {
                model <- lm(res ~ 1 + cbind(x[,cItems0], x[,cItems1]))
              } else {
                model <- lm(res ~ 0 + cbind(x[,cItems0], x[,cItems1]))
              }
            }
            
            coefs <- model$coefficients
            if (type == "const") coefs <- rbind(coefs[-1,], coefs[1,])
            
            if (type == "const") {
              w[c(cItems0,nrow(w)),,lvl] <<- matrix(coefs[1,], length(cItems0)+1, ncol(w), byrow = T)
              w[c(cItems1,nrow(w)),,lvl] <<- matrix(coefs[2,], length(cItems1)+1, ncol(w), byrow = T)
            } else {
              w[cItems0,,lvl] <<- matrix(coefs[1,], length(cItems0), ncol(w), byrow = T)
              w[cItems1,,lvl] <<- matrix(coefs[2,], length(cItems1), ncol(w), byrow = T)
            }
            
            res <- model$residuals / 2
            
            # rerun the function on a half if the length of that half is greater than 1
            if(length(cItems0) > 1 & lvl < max.lvl) {
              recurFun(cov, cItems0, res, max.lvl = max.lvl)
            }
            if(length(cItems1) > 1 & lvl < max.lvl) {
              recurFun(cov, cItems1, res, max.lvl = max.lvl)
            }
            
          }
          
          # run recursion function
          recurFun(cov, sortIx, y, max.lvl)
          final_w <- apply(w, c(1:2), sum, na.rm = T)
          return(final_w)
        }
        
        coefs <- getRecBipart(cov, clustOrder, max.lvl)
        
        for (i in 1:VAR$K) {
          VAR$varresult[[i]]$coefficients <- coefs[,i]
        }
        
        return(VAR)
      }
