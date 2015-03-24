predict.graf.raster <-
  function (object, x, type, CI, maxn, ...) {
    
    # the following adapted from dismo:::predict
    
    variables <- colnames(m$x)
    
    if (is.null(CI)) {
      
      # if CIs aren't required
      out <- raster(x)
      
    } else {
      
      # if CIs are required
      out <- brick(x,
                   nl = 3)
      
    }
    
    ncols <- ncol(out)
    
    firstrow <- 1
    firstcol <- 1
    
    if (!canProcessInMemory(out, 3)) {
      
      filename <- rasterTmpFile()
      out <- writeStart(out, filename = filename, ...)
      inMemory <- FALSE
      
    } else {
      
      v <- array(NA, dim = c(nrow(out), ncol(out), nlayers(out)))
      inMemory <- TRUE
      
    }
    
    tr <- blockSize(out, n = nlayers(x) + 2)
    pb <- pbCreate(tr$n)
    
    for (i in 1:tr$n) {
      
      rr <- firstrow + tr$row[i] - 1
      
      rowvals <- getValuesBlock(x,
                                row = rr,
                                nrows = tr$nrows[i], 
                                firstcol,
                                ncols)
      
      rowvals <- rowvals[, variables, drop = FALSE]
      
      res <- matrix(NA,
                    nrow = nrow(rowvals),
                    ncol = nlayers(out))
      
      rowvals <- na.omit(rowvals)
      
      if (length(rowvals) > 0) {
        
        newdata <- data.frame(rowvals)
        
        # loop through and convert factors to factors
        for (col in ncol(newdata)) {
          
          if (is.factor(m$obsx[, col])) {
            
            newdata[, col] <- factor(newdata[, col])
            
          }
          
        }
        
        # predict to this batch of pixels
        p <- predict.graf(m,
                          newdata = newdata,
                          type = type,
                          CI = CI,
                          maxn = maxn,
                          ...)
        
        naind <- as.vector(attr(rowvals, "na.action"))
        
        if (!is.null(naind)) {
          
          res[-naind, ] <- p
          
        } else {
          
          res <- p
          
        }
        
      }
      
      if (inMemory) {
        
        res <- array(res, dim = c(nrow(out), tr$nrows[i], nlayers(out)))
        cols <- tr$row[i]:(tr$row[i] + dim(res)[2] - 1)
        v[, cols, ] <- res
        
      } else {
        
        out <- writeValues(out, res, tr$row[i])
        
      }
      
      pbStep(pb, i)
      
    }
    
    pbClose(pb)
    
    if (inMemory) {
      
      for (i in 1:nlayers(out)) {
      
        out <- setValues(out,
                         as.vector(v[, , i]),
                         layer = i)
  
      }
      
    } else {
      
      out <- writeStop(out)
      
    }
    
    names(out) <- colnames(p)
    
    return (out)
    
  }