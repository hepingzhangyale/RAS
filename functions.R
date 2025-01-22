screen_forward_max_region <- function(geno, pgs.mat, this.df, num_signals,
                               start.point = 1, save.directory = "simulation", this.chrome = 1,
                               min_window_size = 5, max_window_size = 100,
                               isSimulation = T, this.repetition = 1, screening_round = 1, isPlot = T, skip1 = 100, skip2 = 5){
  
  N <- ncol(geno)  
  n <- nrow(geno)

  overall.min.pvalues <- c()
  this.seq <- seq(1,N, by = skip1)

  p.values <- rep(NA,screening_round*(length(this.seq)))
  

  
  for (kk in 1:screening_round){

    big.count <- 0
    
    
    for (this.start in this.seq){
      big.count <- big.count + 1
      
      sub.seq <- seq(min_window_size, max_window_size , by = this.skip2)
      sub.seq <- c(0, sub.seq)
      if (sub.seq[length(sub.seq)] != max_window_size){
        sub.seq <- c(sub.seq, max_window_size)
      }
      
      sub.p.values <- rep(NA,length(sub.seq))
      
      count <- 0
      this.left.0 <- this.start
      this.right.0 <- this.start

      for (ws in sub.seq){
        count <- count + 1
        

        
        if (kk == 1 && count == 1){
          
          this.pgs <- pgs.mat[,(this.start)]

        }else{
          
       
          this.left.1 <- this.start - ws + 1
          this.left.1[this.left.1 < 1] <- 1
          
          this.right.1 <- this.start + ws - 1
          this.right.1[this.right.1 > N] <- N
          

          if (this.left.1 != this.left.0){
            this.left.seq <- this.left.1:(this.left.0-1)
              
            for (j in this.left.seq){
              this.pgs <- this.pgs + pgs.mat[,j]

            }
          }
          
          if (this.right.1 != this.right.0){
            this.right.seq <- (this.right.0+1):this.right.1
            
            for (j in this.right.seq){
              this.pgs <- this.pgs + pgs.mat[,j]
              
            }
          }          
          
          this.left.0 <- this.left.1
          this.right.0 <- this.right.1

          
          
        }
        
        this.df$this.pgs <- this.pgs
        
        if (is.continuous){ 
          fit <- lm(phenotype2 ~ this.pgs + this.x + sex + age + age_squared + age_sex +
                      pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = this.df) 
          
        }else{
          fit <- glm(phenotype2 ~ this.pgs + this.x + sex + age + age_squared + age_sex +
                       pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, family = binomial(), data = this.df) 
        }        
        
        if (nrow(summary(fit)$coefficients) > 1){
          sub.p.values[(kk - 1)*length(this.seq) + count] <- summary(fit)$coefficients[2,4]
        }else{
          print("warning")
          sub.p.values[(kk - 1)*length(this.seq) + count] <- 1
        }
        
        
        
        
      }
      
      p.values[big.count] <- min(sub.p.values)
      
      
    }
    
    
    
    
    
  }
  
  
  if(isPlot){
    if (isSimulation){
      
      
      
      pdf(file = paste0(save.directory,"/forward-",this.repetition,"-",kk,".pdf"))
      
      
      plot(this.seq,-log(p.values[((kk-1)*length(this.seq) + 1):(kk*length(this.seq))],base = 10), ylab = "-log(p.values)", type = 'l',
           xlab = "Position", main =paste0("Forward, Round ",kk))
      if(num_signals >= 1){
        abline(v = signal.starts, col = "red",  lty = 2)
        abline(v = signal.starts+signal.window.size-1, col = "red",  lty = 2)
      }
      
    }else{
      
      
      pdf(file = paste0(save.directory,"/BPM_chr",this.chrome,"-",this.start,"_forward-",this.repetition,"-",kk,".pdf"))
      
      plot(this.seq,-log(p.values[((kk-1)*length(this.seq)  + 1):(kk*length(this.seq))],base = 10), ylab = "-log(p.values)", type = 'l',
           xlab = "Position", main =paste0("BPM Chrome-",this.chrome, " Starting from ", this.start, "-- Forward, Round ",kk))
      
      
    }
    dev.off()
    
  }
  
  
  return(-log(p.values,base = 10))
}





# For changepoint detection
get_break_points <- function(x,y,t){
  this.df <- data.frame(y = y[1:t], x = x[1:t])
  fit_lm <- lm(y ~ x, data = this.df)  # intercept-only model

  fit_segmented <- try({
    suppressWarnings(segmented(fit_lm, seg.Z = ~x, npsi = 1))
  }, silent = TRUE)
  
  if (inherits(fit_segmented, "try-error")) {
    return(list(break.points = NULL, p.values = 1, slope.left = NULL, slope.right = NULL))
    
  }
  
  if (is.null(fit_segmented$psi)){
    return(list(break.points = NULL, p.values = 1, slope.left = NULL, slope.right = NULL))
  }
  
  this.slopes <- slope(fit_segmented)$x[,1]
  
  tau_hat <- round(fit_segmented$psi[1,2])
  
  if (F){
    
    this.p.value <- get_p_value(x,y,t, tau_hat)
    
  }else{
    
    this.test <- davies.test(fit_lm)
    this.p.value <- this.test$p.value
  }
  

  
  return(list(break.points = tau_hat, p.values = this.p.value, slope.left = this.slopes[1], slope.right = this.slopes[2]))
}


slope_test <- function(x,y,lower.tail){
  model <- lm(y ~ x-1)
  # lower.tail = F: test H0: beta1 >= 0, H1: beta1 < 0
  # Summary of the model to see coefficients
  # summary(model)
  
  # Perform one-tailed t-test on the slope (beta1)
  # H0: beta1 >= 0, H1: beta1 < 0
  
  t_statistic <- coef(summary(model))["x", "t value"]
  
  p_value <- pt(t_statistic, df = df.residual(model), lower.tail = lower.tail)
  p_value
}



# Changepoint detection algorithm

changePoint_detection_window_scanning <- function(x, y, p.values.threshold = 0.01, min.length = 10, skip = 1, window_size = 3000,
                                                  slope_check_window_size = 30, slope.p.values.threshold = 1e-8,
                                                  slope.p.values.threshold.left = 1e-10, slope.p.values.threshold.right = 1e-20){
  tau_hats <- c()
  p.values <- c()
  slope.left <- c()
  
  slope.right <- c()
  end.index <- c()
  all.changepoints <- c()
  all.p.values <- c()
  slope.angle <- c()
  start.t <- min.length - window_size

  while(T){

    
    this.seq <- seq(start.t, length(x) - min.length + 1, by = skip)

    if (this.seq[length(this.seq)] != length(x) - min.length + 1){
      this.seq <- c(this.seq, length(x) - min.length + 1)
    }
    
    for (start.t in this.seq){
      t <- window_size
      
      if (start.t < 0){
        t <- window_size + start.t
        start.t <- 1
      }
      
      if (start.t > (length(x) - window_size + 1)){
        t <- length(x) - start.t + 1
      }
      
        

      this.result <- get_break_points(x[start.t:length(x)],y[start.t:length(x)],t)

      all.changepoints <- c(all.changepoints, this.result$break.points)
      all.p.values <- c(all.p.values, this.result$p.values)
      
      v1 <- this.result$slope.left > 0
      v2 <- this.result$slope.right < 0

      
      if (this.result$p.values <= p.values.threshold){
        
        if(v1&& v2){
          
          this.tau <- this.result$break.points

          this.lower.index <- max(1, this.tau-slope_check_window_size)
          this.upper.index <- min(length(x),this.tau + slope_check_window_size)
          this.pvalue1 <- slope_test(x[this.lower.index:this.tau] -x[this.tau], y[this.lower.index:this.tau] - y[this.tau], F)
          
          this.pvalue2 <- slope_test(x[this.tau:this.upper.index] -x[this.tau], y[this.tau:this.upper.index] - y[this.tau], T)
          
          
          if (this.pvalue1 < slope.p.values.threshold.left && this.pvalue2 < slope.p.values.threshold.right){
            cat("ChangePoint Detected at ", this.result$break.points, " with p-value", 
                this.result$p.values, ", Left-slope ",this.result$slope.left, "Right-slope ", 
                this.result$slope.right, " \n")
            cat(this.pvalue1, " ", this.pvalue2, "\n")
            tau_hats <- c(tau_hats, this.result$break.points)
            p.values <- c(p.values, this.result$p.values)
            slope.left <- c(slope.left, this.result$slope.left)

            slope.right <- c(slope.right, this.result$slope.right)
            
            angle_degrees1 <- (atan(this.result$slope.left) * 180) / pi
            angle_degrees2 <- (atan(this.result$slope.right) * 180) / pi
            theta_degrees <- angle_degrees2 - angle_degrees1 + 180
            slope.angle <- c(slope.angle, theta_degrees)
            
            
            
            end.index <- c(end.index, start.t + t - 1)
            
            start.t <- this.result$break.points
            
            break
          }else{
            all.p.values[length(all.p.values)] <- 1
            
          }
        }else{
          all.p.values[length(all.p.values)] <- 1
          
        }
      }
    }
    
    if (start.t == this.seq[length(this.seq)]){
      break
      
    }
    
  }
  return(list(tau_hats = tau_hats, p.values = p.values, slope.left = slope.left, slope.right = slope.right, all.changepoints = all.changepoints,
              all.p.values = all.p.values, slope.angle = slope.angle, previous_tau_hats = tau_hats))
}


second_scanning <- function(this.result, this.start = 1, this.skip = 30, second_window_size = 50, p.value.threshold = 1e-10){
  
  unique_all.changepoints <- unique(this.result$all.changepoints)
  min_all.p.values <- sapply(unique_all.changepoints, function(x) min(this.result$all.p.values[this.result$all.changepoints == x]))
  tau_hats <- c()
  left.slopes <- c()
  right.slopes <- c()
  
  all.changepoints <- unique_all.changepoints
  all.p.values <- min_all.p.values
  
  all.p.values <- -log(all.p.values, base = 10)

  this.count <- 0
  
  for (tau_hat in this.result$tau_hats){
    
    this.count <- this.count + 1
    
    this.df <- data.frame(y = y[tau_hat:min((tau_hat + second_window_size),length(y))], 
                          x = x[tau_hat:min((tau_hat + second_window_size),length(x))])

    
    fit_lm <- lm(y ~ x, data = this.df)  # intercept-only model

    if(dim(this.df)[1] >= 4){
      this.test1 <- davies.test(fit_lm)
    }else{
      this.test1 <- list()
      this.test1$p.value <- 1.0
    }

    this.df <- data.frame(y = y[max((tau_hat - second_window_size),1):tau_hat], 
                          x = x[max((tau_hat - second_window_size),1):tau_hat])
    
    fit_lm <- lm(y ~ x, data = this.df)  # intercept-only model

    if(dim(this.df)[1] >= 4){
      this.test2 <- davies.test(fit_lm)
    }else{
      this.test2 <- list()
      this.test2$p.value <- 1.0
    }
    
    cat(this.test1$p.value, this.test2$p.value, "\n")
    
    
    if (this.test1$p.value < p.value.threshold || this.test2$p.value < p.value.threshold){
      fit_segmented <- try({
            suppressWarnings(segmented(fit_lm, seg.Z = ~x, npsi = 1))
      }, silent = TRUE)
      
      left.slopes <- c(left.slopes, this.result$slope.left[which(tau_hat == this.result$tau_hats)])
      right.slopes <- c(right.slopes, this.result$slope.right[which(tau_hat == this.result$tau_hats)])
      tau_hats <- c(tau_hats, tau_hat)
      
    }else{
      all.p.values[which(tau_hat == all.changepoints)] <- 0
      all.p.values[which(this.result$previous_tau_hats[this.count] == all.changepoints)] <- 0
      
    }
    
  }



    this.remove <- which(y[tau_hats] <= 2.5)
  
  
  if (length(this.remove) >= 1){
        tau_hats <- tau_hats[-this.remove]
  }

  all.changepoints <- this.start + (all.changepoints - 1)*this.skip
  tau_hats <- this.start + (tau_hats - 1)*this.skip
  return(list(all.changepoints = all.changepoints, tau_hats = tau_hats, all.p.values = all.p.values,
              left.slopes = left.slopes, right.slopes = right.slopes))
}


get_signal_points <- function(num_signals,signal.starts,signal.window.size){
  signal.points <- c()
  if (num_signals >= 1){
    for (i in 1:num_signals){
      signal.points <- c(signal.points, signal.starts[i], signal.starts[i] + signal.window.size)
    }
  }
  return(signal.points)
}


detect_FP <- function(tau_hats, signal.points, tol = 100){
  FP <- 0
  if(length(tau_hats) >= 1){
    for (i in 1:(length(tau_hats))){
      FP <- FP + 1
      if (length(signal.points) >= 1){
        for (j in 1:(length(signal.points)/2)){
          this.diff <- abs(tau_hats[i] - c(signal.points[2*j-1]:signal.points[2*j]))
          if (min(this.diff) <= tol){
            FP <- FP - 1
            break
          }
          
        }
      }
      
      
    }
  }
  return(FP)
}





detect_FN <- function(tau_hats, signal.points, tol = 100){
  FN <- 0
  if(length(signal.points) >= 1){
    for (i in 1:(length(signal.points)/2)){
      FN <- FN + 1
      if (length(tau_hats) >= 1){
        for (j in 1:(length(tau_hats))){
          this.diff <- abs(tau_hats[j] - c(signal.points[2*i-1]:signal.points[2*i]))
          if (min(this.diff) <= tol){
              FN <- FN - 1
              break
          }
          
          
          
        }
      }
      
      
    }
  }
  return(FN)
}


get_local_maximum <- function(y, x0, window.size = 50){
  
  if (F){
    
  if (y[x0 - 1] < y[x0] && y[x0 + 1] < y[x0]) {
    return(x0)
  }
  
  if(y[x0 - 1] > y[x0]){
    while(T){
      if (x0 == 1){
        return(x0)
      }
      
      x0 <- x0-1
      if (y[x0 - 1] < y[x0] && y[x0 + 1] < y[x0]) {
        return(x0)
      }
      

    }
    
    
  }
  
  if(y[x0 + 1] > y[x0]){
    while(T){
      if (x0 == length(y)){
        return(x0)
      }
      
      x0 <- x0 + 1
      
      if (y[x0 - 1] < y[x0] && y[x0 + 1] < y[x0]) {
        return(x0)
      }
      
      
    }
    
    
  }
  
  }else{
    
    lower.bound <- max(1, x0 - window.size)
    upper.bound <- min(length(y), x0 + window.size)
    
    x0 <- which.max(y[lower.bound:upper.bound]) + lower.bound - 1
    return(x0)
    
  }
  
  
  
  
}




