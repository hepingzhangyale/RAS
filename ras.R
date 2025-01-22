library("segmented")
source("functions.R")


save.directory <- "./result"
read.directory <- "./genotype"

if (!dir.exists(save.directory)) {
  dir.create(save.directory, recursive = TRUE)
}

set.seed(123)


# Create dummy genotype matrix
geno <- matrix(sample(0:2, 100 * 1000, replace = TRUE, prob = c(0.8, 0.1, 0.1)), nrow = 100, ncol = 1000)
rownames(geno) <- paste("DUMMY", 1:100, sep = "")

# Create dummy confounder matrix
ID <- rownames(geno)
sex <- sample(c("Male", "Female"), 100, replace = TRUE)
age <- sample(20:60, 100, replace = TRUE)
age_squared <- age^2
age_sex <- age * as.numeric(sex == "Male")
pc1 <- rnorm(100)
pc2 <- rnorm(100)
pc3 <- rnorm(100)
pc4 <- rnorm(100)
pc5 <- rnorm(100)
pc6 <- rnorm(100)
pc7 <- rnorm(100)
pc8 <- rnorm(100)
pc9 <- rnorm(100)
pc10 <- rnorm(100)

confounder.mat <- data.frame(ID, sex, age, age_squared, age_sex, pc1, pc2, pc3, pc4, pc5, pc6, pc7, pc8, pc9, pc10)

# Create dummy phenotype matrix
is.continuous <- T ## Trait type

if (is.continuous){
  phenotype <- data.frame(ID = ID, phenotype = rnorm(100))
}else{
  phenotype <- data.frame(ID = ID, phenotype = rbinom(100,1, prob = 0.2))
}


## Configurations
this.chrome <- 1
N <- ncol(geno)  
n <- nrow(geno)
is.visualize <- T ## Visualization
num_rep <- 2 ## Screening Round
this.start <- 1 ## Starting SNP
this.skip1 <- 10 ## Pivotal SNP Density
this.skip2 <- 20 ## Window density 

this.seq <- seq(1,(N+1-this.start), by = this.skip1)
full.p.values <- rep(0, length(this.seq))



for(this.rep in 1:num_rep){
  cat("==============================\n","Begin Repetition ", this.rep, "\n")
  
  this.sample <- sample(dim(geno)[1], dim(geno)[1]/2, replace = F)
  this.leftout <- setdiff(seq(dim(geno)[1]), this.sample)

  
  if (T){ ## Data Splitting
    
    this.sample <- sort(this.sample)
        
    if (is.continuous){
      this.df <- merge(confounder.mat,phenotype, by = "ID")  
      lm0 <- lm(phenotype ~ sex + age + age_squared + age_sex + 
                  pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = this.df)
      
      phenotype$phenotype <- as.numeric(lm0$residuals)
      
      
      phenotype1 <- phenotype$phenotype[this.sample]
      phenotype2 <- phenotype$phenotype[this.leftout]
    }else{
      
      phenotype1 <- phenotype$phenotype[this.sample]
      phenotype2 <- phenotype$phenotype[this.leftout]
    }
    
  }
  
  ## Calculating effect sizes
  
  coef.mat <- matrix(NA, nrow = N,ncol = 4)
  cat("Starting GWAS ... \n")
  
  for(i in 1:N){
    if(i %% 10000 == 0){print(i)}
    this.x <- geno[this.sample,i]
    
    if (is.continuous){
      
      fit <- lm(phenotype1[is.na(this.x) == F] ~ na.omit(this.x))
      
    }else{
      
      this.df.sub <- this.df[this.sample,]
      this.df.sub$this.x <- this.x
      this.df.sub$phenotype1 <- phenotype1
      clean_df <- na.omit(this.df.sub)  # Removes all rows with any NA values in the specified columns
      
      fit <- lm(phenotype1 ~ this.x + sex + age + age_squared + age_sex +
                  pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, data = clean_df)
      
    }
    coef.mat[i,] <- summary(fit)$coefficients[2,]
  }
  
  colnames(coef.mat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  pgs.weights <- coef.mat[,1]

  saveRDS(coef.mat,paste0(save.directory,"/chr-",this.chrome,"_coef_mat-",this.rep,".rds"))
  

  geno[is.na(geno)] <- 0
  pgs.mat <- matrix(data = NA, nrow = length(this.leftout), ncol = dim(geno)[2])
  for (i in 1:length(this.leftout)){
    pgs.mat[i,] <- geno[this.leftout[i],] * pgs.weights
  }
  

  cat("Starting Scanning ... \n")
 

  this.df.sub <- this.df[this.sample,]
  this.df.sub$phenotype2 <- phenotype2
  
  return.p.values <- screen_forward_max_region(geno, pgs.mat, this.df.sub, -1,
                                   start.point = 1, save.directory, this.chrome,
                                   min_window_size = 5, max_window_size = 100,
                                  isSimulation = FALSE, this.repetition = this.rep, screening_round = 1, isPlot = F, skip1 = this.skip1,
                                  skip2 = this.skip2)

  
  full.p.values <- full.p.values + return.p.values
  mean.p.values <- full.p.values / this.rep

}

## Changepoint Detection

x <- seq(length(mean.p.values))
y <- mean.p.values
cat("Starting Detection ... \n")

this.result <- changePoint_detection_window_scanning(x, y, p.values.threshold = 1e-8, min.length = 10, skip = 1, 
                                                     window_size = 100, slope_check_window_size = 30, slope.p.values.threshold.left = 1e-10,
                                                     slope.p.values.threshold.right = 1e-10)
this.result$tau_hats <- sapply(this.result$tau_hats, get_local_maximum, y=y)


this.result <- second_scanning(this.result, second_window_size = 100, this.skip = this.skip1, p.value.threshold = 1e-15)



all.changepoints <- this.result$all.changepoints
all.p.values <- this.result$all.p.values

this.FP <- NULL
this.FN <- NULL

x <- this.seq

saveRDS(list(x = x, all.changepoints = all.changepoints, all.p.values = all.p.values, 
             tau_hats = this.result$tau_hats, left.slopes = this.result$left.slopes,
             right.slopes = this.result$right.slopes, FP = this.FP, FN = this.FN), 
        file = paste0(save.directory, "/final_detection_result-", this.chrome,".RDS"))



if(T){

  pdf(file = paste0(save.directory,"/chr-", this.chrome, "-cp-plot.pdf"))
    
  
  
  plot(x, y, type = 'l', main = "Detected Changepoint with p-value threshold = 1e-8")
  abline(v = this.result$tau_hats, col = "red", lty = 2)  
  
  dev.off()
}


## Visualization

if(is.visualize){


  pdf(file = paste0(save.directory,"/chr-", this.chrome, "-cp-p-values-plot.pdf"))
    
  par(mar=c(5, 4, 4, 5) + 0.1)  
  
  plot(x, y, type = 'l', main = "The Existence Test p-value at Potential Change Points", xlim = range(x), ylab = "RAS", xlab = "SNP Index",
       ylim = c(-max(y),max(y)))

  par(new=TRUE)
  
  plot(all.changepoints, all.p.values, col='red', pch=17, ylab='', xlab='', xaxt='n', yaxt='n', ylim=c(min(all.p.values), max(all.p.values)*2), xlim = range(x))

  axis(4)
  mtext('-log(p.values) (ChangePoint)', side=4, line=3)
  segments(x0=all.changepoints, y0=all.p.values, x1=all.changepoints, y1=0, col='red', lty=2)
  abline(h = 8,lty = 2, col = "blue")
  dev.off()

}


