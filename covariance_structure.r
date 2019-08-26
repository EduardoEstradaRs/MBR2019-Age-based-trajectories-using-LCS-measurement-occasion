# R code for the statistical models used in:

# Estrada, Hamagami, & Ferrer (2019).
# Estimating age-based developmental trajectories using latent change score models
# based on measurement occasion.
# Multivariate Behavioral Reseach.
# https://doi.org/10.1080/00273171.2019.1647822

## RECOVERING THE AGE-BASED VARIANCE-COVARIANCE STRUCTURE

# Load required libraries
if (!"package:minpack.lm" %in% search()) {require(minpack.lm)}
if (!"package:lavaan" %in% search()) {require(lavaan)}

# Load data set
d <- myData
# d is a data frame in wide format. Variables of interest:
# age1: Age at first measurement occasion (value of 0 indicates 5 years of age)
# oy1 oy2 oy3: observed values for Y at the three time points
# cohort: Cohort or grouping variable



## Data pre-processing --------------------------------------
# Nonlinear transformation of age1
d$age1NL <- d$age1 ** (1/2) # Model 8
#d$age1NL <- d$age1 ** (2)   # Model 9
#d$age1NL <- exp(-1*d$age1)  # Model 12
#...

# Model 13 (two-step solution)
model13 <- nlsLM(data=d, formula = oy1 ~ a - (a-b)*exp(c*age1*-1),
                 start = list(a = 30, b = 10, c=.2),
                 control = nls.control(maxiter = 1000, minFactor = 1/2048))

d$hatY1 <- predict(model13, newdata=data.frame(age1 = d$age1) )




## Specification of LCS-SEM for recovering THE COVARIANCE STRUCTURE ------------------
# Specification of models
oneStepModelMG = "
# Declaring latent level
y1	=~ 1* oy1
y2	=~ 1* oy2
y3	=~ 1* oy3

# Auto-regression
y2	~ 1* y1
y3	~ 1* y2

# Define latent change
dy2	 =~ 1* y2
dy3	 =~ 1* y3

# Auto-proportions
dy2	 ~ b_y * y1
dy3	 ~ b_y * y2

# Latent intercept and slope
y0 =~ 1 * y1
yA =~
1*dy2  + 1*dy3

# Latent means, variances and covariances
y0 ~ y0mn * 1
yA ~ yAmn * 1

y0 ~~ y0    ## This is the only parameter varying across cohorts
yA ~~ yAV   * yA
y0 ~~ y0Acv* yA # Covariance

# Errors
dy2  ~~ 0 * dy2
dy3  ~~ 0 * dy3

# Covariance between dynamic errors
y1  ~~ 0 * y1
y2  ~~ 0 * y2
y3  ~~ 0 * y3

oy1  ~~ MerY * oy1
oy2  ~~ MerY * oy2
oy3  ~~ MerY * oy3

# Effect of age on initial level
y0 ~ a_L * age1
y0 ~ a_NL * age1NL
"

# For the 2 step model (model 13), we just need to replace the predictors of y0
twoStepModelMG <- gsub("y0 ~ a_L [*] age1\ny0 ~ a_NL [*] age1NL",
                     "y0 ~ a_NL * hatY1", oneStepModelMG)

# Print model specifications
cat(oneStepModelMG)
cat(twoStepModelMG)




## Estimation of the Multi-group LCS models --------------------
# The model below may lead to non-positive definite matrices for some cohorts
# Generally, this is not a problem for recovering the covariance structure
fit1stepMG <- lavaan(model = oneStepModelMG, data = d,
                     group = "cohort",
                     group.equal = c("loadings", "intercepts","means",
                                     "regressions","residuals")  )

fit2stepMG <- lavaan(model = twoStepModelMG, data = d,
                     group = "cohort",
                     group.equal = c("loadings", "intercepts","means",
                                     "regressions","residuals")  ) 


## Recover THE COVARIANCE STRUCTURE ----------------
cv1 <- inspect(fit1stepMG, what = "cov.lv")
cv2 <- inspect(fit2stepMG, what = "cov.lv")


## 1. Extract covariance for y1-y3 and rename it
getCvs <- function(kCoh, cv, tLag = 2L) {
  out <- cv[[kCoh]]
  out <- out[rownames(out) %in% c("y1", "y2", "y3"),
             colnames(out) %in% c("y1", "y2", "y3")]
  colnames(out) <- rownames(out) <- paste0("y", c(kCoh, kCoh+tLag, kCoh+2*tLag)) 
  return(out)   }

cv1 <- lapply(seq_along(cv1), getCvs, cv=cv1)
cv2 <- lapply(seq_along(cv2), getCvs, cv=cv2)



## 2. Create table of moments in the age-based matrix
maxTime <- max(d$cohort)+4
timepoints <- paste0("y", seq(1, maxTime ))
combTP <- expand.grid(col=timepoints, row=timepoints)
combTP <- combTP[c("row", "col")]
combTP$name <- paste(combTP$row, combTP$col, sep="-")


# Extract one moment from one cohort
get1moment <- function(momt, cohIn, combTP) {
  cMoment <- as.matrix(combTP[momt,c("row", "col")])
  out <- cohIn[row.names(cohIn)==cMoment[[1]], colnames(cohIn)==cMoment[[2]]   ]
  out <- ifelse(length(out)!=1, NA, out) 
  return(out) }

# Extract all moments from one cohort
extractMmt <- function(cCoh, cmbTP = combTP) {
  cMoments <- sapply(seq_len(nrow(cmbTP)), get1moment, cohIn=cCoh, combTP=cmbTP ) 
  names(cMoments) <- cmbTP$name
  return(cMoments) 
}


## Obtain and average estimates from all cohorts
allMoments1 <- sapply(cv1, extractMmt, cmbTP=combTP)
allMoments2 <- sapply(cv2, extractMmt, cmbTP=combTP)
combTP$momtMeans1 <- rowMeans(allMoments1, na.rm = TRUE)
combTP$momtMeans2 <- rowMeans(allMoments2, na.rm = TRUE)
combTP$momtMeans1[is.nan(combTP$momtMeans1)] <- NA
combTP$momtMeans2[is.nan(combTP$momtMeans2)] <- NA




## 3. Reconstruct the age-based covariance matrix
expMat1 <- matrix(NA, maxTime,maxTime)
colnames(expMat1) <- rownames(expMat1) <- timepoints
expMat2 <- expMat1

for (cMomt in seq_len(nrow(combTP)) ) {
  expMat1[combTP[cMomt,"row"],combTP[cMomt,"col"]] <- combTP[cMomt, "momtMeans1"]
  expMat2[combTP[cMomt,"row"],combTP[cMomt,"col"]] <- combTP[cMomt, "momtMeans2"] 
  }

expMat1[upper.tri(expMat1)] <- NA
expMat2[upper.tri(expMat2)] <- NA



