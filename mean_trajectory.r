# R code for the statistical models used in:

# Estrada, Hamagami, & Ferrer (2019, in press).
# Estimating age-based developmental trajectories using latent change score models
# based on measurement occasion.
# Multivariate Behavioral ReseachPsychological Methods.


## RECOVERING THE AGE-BASED MEAN TRAJECTORY

# Load required libraries
if (!"package:minpack.lm" %in% search()) {require(minpack.lm)}
if (!"package:lavaan" %in% search()) {require(lavaan)}

# Load data set
d <- myData
# Data frame in wide format. Variables of interest:
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



## Specification of LCS-SEM for recovering THE MEAN TRAJECTORY ------------------
# Specification of models
oneStepModel = "
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

y0 ~~ y0V   * y0
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
twoStepModel <- gsub("y0 ~ a_L [*] age1\ny0 ~ a_NL [*] age1NL",
                     "y0 ~ a_NL * hatY1", oneStepModel)

# Print model specifications
cat(oneStepModel)
cat(twoStepModel)



## Estimation of the LCS models for recovering THE MEAN TRAJECTORY -----------------------
fit1step <- lavaan(model = oneStepModel, data = d)
fit2steps <- lavaan(model = twoStepModel, data = d)

summary(fit1step)
summary(fit2steps)



## Computation of model expectations FOR THE MEAN TRAJECTORY -----------
# Given a particular age at t1
params1step <- parameterEstimates(fit1step)$est
names(params1step) <- parameterEstimates(fit1step)$label

params2steps <- parameterEstimates(fit2steps)$est
names(params2steps) <- parameterEstimates(fit2steps)$label


# Age at first measurement occasion
ageT1 <- 1 # (the value 0 indicates 5 years of age)


# Expectation from the 1 step model
y1 <- params1step["y0mn"] + params1step["a_L"]*ageT1 + params1step["a_NL"]*ageT1
y2 <- y1 + params1step["yAmn"] + params1step["b_y"] * y1
y3 <- y2 + params1step["yAmn"] + params1step["b_y"] * y2

oneStepExpects <- c(ageT1, y1, y2, y3)
names(oneStepExpects) <- c("age1", "y1", "y2", "y3")


# Expectation from the 2 step model
hatY1 <- predict(model13, newdata=data.frame(age1 = ageT1) )

y1 <- params2steps["y0mn"] + params2steps["a_NL"]*hatY1
y2 <- y1 + params2steps["yAmn"] + params2steps["b_y"] * y1
y3 <- y2 + params2steps["yAmn"] + params2steps["b_y"] * y2

twoStepExpects <- c(ageT1, y1, y2, y3)
names(twoStepExpects) <- c("age1", "y1", "y2", "y3")
rm(ageT1, y1, y2, y3, hatY1)

rbind(oneStepExpects,twoStepExpects)


