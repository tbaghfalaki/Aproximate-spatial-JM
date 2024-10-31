# Clear the workspace by removing all objects
rm(list = ls())

# Load necessary libraries for statistical modeling, spatial analysis, and visualization
library(mvtnorm)
library(INLA)         # Integrated Nested Laplace Approximations for Bayesian inference
library(extraDistr)   # Provides additional probability distributions
library(splines)      # For spline-based smoothing
library(SpatialEpi)   # Spatial epidemiology methods
library(cplm)         # Compound Poisson GLM
library(survival)     # Survival analysis tools
library(ggplot2)      # Data visualization
library(dplyr)        # Data manipulation
library(ggfortify)    # ggplot2-based data visualization for statistical models
library(R2OpenBUGS)   # Interface for Bayesian modeling using OpenBUGS

# Set the working directory to the location of the data files
setwd("/Users/taban/Desktop/Taban_areal/")

# Initialize empty lists to store joint, longitudinal, and survival data
joint <- list()
longdat <- list()
survdat <- list()

# Load the data from a text file into a data frame
Data <- read.table("Data_new.txt", header = TRUE)

# Display column names and first few rows of the data for inspection
names(Data)
attach(Data)
head(Data)

# Additional library loaded for mixed model correlation functions
library(CorrMixed)

# Data preprocessing: Convert 'timecd4' from days to years and take square root of 'cd4' for normalization
Data$timecd41 <- timecd4 / 365.25
Data$cd4s <- sqrt(Data$cd4)

# Generate a spaghetti plot of CD4 measurements over time for each patient
Spaghetti.Plot(
  Dataset = Data, Outcome = cd4s, Id = patient, Time = timecd41,
  Add.Mean = FALSE, Add.Median = FALSE, Col = "blue",
  xlab = "Time (Year)", ylab = expression(paste(sqrt(CD4)))
)

# Define some useful variables for the longitudinal analysis
n2 <- length(unique(patient))  # Number of unique patients
Nobs <- length(patient)        # Total number of observations

# Load spatial data for Brazil, representing geographic regions or states
Fox_Lattice <- sf::st_read("Brasil.shp")
plot(Fox_Lattice)
# Initialize vectors to store demographic and clinical data by patient
ZONE_CODE <- Surv0 <- Surv1 <- Cen <- Gender <- Race <- Age <- POI <- c()
m <- c()
TIME <- CD4 <- matrix(NA, n2, 17)  # Matrices to store longitudinal data

# Loop through each unique patient to extract and organize their data
r <- 0
for (k in unique(patient)) {
  r <- r + 1
  m[r] <- length(states[patient == k])           # Number of observations for this patient
  ZONE_CODE[r] <- as.numeric(states[patient == k][1]) # State code for patient's location
  Surv0[r] <- as.numeric(dtini[patient == k][1])      # Start time of follow-up
  Surv1[r] <- as.numeric(dtend[patient == k][1])      # End time of follow-up
  Cen[r] <- as.numeric(max(Censure[patient == k]))    # Censoring indicator
  Gender[r] <- as.numeric(max(sex[patient == k]))     # Gender
  Race[r] <- as.numeric(max(race[patient == k]))      # Race
  Age[r] <- as.numeric(max(age[patient == k]))        # Age at baseline
  POI[r] <- as.numeric(max(prevoi[patient == k]))     # Previous opportunistic infection
  TIME[r, 1:length(timecd4[patient == k])] <- timecd4[patient == k] / 365.25 # Time in years
  CD4[r, 1:length(timecd4[patient == k])] <- cd4[patient == k]               # CD4 counts
}

# Summarize age variable for reference
mean(Age)
max(Age)

# Calculate survival time as the difference between start and end times
Surv <- Surv1 - Surv0
ID <- rep(c(1:n2), times = m)  # Patient IDs for each observation

# Display distribution of ZONE_CODE to understand data across different regions
table(ZONE_CODE)

# Load the 'spdep' package for spatial dependency analysis
require(spdep)

# Create a neighborhood list for spatial contiguity based on regions in the spatial object
Lattice_Temp <- poly2nb(Fox_Lattice)

# Create the adjacency matrix in INLA format, save it as "Lattice.graph"
nb2INLA("Lattice.graph", Lattice_Temp)

# Set the location of the adjacency matrix for INLA to use
Lattice.adj <- paste(getwd(), "/Lattice.graph", sep = "")

# Configure INLA options to disable scaling
inla.setOption(scale.model.default = FALSE)

# Load the adjacency matrix as a graph for spatial analysis
H <- inla.read.graph(filename = "Lattice.graph")

# Plot adjacency matrix to visualize spatial dependencies
image(inla.graph2matrix(H), xlab = "", ylab = "")

# Generate knot points for a B-spline basis on time
knots <- quantile(timecd4 / 365.25, prob = c(.25, .50, .75))
SS <- bsp(timecd4 / 365.25, k = 5)$Z  # Generate B-spline basis with 5 knots

# Create lists to store longitudinal and survival data for modeling
longdat1 <- survdat1 <- list()

# Populate longitudinal data list with CD4 count, time, and demographic data
longdat1$y <- sqrt(cd4)                    # Square root of CD4 counts
longdat1$TIME <- timecd4 / 365.25          # Time in years
longdat1$Gender <- sex
longdat1$Race <- race
longdat1$Age <- (age - mean(Age)) / sd(Age) # Age standardized
longdat1$POI <- prevoi
longdat1$ID <- ID

# Populate survival data list with censoring indicator, survival time, and demographic data
survdat1$CENSOR <- Cen
survdat1$SURVTIME <- Surv / 365.25         # Survival time in years
survdat1$Gender <- Gender
survdat1$Race <- Race
survdat1$Age <- (Age - mean(Age)) / sd(Age) # Age standardized
survdat1$POI <- POI

# Assign populated lists to main data objects for joint modeling
longdat <- longdat1
survdat <- survdat1

# Define the number of observations for longitudinal and survival datasets
n1 <- length(longdat1$y)
n2 <- length(survdat1$CENSOR)


################## ################## ################## ###############
# Prepare joint data objects for longitudinal and survival analysis
y.long <- c(longdat$y) # Extract longitudinal outcome (CD4 count)
y.surv <- inla.surv(time = c(rep(NA, n1), survdat$SURVTIME), event = c(rep(NA, n1), survdat$CENSOR)) # Define survival data for INLA
Yjoint <- list(y.long, y.surv) # Combine longitudinal and survival outcomes

# Fit a Kaplan-Meier survival curve to visualize survival data
km_trt_fit <- survfit(Surv(SURVTIME, CENSOR) ~ 1, data = survdat1)

# Extract survival curve data for ggplot
km_data <- data.frame(
  time = km_trt_fit$time,
  surv = km_trt_fit$surv,
  lower = km_trt_fit$lower,
  upper = km_trt_fit$upper
)

# Plot the Kaplan-Meier curve with ggplot, showing survival probability over time
ggplot(km_data, aes(x = time, y = surv)) +
  geom_step(color = "blue") + # Survival curve in blue
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) + # Shadow effect for confidence interval
  ggtitle("Kaplan-Meier Survival Curve") +
  xlab("Time (Year)") +
  ylab("Survival Probability")

# Define linear covariates for joint model, including both longitudinal and survival predictors
linear.covariate <- data.frame(
  mu = as.factor(c(rep(1, n1), rep(2, n2))), # Indicator for model component (1=longitudinal, 2=survival)
  l.TIME = c(longdat$TIME, rep(0, n2)),      # Longitudinal time variable
  l.Age = c(longdat$Age, rep(0, n2)),        # Age for longitudinal part
  l.Gender = c(longdat$Gender, rep(0, n2)),  # Gender for longitudinal part
  l.POI = c(longdat$POI, rep(0, n2)),        # Previous opportunistic infection for longitudinal part
  s.Age = c(rep(0, n1), survdat$Age),        # Age for survival part
  s.Gender = c(rep(0, n1), survdat$Gender),  # Gender for survival part
  s.POI = c(rep(0, n1), survdat$POI)         # Previous opportunistic infection for survival part
)

# Extend spline basis matrix for longitudinal model
SS2 <- matrix(0, n2, length(SS[1, ]))
linear.covariate$l.SS <- rbind(SS, SS2) # Combine longitudinal splines and zero-padded splines for survival

# Define random covariates for joint model
random.covariate <- list(
  U11 = c(longdat1$ID, rep(NA, n2)),           # Random intercept for longitudinal
  U21 = c(longdat1$ID + n2, rep(NA, n2)),      # Random slope for longitudinal
  U12 = c(rep(NA, n1), 1:n2),                  # Random intercept for survival
  U22 = c(rep(NA, n1), n2 + (1:n2)),           # Random slope for survival
  ZONE_CODE = c(rep(NA, n1), ZONE_CODE)        # Spatial random effect (region-based)
)

# Combine linear and random covariates, and outcome data for INLA modeling
joint.data <- c(linear.covariate, random.covariate)
joint.data$Y <- Yjoint

# Define the joint model formula with fixed effects, random effects, and spatial effects
## b0+b1t, gamma1*b0+gamma2*b1, v_k
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11, model = "iid2d", n = 2 * n2) +       # Random intercept for longitudinal component
  f(U21, l.TIME, copy = "U11") +              # Random slope for longitudinal, copying structure of U11
  f(U12, copy = "U11", fixed = FALSE) +       # Random intercept for survival, unlinked from longitudinal
  f(U22, copy = "U11", fixed = FALSE) +       # Random slope for survival, unlinked from longitudinal
  f(ZONE_CODE, model = "besag", graph = Lattice.adj, # Spatial effect using Besag model
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.1)))) # Penalized complexity prior

# Record start time to measure model execution duration
start_time0 <- proc.time()

# Fit joint model using INLA with Gaussian distribution for longitudinal data and Weibull for survival data
joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"),
                   data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE), # Compute model evaluation metrics
                   control.family = list(list(), list()) # Separate control parameters for each family
)

# Print summary of the INLA model results
summary(joint.inla)

print("End_INLA")

# Calculate execution time in minutes
end_time0 <- proc.time()
execution_time0 <- end_time0 - start_time0
elapsed_time_in_minutes0 <- execution_time0["elapsed"] / 60

# Extract fixed effects estimates
joint.inla$summary.fixed[, 1] # Point estimates of fixed effects
joint.inla$summary.fixed[, c(3, 5)] # 95% credible intervals for fixed effects

# Extract hyperparameters (e.g., random effects variance)
joint.inla$summary.hyperpar[, 1] # Point estimates for hyperparameters
joint.inla$summary.hyperpar[, c(3, 5)] # 95% credible intervals for hyperparameters

# Compute the log pseudo-marginal likelihood (LPML) for model comparison
LPML1 <- sum(log(joint.inla$cpo$cpo))

# Output model quality metrics for model evaluation
joint.inla$waic$waic # WAIC (Widely Applicable Information Criterion)
joint.inla$dic$dic   # DIC (Deviance Information Criterion)



#################################################
# Full model with spatial and temporal random effects
## b0+b1t, 0, v_k
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11, model = "iid2d", n = 2 * n2) +                # Random intercept
  f(U21, l.TIME, copy = "U11") +                       # Random slope for time linked to U11
  f(ZONE_CODE, model = "besag", graph = Lattice.adj,   # Spatial effect with Besag model
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.1))))

# Fit the model using INLA
joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"), data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list()))
summary(joint.inla)

# Calculate LPML
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

# Retrieve WAIC and DIC values
joint.inla$waic$waic
joint.inla$dic$dic

#################################################
# Model without spatial effect (temporal random effects only)
## b0+b1t, 0, 0
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11, model = "iid2d", n = 2 * n2) + 
  f(U21, l.TIME, copy = "U11")

joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"), data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list()))
summary(joint.inla)

# Calculate LPML, WAIC, and DIC
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
joint.inla$waic$waic
joint.inla$dic$dic

################################################# 
# Fixed-effects-only model (no random effects)
## 0, 0, 0
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI

joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"), data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list()))
summary(joint.inla)

# Calculate LPML, WAIC, and DIC
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
joint.inla$waic$waic
joint.inla$dic$dic

################################################# 
# Model with individual random intercept only
## b0, 0, 0
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11)

joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"), data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list()))
summary(joint.inla)

# Calculate LPML, WAIC, and DIC
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
joint.inla$waic$waic
joint.inla$dic$dic

################################################# 
# Model with additional random effect for survival component (no spatial effect)
## b0, gamma1*b0, 0
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11) + 
  f(U12, copy = "U11", fixed = FALSE)

joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"), data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list()))
summary(joint.inla)

# Calculate LPML, WAIC, and DIC
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
joint.inla$waic$waic
joint.inla$dic$dic

################################################# 
# Full model with random effects and spatial dependency
## b0, gamma1*b0, v_k
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11) + 
  f(U12, copy = "U11", fixed = FALSE) + 
  f(ZONE_CODE, model = "besag", graph = Lattice.adj,
    hyper = list(prec = list(prior = "pc.prec", param = c(0.1, 0.1))))

joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"), data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list()))
summary(joint.inla)

# Calculate LPML, WAIC, and DIC
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
joint.inla$waic$waic
joint.inla$dic$dic

#################################################
# Base model with Spatial correlation based on ZONE_CODE using a Besag model
## 0, 0, v_k
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(ZONE_CODE, model = "besag", graph = Lattice.adj, hyper = list(prec = list(prior = "pc.prec", param = c(01, 0.1))))

# Fit the model using INLA
joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"),
                   data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list())
)

# Summarize the results
summary(joint.inla)

# Calculate Log Pseudo Marginal Likelihood (LPML)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

# Extract WAIC and DIC for model comparison
joint.inla$waic$waic
joint.inla$dic$dic

#################################################
# Add a random intercept (U11) and a time-dependent random slope (U21)
# Random intercept allows for individual variation, while U21 captures time effects
## b0+b1t, gamma1*b0, v_k
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11, model = "iid2d", n = 2 * n2) +       # Random intercept
  f(U21, l.TIME, copy = "U11") +              # Time-varying random slope
  f(U12, copy = "U11", fixed = FALSE) +       # Another random slope copy of U11
  f(ZONE_CODE, model = "besag", graph = Lattice.adj, hyper = list(prec = list(prior = "pc.prec", param = c(01, 0.1))))

# Fit the model using INLA
joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"),
                   data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list())
)

# Summarize the results
summary(joint.inla)

# Calculate LPML
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

# Extract WAIC and DIC
joint.inla$waic$waic
joint.inla$dic$dic

#################################################
# Add another independent random slope (U22)
# This allows for additional variability in the slopes
## b0+b1t, gamma2*b1, v_k
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11, model = "iid2d", n = 2 * n2) +
  f(U21, l.TIME, copy = "U11") +
  f(U22, copy = "U11", fixed = FALSE) +       # Additional random slope
  f(ZONE_CODE, model = "besag", graph = Lattice.adj, hyper = list(prec = list(prior = "pc.prec", param = c(01, 0.1))))

# Fit the model using INLA
joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"),
                   data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list())
)

# Summarize the results
summary(joint.inla)

# Calculate LPML
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

# Extract WAIC and DIC
joint.inla$waic$waic
joint.inla$dic$dic

#################################################
# Include additional random effects but without the spatial structure
# Focus on the impact of both random slopes and intercepts
## b0+b1t, gamma1*b0+gamma2*b1, 0
formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11, model = "iid2d", n = 2 * n2) +
  f(U21, l.TIME, copy = "U11") +
  f(U12, copy = "U11", fixed = FALSE) +       # Random slope
  f(U22, copy = "U11", fixed = FALSE)          # Additional random slope

# Fit the model using INLA
joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"),
                   data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list())
)

# Summarize the results
summary(joint.inla)

# Calculate LPML
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

# Extract WAIC and DIC
joint.inla$waic$waic
joint.inla$dic$dic

#################################################
# Final model: Add complexity by combining all elements
# This model captures both random slopes, intercepts, and spatial effects
## b0+b1t, gamma1*b0+gamma2*b1, v_k

formula <- Y ~ mu - 1 + l.Age + l.Gender + l.POI + l.SS + s.Age + s.Gender + s.POI +
  f(U11, model = "iid2d", n = 2 * n2) +
  f(U21, l.TIME, copy = "U11") +
  f(U12, copy = "U11", fixed = FALSE) +
  f(U22, copy = "U11", fixed = FALSE) +
  f(ZONE_CODE, model = "besag", graph = Lattice.adj, hyper = list(prec = list(prior = "pc.prec", param = c(01, 0.1))))

# Fit the model using INLA
joint.inla <- inla(formula,
                   family = c("gaussian", "weibullsurv"),
                   data = joint.data,
                   control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                   control.family = list(list(), list())
)

# Summarize the results
summary(joint.inla)
A=joint.inla$summary.hyperpar[,c(1,2,3,5,4)]
xtable::xtable(A,digits = 3)

# Calculate LPML
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

# Extract WAIC and DIC
joint.inla$waic$waic
joint.inla$dic$dic


# Load necessary libraries for mapping and visualization
library(sf)        # For handling spatial data
library(ggplot2)   # For creating plots
library(viridis)   # For color scales
library(gridExtra) # For arranging multiple plots


####################### Map #########################

# Read spatial data (shapefile for Brazil)
nc <- st_read("Brasil.shp", quiet = TRUE)

# Prepare to plot multiple figures in a single plot area
par(mfrow = c(1, 2))

# Alternative path for the shapefile (if necessary)
# nc <- st_read("//Users//taban//Desktop//Rui_Martins//Brasil.shp", quiet = TRUE)

# Create a vector v0 containing predicted values for each zone
v0 <- exp(joint.inla$summary.random$ZONE_CODE[, 2]) # Exponentiating to get the predicted values

# Plot the first map with predicted values by ZONE_CODE
p1 <- ggplot(data = nc, aes(fill = v0)) +  # Use ggplot2 for spatial plotting
  geom_sf() +                             # Add spatial features
  scale_fill_viridis() +                  # Use viridis color scale for better visualization
  theme_bw() +                            # Use a clean black-and-white theme
  theme(legend.title = element_blank())    # Remove legend title for clarity

###################################
# Aggregate predictions for the linear predictor
v <- aggregate(joint.inla$summary.linear.predictor[(n1 + 1):(n1 + n2), 1] ~ ZONE_CODE, FUN = mean)[, 2]

# Create a new vector for visual representation
v1 <- rep(1, 27)                        # Initialize vector with default values
v1[1:12] <- v[1:12]                     # Assign predictions for the first 12 zones
v1[15:27] <- v[13:25]                   # Assign predictions for zones 13-25
v1[13:14] <- mean(v1[-c(13:14)])       # Calculate mean for the zones 13 and 14

# Plot the second map with aggregated predictions
p2 <- ggplot(data = nc, aes(fill = exp(v1))) +  # Exponentiating for visual scale
  geom_sf() +                                    # Add spatial features
  scale_fill_viridis() +                         # Use viridis color scale
  theme_bw() +                                   # Use a clean theme
  theme(legend.title = element_blank())           # Remove legend title for clarity

# Arrange both plots side by side for comparison
grid.arrange(p2, p1, ncol = 2)


#############

#############################
# Extracting hyperparameter estimates from the joint INLA model
r1 = joint.inla$summary.hyperpar$mean[7]  # Extract the first random effect hyperparameter
r2 = joint.inla$summary.hyperpar$mean[8]  # Extract the second random effect hyperparameter

# Constructing the covariance matrix for the random effects
D = matrix(c(joint.inla$summary.hyperpar$mean[2]^-1,
             joint.inla$summary.hyperpar$mean[2]^-1 * joint.inla$summary.hyperpar$mean[3]^-1 * joint.inla$summary.hyperpar$mean[4],
             joint.inla$summary.hyperpar$mean[2]^-1 * joint.inla$summary.hyperpar$mean[3]^-1 * joint.inla$summary.hyperpar$mean[4],
             joint.inla$summary.hyperpar$mean[3]^-1), 2, 2)

# Setting a weight for combining the identity matrix with D
a = 1
D1 = a * diag(1, 2) + (1 - a) * D  # Adjusted covariance matrix

# Generating random variates from a multivariate normal distribution
b = rmvnorm(n2, rep(0, 2), D1)  # Simulating random effects

# Extracting another hyperparameter for variance
tau = joint.inla$summary.hyperpar$mean[5]  # Extract the variance hyperparameter
Q = matrix(tau, H$n, H$n)  # Initialize a covariance matrix based on tau
for (kk in 1:H$n) {  # Loop through each element to adjust the diagonal
  Q[kk, kk] = tau * (H$nnbs[kk])  # Modify diagonal based on neighborhood structure
}

# Regularizing the covariance matrix
a = 0.01
Q1 = a * diag(1, H$n) + (1 - a) * Q  # Regularized covariance matrix
u = c()  # Initialize a vector for random effects
u = rmvnorm(1, rep(0, H$n), solve(Q1))  # Simulating random effects using the regularized covariance matrix

# Calculating the hazard function (lambda) for each patient
lambda = c()  # Initialize a vector for hazard rates
for (patientnr in 1:n2) { 
  # Constructing the covariate vector for each patient
  X = c(1, linear.covariate$s.Age[n1 + patientnr], linear.covariate$s.Gender[n1 + patientnr], linear.covariate$s.POI[n1 + patientnr])
  llll = length(joint.inla$summary.fixed$mean)  # Length of fixed effects estimates
  # Calculate hazard rate using fixed effects, random effects, and zone-specific effects
  lambda[patientnr] = c(joint.inla$summary.fixed$mean[2], joint.inla$summary.fixed$mean[(llll - 2):llll]) %*% X +
    r1 * b[patientnr, 1] + r2 * b[patientnr, 2] + u[ZONE_CODE[patientnr]]
}

# Plotting setup for survival analysis
par(mfrow = c(1, 1))  # Reset plotting layout

# Calculating survival probabilities using the Weibull distribution
surv2 = 1 - pweibull(survdat1$SURVTIME, r, exp(-lambda / r))  # Compute survival probabilities
# res.coxsnell represents the Cox-Snell residuals
res.coxsnell <- -log(surv2)  # Compute Cox-Snell residuals

# Fitting a Kaplan-Meier survival curve to the Cox-Snell residuals
km <- survfit(Surv(res.coxsnell, survdat1$CENSOR, type = 'right') ~ 1)
# Plotting the Kaplan-Meier survival function
plot(km, xlab = "Cox-Snell Residuals", xlim = c(0, 1.2), main = "Survival Function of Cox-Snell Residuals",
     mark.time = F, conf.int = F, lwd = 2)  # Setting plot parameters

# Plotting a reference line for the survival function of an exponential distribution
x <- seq(0, 10, 0.01)  # Generating a sequence for x-axis
lines(x, 1 - pexp(x), "l", lwd = 1, lty = 2)  # Adding exponential survival function for comparison






########
# Extract predicted values from the model
predl = joint.inla$summary.fitted.values$mean[1:n1]
preds = joint.inla$summary.fitted.values$mean[(1+n1):(n2+n1)]

# Initialize a matrix to store predictions for each patient
predL = matrix(NA, n2, 17)

# Populate the prediction matrix with values based on patient IDs
r = 0
for (k in unique(patient)) {
  r = r + 1
  predL[r, 1:length(timecd4[patient == k])] = predl[patient == k]
}

# Define specific patients for plotting
patients = c(12, 230, 489)

# Set up the custom plotting layout: 3 rows and 2 columns with the first column wider
layout(matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE), widths = c(2, 1)) 

# Adjust margins for the plots to reduce spacing
par(mar = c(2, 4, 2, 1))  # Bottom, Left, Top, Right margins

# Extract hyperparameter estimates
alpha = joint.inla$summary.hyperpar$mean[1]
j_est = joint.inla$summary.hyperpar$mean[2:6]
f_est = joint.inla$summary.fixed$mean

# Loop through each specified patient to generate plots
for (i in 1:length(patients)) {
  patientnr = patients[i]
  dataHi = Data[Data$patient == patientnr, ]
  
  # Plot CD4 trajectory for the current patient
  plot(TIME[patientnr, ], sqrt(CD4[patientnr, ]),
       ylab = "CD4 count", 
       xlab = "Time (months)", 
       type = "l",
       xlim = c(0, 5), 
       ylim = c(0, 45),
       main = paste("CD4 trajectory - patient", patientnr))
  
  # Add predictions to the plot
  lines(TIME[patientnr, ], predL[patientnr, ], col = "blue", lty = 2)
  
  # Add legend only for the first patient, without box
  if (i == 1) {
    legend("topright", c("Observed", "Prediction"), lty = 1:2, col = c(1, "blue"), bty = "n")
  }
  
  # Calculate lambda for the survival curve
  r = joint.inla$summary.hyperpar$mean[2]
  X = c(1, linear.covariate$s.Age[n1 + patientnr], 
        linear.covariate$s.Gender[n1 + patientnr], 
        linear.covariate$s.POI[n1 + patientnr])
  
  llll = length(joint.inla$summary.fixed$mean)
  lambda = c(joint.inla$summary.fixed$mean[2], 
             joint.inla$summary.fixed$mean[(llll - 2):llll]) %*% X
  
  # Generate the survival curve
  pred_time = seq(0, 5, by = 0.01)
  bbbb = (as.numeric(pred_time) / as.numeric(-lambda))^r
  
  # Plot survival probability for the current patient
  plot(pred_time, exp(-bbbb),
       type = "l", 
       ylab = "Survival probability", 
       xlab = "Time (months)",
       main = paste("Survival curve - patient", patientnr))
  
  # Add a horizontal line at y = 0.5
  abline(h = 0.5, col = "red")
}



