##### Marine conservation financial analysis
##### Start date of analysis: 9 April 2025
##### Author: Ethan Ignoffo, ethanignoffo@berkeley.edu

# Creates summary statistics for the funding allocation of ESA-listed marine species
# Cleans Marine dataset for use in analysis
# Analyzes the relationships between federal funding and population trends for ESA-listed marine species
# Analyzes the relationships between length of listing and federal funding allotment for ESA-listed marine species

##### Load libraries 
library(readxl)
library(sf)
library(tidyverse)
library(dplyr)
library(psych)
library(bestglm)
library(ResourceSelection)
library(caret)
library(lmtest)
library(ggplot2)
library(ggthemes)
library(scales)
library(wesanderson)
library(hrbrthemes)
library(stargazer)
library(htmltools)
library(webshot)

rm(list = ls())

##### Read in the data 
marine <- read_excel("/Users/ethanignoffo/Desktop/Thesis/Data/Marine.xlsx", sheet = 'DataSet')
str(marine)

##### Clean the data 
# Remove unnecessary columns 
marine <- marine %>% select(-c(RecoveryPlanYear, NSShortPop, IUCNExtantRange))
# Convert columns to numeric
marine$FedFundingTotal <- as.numeric(marine$FedFundingTotal) 
marine$YearsSinceListing <- as.numeric(marine$YearsSinceListing)
marine$IUCNThreats <- as.numeric(marine$IUCNThreats)
# Convert columns to binary 
marine$IUCNSkewPop <- ifelse(marine$IUCNSkewPop == "Increasing", 1, 0)
marine$ListingStatus <- ifelse(marine$ListingStatus == "Endangered", 1, 0)
# Standardize numeric columns 
numeric_cols <- c("FedFundingTotal", "YearsFunded", "ListingYear",
                  "YearsSinceListing", "NOAARegions", "IUCNThreats")
marine_scaled <- marine
marine_scaled[numeric_cols] <- scale(marine_scaled[numeric_cols])

############### Summary Statistics

##### 1. Test Pearson's correlation between continuous co-variates
marine_corr <- marine_scaled
marine_corr <- marine_corr %>% select(c("FedFundingTotal", "YearsFunded", 
                                        "ListingYear", "YearsSinceListing", 
                                        "IUCNThreats", "NOAARegions")) # Select continuous/numeric variables
marine_corr <- marine_corr %>% drop_na() # Drop N/As within numeric subset

# Produce manual correlation matrix & test for correlation
pairs(marine_corr, pch = 19)

# Add correlation coefficients to plot 
panel.cor <- function(x,y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0,1,0,1))
  r <- round(cor(x,y, method = "pearson"), digits = 2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(marine_corr, lower.panel = panel.cor)

# Check correlation ellipses 
pairs.panels(marine_corr,
             method = "pearson", # correlation method
             hist.col = wes_palette("FantasticFox1")[3],
             density = TRUE, # show density plots
             ellipses = TRUE # show ellipses
)

## Problematic correlations identified between YearsFunded, ListingYear, and YearsSinceListing. 
## All other correlations identified as not problematic with |r| < 0.7

##### 2. Test Point Biserial correlations for correlations between continuous and binary variables
# Conduct Point Biserial correlation tests for FedFundingTotal
cor.test(marine_scaled$RecoveryPlanBinary, marine_scaled$FedFundingTotal) # There is no correlation between federal funding and recovery plans (p>0.05)
cor.test(marine_scaled$ListingStatus, marine_scaled$FedFundingTotal) # There is no correlation between federal funding and listing status (p>0.05)
cor.test(marine_scaled$IUCNSkewPop, marine_scaled$FedFundingTotal) # There is no correlation between federal funding and population trend (p>0.05)
cor.test(marine_scaled$LandUse, marine_scaled$FedFundingTotal) # There is no correlation between federal funding and land use (p>0.05)

# Conduct Point Biserial correlation tests for YearsFunded
cor.test(marine_scaled$RecoveryPlanBinary, marine_scaled$YearsFunded) # There is no correlation between years funded and recovery plans (p>0.05) 
cor.test(marine_scaled$ListingStatus, marine_scaled$YearsFunded) # There is a moderate positive correlation between years funded and listing status (p=3.741e-05)
cor.test(marine_scaled$IUCNSkewPop, marine_scaled$YearsFunded) # There is a moderate positive correlation between years funded and population trends (p=6.791e-06)
cor.test(marine_scaled$LandUse, marine_scaled$YearsFunded) # There is a weak positive correlation between years funded and land use (p=0.008704)

# Conduct Point Biserial correlation tests for ListingYear
cor.test(marine_scaled$RecoveryPlanBinary, marine_scaled$ListingYear) # There is no correlation between listing year and recovery plans (p>0.05)
cor.test(marine_scaled$ListingStatus, marine_scaled$ListingYear) # There is a moderate negative correlation between listing year and listing status (p=0.0002296)
cor.test(marine_scaled$IUCNSkewPop, marine_scaled$ListingYear) # There is a moderate negative correlation between listing year and population trends (p=8.247e-06)
cor.test(marine_scaled$LandUse, marine_scaled$ListingYear) # There is a weak negative correlation between listing year and land use (p=0.01685)

# Conduct Point Biserial correlation tests for YearsSinceListing
cor.test(marine_scaled$RecoveryPlanBinary, marine_scaled$YearsSinceListing) # There is no correlation between the years since listing and recovery plans (p>0.05)
cor.test(marine_scaled$ListingStatus, marine_scaled$YearsSinceListing) # There is a moderate positive correlation between the years since listing and listing status (p=0.0004297)
cor.test(marine_scaled$IUCNSkewPop, marine_scaled$YearsSinceListing) # There is a moderate positive correlation between the years since listing and population trends (p=1.055e-05)
cor.test(marine_scaled$LandUse, marine_scaled$YearsSinceListing) # There is a weak positive correlation between the years since listing and land use (p=0.02132)

# Conduct Point Biserial correlation tests for IUCNThreats
cor.test(marine_scaled$RecoveryPlanBinary, marine_scaled$IUCNThreats) # There is no correlation between IUCN threats and recovery plans (p>0.05)
cor.test(marine_scaled$ListingStatus, marine_scaled$IUCNThreats) # There is a weak negative correlation between IUCN threats and listing status (p=0.03196)
cor.test(marine_scaled$IUCNSkewPop, marine_scaled$IUCNThreats) # There is no correlation between IUCN threats and population trends (p>0.05)
cor.test(marine_scaled$LandUse, marine_scaled$IUCNThreats) # There is no correlation between IUCN threats and land use (p>0.05)

# Conduct Point Biserial correlation tests for NOAARegions
cor.test(marine_scaled$RecoveryPlanBinary, marine_scaled$NOAARegions) # There is no correlation between NOAA regions and recovery plans (p>0.05)
cor.test(marine_scaled$ListingStatus, marine_scaled$NOAARegions) # There is a weak positive correlation between NOAA regions and listing status (p=0.003139)
cor.test(marine_scaled$IUCNSkewPop, marine_scaled$NOAARegions) # There is a weak positive NOAA regions and population trends (p=0.005021)
cor.test(marine_scaled$LandUse, marine_scaled$NOAARegions) # There is no correlation between NOAA regions and land use (p>0.05)

##### 3. Test Chi-2 correlations for categorical non-binary variables
marine_corr2 <- marine_scaled %>% select(c("ScientificName", "CommonName", "SpeciesType")) # Select categorical non-binary variables

# Conduct for Chi-2 tests for all categorical non-binary variables
# Conduct for all variables 
chisq.test(marine_corr2$ScientificName, marine_corr2$CommonName) # There is no correlation between common name and scientific name (p>0.05)
chisq.test(marine_corr2$ScientificName, marine_corr2$SpeciesType) # There is no correlation between species type and scientific name (p>0.05)
chisq.test(marine_corr2$CommonName, marine_corr2$SpeciesType) # There is no correlation between species type and common name (p>0.05)

## There is no evidence of significant correlations between the categorical non-binary variables in this analysis

############### ANOVA Exploratory Analysis

# Define dataset for ANOVA
marine_aov <- na.omit(marine)

# Convert column back to categorical for ANOVA 
marine_aov$ListingStatus <- ifelse(marine_aov$ListingStatus == 1, "Endangered", "Threatened")

# ANOVA 1: Species type on federal funding
a1 <- aov(FedFundingTotal ~ SpeciesType, data = marine_aov)
summary(a1)  # Fail to reject null (p>0.05)

# ANOVA 2: Listing status on federal funding
a2 <- aov(FedFundingTotal ~ ListingStatus, data = marine_aov)
summary(a2) # Fail to reject null (p>0.05)

# ANOVA 3: Land use on federal funding
a3 <- aov(FedFundingTotal ~ LandUse, data = marine_aov)
summary(a3) # Fail to reject null (p>0.05)

############### Population trend analysis (Logistic regression)

# Define dataset for regression
marine_log <- marine_scaled %>% select(c("IUCNSkewPop", "FedFundingTotal", "YearsFunded", "ListingYear",
                                         "YearsSinceListing", "RecoveryPlanBinary", "ListingStatus",
                                         "IUCNThreats", "NOAARegions", "LandUse"))
# Omit N/As 
marine_log <- na.omit(marine_log)

##### 1. Fit the models
# Preliminary logistic regression: all covariates included
glm.base <- glm(IUCNSkewPop ~ FedFundingTotal + YearsFunded + 
              ListingStatus + RecoveryPlanBinary + IUCNThreats + 
              NOAARegions + LandUse,
            family = binomial,
            data = marine_log)
summary(glm.base)

# Backwards step-wise selection: removed insignificant covariates manually
glm.backwards <- glm(IUCNSkewPop ~ YearsFunded + ListingStatus,
                     family = binomial, 
                     data = marine_log)
summary(glm.backwards)

# Compare AICs from preliminary and backwards stepwise 
AIC(glm.base, glm.backwards) # Base = 42.78300, Backwards = 39.01822

# Predict outcome based on model coefficients
glm.probs <- predict(glm.base, type = "response")
glm.probs[1:10]
glm.pred <- rep(0, length(glm.probs))
glm.pred[glm.probs > 0.5] = 1
table(glm.pred, marine_log$IUCNSkewPop)

# Best subset selection 
marine_log <- as.data.frame(marine_log) # Convert to true dataframe
best.logit <- bestglm(marine_log,
                      IC = "AIC", # Information criteria for
                      family = binomial,
                      method = "exhaustive")
summary(best.logit$BestModel) # AIC

# Best model is as follows
glm.fits <- glm(IUCNSkewPop ~ FedFundingTotal + YearsFunded,
                family = binomial,
                data = marine_log)
summary(glm.fits) # Null deviance = 56.70, Residual deviance = 34.78

# Compare AICs from all 3 
AIC(glm.base, glm.backwards, glm.fits) # Base = 42.78300, Backwards = 39.01822, Best = 40.78036

# The model of best fit reduces the AIC value compared to the base model, while still including FedFundingTotal in the glm

##### 2. Test the models
# Hoslem tests 
hl_best <- hoslem.test(glm.fits$y, fitted(glm.fits), g = 10) # Hosmer-and-Lemeshow goodness of fit (GOF) test
hl_best # No evidence of poor fit (p=0.4742)

# Confusion matrix 
confmtx <- confusionMatrix(data = factor(glm.pred), reference = factor(marine_log$IUCNSkewPop))
confmtx # True positives and true negatives outnumber false positives and false negatives

############### Total federal funding analysis (Log-linear regression)

# Define dataset for regression: utilizing unscaled data for this regression
marine_lm <- marine %>% select(c("YearsFunded", "ListingStatus", 
                                 "RecoveryPlanBinary", "IUCNThreats",
                                 "NOAARegions", "LandUse", "FedFundingTotal"))

# Omit N/As
marine_lm <- na.omit(marine_lm)

##### 1. Fit the models
# Preliminary linear regression: all covariates included
lm.base <- lm(FedFundingTotal ~ YearsFunded + ListingStatus + 
                RecoveryPlanBinary + IUCNThreats + NOAARegions + LandUse,
              data = marine_lm)

summary(lm.base) # Estimations are extremely large, will continue with log transformation of dependent variable 

# Log transform dependent variable 
marine_lm$FedFundingTotal <- log(marine_lm$FedFundingTotal)

# Preliminary linear regression: dependent variable log transformed
lm.base.log <- lm(FedFundingTotal ~ YearsFunded + ListingStatus + 
                RecoveryPlanBinary + IUCNThreats + NOAARegions + LandUse,
              data = marine_lm)

summary(lm.base.log) # Estimations are far more reasonable after log transformation

# Best subset selection 
marine_lm <- as.data.frame(marine_lm) # Convert to true dataframe
best.logit2 <- bestglm(marine_lm,
                       family = gaussian, # Using gaussian as response variable has been log transformed
                      IC = "AIC", # Information criteria for
                      method = "exhaustive")
summary(best.logit2$BestModel)

# Best model is as follows
lm.fits <- lm(FedFundingTotal ~ YearsFunded + RecoveryPlanBinary, 
              data = marine_lm)
summary(lm.fits)

# Pull AIC from models
AIC(lm.base, lm.base.log, lm.fits) # Base = 2332.3832, Base Log = 188.1454, Best = 186.5667

##### 2. Diagnostic Plots
plot(lm.fits) # Print out diagnostic plots linear model of best fit

############### Generate Plots

##### 1. ANOVA Plots
# Sample: Species Tyoe
sample1 = marine_aov %>% group_by(SpeciesType) %>% summarize(num=n())


# Plot: FedFundingTotal by Species Type
marine_aov %>%
  left_join(sample1) %>%
  mutate(myaxis = paste0(SpeciesType, "\n", "n=", num)) %>%
  ggplot(aes(x = myaxis, 
             y = FedFundingTotal, 
             fill = SpeciesType)) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.1, 
               color = "black", 
               alpha = 0.2) +
  scale_fill_manual(
    values = wes_palette("FantasticFox1", 
                         n = length(unique(marine_aov$SpeciesType)), 
                         type = "continuous")) +
  theme_light() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  ggtitle("Distribution of Federal Funding by Species Type") +
  xlab("Species Type") +
  ylab("Total Federal Funding") +
  scale_y_log10(labels = label_dollar(scale = 1e-6, suffix = "M"))

# Sample: Listing Status
sample2 = marine_aov %>% group_by(ListingStatus) %>% summarize(num=n())

# Plot: FedFundingTotal by ListingStatus
marine_aov %>%
  left_join(sample2) %>%
  mutate(myaxis = paste0(ListingStatus, "\n", "n=", num)) %>%
  ggplot(aes(x = myaxis, 
             y = FedFundingTotal, 
             fill = factor(ListingStatus))) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.1, 
               color = "black", 
               alpha = 0.2) +
  scale_fill_manual(
    values = wes_palette("FantasticFox1",
                         n = length(unique(marine_aov$ListingStatus)), 
                         type = "continuous")) +
  theme_light() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)) +
  ggtitle("Distribution of Federal Funding by Listing Status") +
  xlab("Listing Status") +
  ylab("Total Federal Funding") +
  scale_y_log10(labels = label_dollar(scale = 1e-6, suffix = "M"))

# Sample: Land Use
sample3 = marine_aov %>% group_by(LandUse) %>% summarize(num=n())

# Plot: FedFundingTotal by LandUse
marine_aov %>%
  left_join(sample3) %>%
  mutate(myaxis = paste0(LandUse, "\n", "n=", num)) %>%
  ggplot(aes(x = myaxis, 
             y = FedFundingTotal, 
             fill = factor(LandUse))) +
  geom_violin(width = 1.4) +
  geom_boxplot(width = 0.1, 
               color = "black", 
               alpha = 0.2) +
  scale_fill_manual(
    values = wes_palette("FantasticFox1",
                         n = length(unique(marine_aov$LandUse)), 
                         type = "discrete")) +
  theme_light() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)) +
  ggtitle("Distribution of Federal Funding by Land Use") +
  xlab("Land Use") +
  ylab("Total Federal Funding") +
  scale_y_log10(labels = label_dollar(scale = 1e-6, suffix = "M"))

##### 1. Logistic Plots
# Define un-scaled data set for plot 
marine_log2 <- marine %>% select(c("IUCNSkewPop", "FedFundingTotal", "YearsFunded", "ListingYear",
                                         "YearsSinceListing", "RecoveryPlanBinary", "ListingStatus",
                                         "IUCNThreats", "NOAARegions", "LandUse"))
# Omit N/As 
marine_log2 <- na.omit(marine_log2)

# Fit un-scaled model
glm.fits.unscaled <- glm(IUCNSkewPop ~ FedFundingTotal + YearsFunded,
                         family = binomial,
                         data = marine_log2)
summary(glm.fits.unscaled)

# Create a new data frame for predictions with a range of YearsFunded values
new_data <- data.frame(
  YearsFunded = seq(min(marine_log2$YearsFunded), 
                    max(marine_log2$YearsFunded), 
                    length.out = 100),
  FedFundingTotal = median(marine_log2$FedFundingTotal)
)

# Get predicted probabilities from the un-scaled model
new_data$pred_prob <- predict(glm.fits.unscaled, 
                              newdata = new_data, 
                              type = "response")

# Create the actual YearsFunded data points showing the observed outcomes
actual_points <- marine_log2 %>%
  select(YearsFunded, IUCNSkewPop) %>%
  mutate(
    # Add a small amount of jitter for better visualization
    prob_jitter = ifelse(IUCNSkewPop == 1, 0.98, 0.02)
  )

# Create the plot: YearsFunded
ggplot() +
  geom_line(data = new_data, 
            aes(x = YearsFunded, 
                y = pred_prob), 
            color = wes_palette("FantasticFox1")[3], 
            size = 1.2) +
  geom_point(data = actual_points, aes(x = YearsFunded, 
                                       y = prob_jitter), 
             shape = 16, 
             color = wes_palette("FantasticFox1")[4],
             size = 2) +  
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "black") +
  geom_hline(yintercept = 1, 
             linetype = "dashed", 
             color = "black") +
  labs(x = "Years Funded", 
       y = "Probability of Increasing Population Trend", 
       title = "Effect of Years Funded on Population Trend Probability") +
  theme_light() +
  scale_y_continuous(limits = c(0, 1))

# Create a new data frame for predictions with a range of FedFundingTotal values
new_data2 <- data.frame(
  FedFundingTotal = seq(min(marine_log2$FedFundingTotal), 
                        max(marine_log2$FedFundingTotal), 
                        length.out = 100),
  YearsFunded = median(marine_log2$YearsFunded)
)

# Get predicted probabilities from the unscaled model
new_data2$pred_prob <- predict(glm.fits.unscaled, 
                               newdata = new_data2, 
                               type = "response")

# Create the actual FedFundingTotal data points showing the observed outcomes
actual_points2 <- marine_log2 %>%
  select(FedFundingTotal, IUCNSkewPop) %>%
  mutate(prob_jitter = ifelse(IUCNSkewPop == 1, 0.98, 0.02)
  )

# Create the plot: FedFundingTotal
ggplot() +
  geom_line(data = new_data2, 
            aes(x = FedFundingTotal, 
                y = pred_prob), 
            color = wes_palette("FantasticFox1")[3], 
            size = 1.2) +
  geom_point(data = actual_points2, aes(x = FedFundingTotal, 
                                       y = prob_jitter), 
             shape = 16, 
             color = wes_palette("FantasticFox1")[4],
             size = 2) +  
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "black") +
  geom_hline(yintercept = 1, 
             linetype = "dashed", 
             color = "black") +
  labs(x = "Total Federal Funding", y = "Probability of Increasing Population Trend",
       title = "Effect of Federal Funding on Population Trend Probability") +
  theme_light() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(labels = label_dollar(scale = 1e-6, suffix = "M"))

############### Generate regression output table 

webshot::install_phantomjs(force = TRUE)

# Extract null and residual deviance for glm models
# Best fit glm
fits_nulldev <- glm.fits$null.deviance
fits_resdev <- glm.fits$deviance
# Base glm 
base_nulldev <- glm.base$null.deviance
base_resdev <- glm.base$deviance

# Create the table
stargazer(glm.base,
          glm.fits,
          lm.base.log,
          lm.fits,
          type = "html",
          title = "Table 1: Logistic & Linear Regression Models",
          out = "regression_table.html",
          omit.stat = c("aic"), # Omit auto-generated AIC line
          add.lines = list(
            c("AIC", round(AIC(glm.base), 2), # Add AIC lines
                     round(AIC(glm.fits), 2),
                     round(AIC(lm.base.log), 2),
                     round(AIC(lm.fits), 2)),
            c("Null Deviance", round(base_nulldev, 2), # Add null deviance lines
                               round(fits_nulldev, 2)),
            c("Residual Deviance", round(base_resdev, 2), # Add residual deviance lines 
                                   round(fits_resdev, 2))
          ))
webshot("regression_table.html", "regression_table.png")