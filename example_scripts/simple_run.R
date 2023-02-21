library(tidyverse)
library(usethis)
library(devtools)
devtools::install_github("TranscriptionFactory/LassoReg", force = T,
                         dependencies = T, quiet = F)

library(LassoReg, attach.required = T)
library(Matrix)

########################################################
# 1. Load data as dataframe with Y values as first column
########################################################

df = as.data.frame(LassoReg::exampledata)

alphaValues = c(0.75, 1.0, 1.25)

########################################################
# 2. Run Lasso with multipliers on the lambda chosen from CV
########################################################

results = LassoReg::LASSO_Grid(df, alphaValues = alphaValues)


########################################################
# 3. Extract chosen features and append to results
########################################################

chosen_vars = LassoReg::extractVars(results)

results$chosen_vars = chosen_vars

saveRDS(results, 'results/results.RDS')

########################################################
# 4. Plot results
########################################################

LassoReg::plotResults(results, df, 'plots/')
