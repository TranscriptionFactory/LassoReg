

required_pkgs = c("tidyverse",
                  "ggrepel", "gridExtra", "devtools",
                  "ggpubr", "stats", "ggpmisc",
                  "glmnet", "e1071", "caret", "randomForest",
                  "tree", "gbm", "matrixStats", "cvAUC",
                  "pROC", "pls",
                  "ROCR", "BAS", "gmp", "matrixcalc",
                   # "LassoReg", 
                  "EssReg")

for (p in required_pkgs) {
  installed = require(p, character.only = TRUE)
  if (p == "Essreg") {
    devtools::install_github("Hanxi-002/EssReg", force = FALSE,
                             dependencies = TRUE)
  }
  # else if (p == "LassoReg") {
  #   # devtools::install_github("TranscriptionFactory/LassoReg", force = FALSE,
  #   #                          dependencies = TRUE)
  # }
  else if (installed == FALSE) {
      install.packages(p, dependencies = TRUE)
  }
  
  library(p, character.only = TRUE)
  
}
