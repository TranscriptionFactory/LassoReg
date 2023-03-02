# load results file output from LASSO_GRID


# extract classification data from each repeat
library(caret)
library(glmnet)
library(tidyverse)

# results = LassoReg::LASSO_Grid(modulePA, alphaValues = alphaValues)


true_confusionMatrices = list()
permuted_confusionMatrices = list()

plottingData = data.frame()

convert_Factor_To_Numeric_To_Run_auc = function(trueY, predY) {
  return(glmnet:::auc(as.numeric(as.character(trueY)), as.numeric(as.character(predY))))
}


for (res in 1:length(results$gridResults)) {
  gridResults = results$gridResults
  for (lambdaValue in 1:length(res)) {

    trueY = factor(gridResults[[res]][[lambdaValue]]$trueY)

    # scramble the permuted
    permutedY = sample(factor(gridResults[[res]][[lambdaValue]]$trueY))

    # svm and rf predictions to factors (need to have same levels as true Y values)
    svm_Y = factor(gridResults[[res]][[lambdaValue]]$append_svm, levels = levels(trueY))
    rf_Y = factor(gridResults[[res]][[lambdaValue]]$append_RF, levels = levels(trueY))

    permuted_svm_Y = factor(gridResults[[res]][[lambdaValue]]$append_svm_permute, levels = levels(permutedY))
    permuted_rf_Y = factor(gridResults[[res]][[lambdaValue]]$append_RF_permute, levels = levels(permutedY))

    # recalculate confusion matrices and AUC

    temp_data = list("run" = res,
      "lambda" = results$lambdas[lambdaValue],
      "svm_cfm" = caret::confusionMatrix(svm_Y, trueY)$overall[['Accuracy']],
      "rf_cfm" = caret::confusionMatrix(rf_Y, trueY)$overall[['Accuracy']],
      "svm_cfm_permute" = caret::confusionMatrix(permuted_svm_Y, permutedY)$overall[['Accuracy']],
      "rf_cfm_permute" = caret::confusionMatrix(permuted_rf_Y, permutedY)$overall[['Accuracy']],
      "svm_auc" = convert_Factor_To_Numeric_To_Run_auc(trueY, svm_Y),
      "rf_auc" = convert_Factor_To_Numeric_To_Run_auc(trueY, rf_Y),
      "svm_auc_permute" = convert_Factor_To_Numeric_To_Run_auc(permutedY, permuted_svm_Y),
      "rf_auc_permute" = convert_Factor_To_Numeric_To_Run_auc(permutedY, permuted_rf_Y))

    plottingData = rbind.data.frame(plottingData, temp_data)
  }
}

# boxplots

plot_auc = ggpubr::ggboxplot(plottingData %>%
                               tidyr::pivot_longer(cols = c("svm_auc", "rf_auc",
                                                            "svm_auc_permute", "rf_auc_permute"),
                                                   names_to = "auc_method", values_to = "auc"),
                             x = "auc_method", y = "auc", fill = "auc_method", palette = "npg",
                             order = c("svm_auc", "svm_auc_permute", "rf_auc", "rf_auc_permute"),
                             facet.by = "lambda", repel = T,
                             add = "point", add.params = list(size = 2), font.label = list(size = 20),
                             title = "SVM or Random Forest Classification (AUC)") +
  ggpubr::stat_compare_means(method = "t.test", comparisons = list(c("svm_auc", "svm_auc_permute"),
                                                                   c("rf_auc", "rf_auc_permute")))

plot_auc

plot_cfm = ggpubr::ggboxplot(plottingData %>%
                               tidyr::pivot_longer(cols = c("svm_cfm", "rf_cfm",
                                                            "svm_cfm_permute", "rf_cfm_permute"),
                                                   names_to = "cfm_accuracy", values_to = "accuracy"),
                             x = "cfm_accuracy", y = "accuracy", fill = "cfm_accuracy", palette = "npg",
                             order = c("svm_cfm", "svm_cfm_permute", "rf_cfm", "rf_cfm_permute"),
                             facet.by = "lambda", repel = T, xlab = "Method",
                             ylab = "Accuracy (Confusion matrix)",
                             add = "point", add.params = list(size = 2),
                             title = "SVM or Random Forest Classification (Confusion Matrix)") +
  ggpubr::stat_compare_means(method = "t.test", comparisons = list(c("svm_cfm", "svm_cfm_permute"),
                                                                   c("rf_cfm", "rf_cfm_permute")))

plot_cfm



