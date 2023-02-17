#' @export
digestResults = function(results) {


  auc_matrix = data.frame("alpha",
                          "svm",
                          "svm_permute",
                          "svm_cfm",
                          "svm_cfm_permute",
                          "rf",
                          "rf_permute",
                          "rf_cfm",
                          "rf_cfm_permute",
                          "numVars", "entry", "sublist")

  for (entry in 1:length(results)) {
    for (sublist in 1:length(results[[entry]])) {
      currentList = results[[entry]][[sublist]]
      #
      auc_matrix[nrow(auc_matrix) + 1, ] = as.double(c(currentList$alpha,
                                                       format(currentList$svm$auc_model, digits = 3),
                                                       format(currentList$svm$auc_model_permute, digits = 3),
                                                       currentList$svm$cfm$overall[['Accuracy']],
                                                       currentList$svm$cfm_permute$overall[['Accuracy']],
                                                       format(currentList$rf$auc_model, digits = 3),
                                                       format(currentList$rf$auc_model_permute, digits = 3),
                                                       currentList$rf$cfm$overall[['Accuracy']],
                                                       currentList$rf$cfm_permute$overall[['Accuracy']],
                                                       length(currentList$chosenFeats), entry, sublist))
    }
  }

  names(auc_matrix) = auc_matrix[1,]

  auc_matrix = as.data.frame(auc_matrix[2:nrow(auc_matrix),])

  return(as.data.frame(lapply(auc_matrix, as.numeric)))
}

#' @export
extractVars = function(results) {

  allvars = list()

  for (index in results$lambdas) {
    vars_across_folds = c()

    for (entry in 1:length(results$gridResults)) {
      vars_across_folds = c(vars_across_folds,
                            results$gridResults[[entry]][[index]]$chosenFeats)
    }

    allvars[[index]] = vars_across_folds
  }
  return(allvars)
}


#' @export
getChullPolygon = function(data) {
  results = data.frame()
  for (group in levels(data$True)) {
    # filter each group
    df = data %>% dplyr::filter(True == group)

    outliers1 = graphics::boxplot(df$Comp1, plot = F)$out
    outliers2 = graphics::boxplot(df$Comp2, plot = F)$out

    chdf = df %>% dplyr::filter(!Comp1 %in% c(outliers1) & !Comp2 %in% c(outliers2))

    boundary = grDevices::chull(chdf)
    BumpX = chdf$Comp1[boundary] #+ 0.1*(df$Comp1[boundary] - mx)
    BumpY = chdf$Comp2[boundary] #+ 0.1*(df$Comp2[boundary] - my)

    # results[[paste0("l",group)]] = list(Comp1 = BumpX, Comp2 = BumpY, True = chdf$True[boundary])
    results = rbind(results, list(Comp1 = BumpX, Comp2 = BumpY, True = chdf$True[boundary]))
  }
  return(results)
}



plotResults = function(resultsdf, orig_df, outpath) {
  results = resultsdf$gridResults
  lambdas = resultsdf$lambdas
  df = orig_df

  auc_matrix = digestResults(results)

  ###### for LASSO ######

  vars = LassoReg::extractVars(results)[[1]] %>% table() %>%
    as.data.frame() %>% dplyr::arrange(desc(Freq)) #%>% filter(Freq > 5)

  # use these to subset the original data

  downselected = df[, names(df) %in% vars$.]
  downselected = cbind.data.frame(df$Group, downselected)
  names(downselected)[1] = "Group"

  plsr_vals = pls::plsr(Group ~ ., data = downselected, scale = T, validation = "CV")

  plsr_df = data.frame("Comp1" = plsr_vals$scores[,1], "Comp2" = plsr_vals$scores[,2],
                       "True" = as.factor(downselected$Group))

  chull_df = LassoReg::getChullPolygon(plsr_df)

  plot_plsr = ggplot2::ggplot(plsr_df, aes(Comp1, Comp2, True)) + ggplot2::theme_bw() +
    ggplot2::geom_polygon(data = data.frame(chull_df),  aes(x = Comp1, y = Comp2, fill = True), alpha = 0.25, inherit.aes = F) +
    ggplot2::geom_point(aes(color = True, fill = True), alpha = 1) +
    ggplot2::labs(x = 'PLS-DA Comp1', y = 'PLS-DA Comp2', title = 'PLS-DA using only Lasso Features') +
    ggplot2::theme(axis.text = element_text(size = 14))


  ggplot2::ggsave(paste0(outpath, "/plots/plot_plsr_lasso.png"), plot_plsr, height = 7, width = 7)


  #########################################################################
  ##### showing true vs permuted auc by alpha for both svm and rf
  # rf
  plot_auc = ggpubr::ggboxplot(auc_matrix %>%
                              tidyr::pivot_longer(cols = c("svm", "rf",
                                                    "svm_permute", "rf_permute"),
                                           names_to = "auc_method", values_to = "auc"),
                            x = "auc_method", y = "auc", fill = "auc_method", palette = "npg",
                            order = c("svm", "svm_permute", "rf", "rf_permute"),
                            facet.by = "alpha", repel = T,
                            add = "point", add.params = list(size = 2), font.label = list(size = 20),
                            title = "SVM or Random Forest Classification (AUC)") +
    ggpubr::stat_compare_means(method = "t.test", comparisons = list(c("svm", "svm_permute"),
                                                             c("rf", "rf_permute")))
  ggplot2::ggsave(paste0(outpath, "/plots/plot_auc_lasso.png"), plot = plot_auc, width = 14, height = 6)


  plot_cfm = ggpubr::ggboxplot(auc_matrix %>%
                             tidyr::pivot_longer(cols = c("svm_cfm", "rf_cfm",
                                                    "svm_cfm_permute", "rf_cfm_permute"),
                                           names_to = "auc_method", values_to = "auc"),
                            x = "auc_method", y = "auc", fill = "auc_method", palette = "npg",
                            order = c("svm_cfm", "svm_cfm_permute", "rf_cfm", "rf_cfm_permute"),
                            facet.by = "alpha", repel = T, xlab = "Method",
                            ylab = "Accuracy (Confusion matrix)",
                            add = "point", add.params = list(size = 2),
                            title = "SVM or Random Forest Classification (Confusion Matrix)") +
    ggpubr::stat_compare_means(method = "t.test", comparisons = list(c("svm_cfm", "svm_cfm_permute"),
                                                             c("rf_cfm", "rf_cfm_permute")))
  ggplot2::ggsave(paste0(outpath, "/plots/plot_cfm_lasso.png"), plot = plot_cfm, width = 16, height = 6)

}
