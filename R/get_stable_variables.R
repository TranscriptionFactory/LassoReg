
#' @export
get_stable_variables = function(lasso_results, freq_threshold = 0.50, max_vars = 100) {

  plot_list = list()
  # repeat for each lambda
  for ( var_list_num in 1:length(lasso_results$vars) ) {
    var_table = lasso_results$vars[[var_list_num]]$chosen_vars_freq

    # convert frequency to true frequency
    var_table$Freq = var_table$Freq / 100.0

    selected_vars = var_table[var_table$Freq > freq_threshold, 1]

    num_cols = round(max_vars / 20)


    max_vars = ifelse(max_vars > length(selected_vars), length(selected_vars), max_vars)


    # make it divisible by 5
    max_vars = max_vars - (max_vars %% 5)

    selected_vars = selected_vars[1:max_vars]

    # get average length of variables for plotting
    average_var_length = mean(stringr::str_length(selected_vars))


    selected_vars = as.data.frame(selected_vars)

    # # append column with positions for plotting variables
    selected_vars$pos = rep(1:num_cols, length.out = max_vars)

    # order so that we have columns grouped together
    selected_vars = selected_vars[order(selected_vars$pos), ]


    num_vars_each_col = max_vars / num_cols
    selected_vars$plot_y = rep(num_vars_each_col:1, length.out = max_vars)

    plt = selected_vars %>% ggplot2::ggplot(., aes(x = factor(pos),
                                                   y = plot_y,
                                                   label = selected_vars)) +
      ggplot2::geom_text() + ggplot2::theme_void() +
      ggplot2::theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
      ggplot2::ggtitle(paste0("Variables Selected by Lasso in ", freq_threshold, "% of folds in repeated cross-fold validation \n",
                              "lambda = ", lasso_results$lambdas[var_list_num])) +
      ggplot2::ylim(0, 20)

    # append to list
    plot_list[[var_list_num]] = plt
  }
  return(plot_list)
}
