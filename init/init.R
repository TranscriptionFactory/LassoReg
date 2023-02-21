#!/usr/bin/env Rscript
args = R.utils::commandArgs(asValues = TRUE, excludeReserved = TRUE, args = TRUE)

library(usethis)
library(devtools)
devtools::install_github("TranscriptionFactory/LassoReg", force = T,
                         dependencies = F, quiet = F)

library(LassoReg, attach.required = F)
############################
# Code for running arguments from command line
############################

# for -d and --datapath and -c and --clusterpath
dpath_args = c("d", "datapath")
cpath_args = c("c", "clusterpath")

# leave cpath as "" if want to run normal Lasso (if not running from command line)
cpath = ""
dpath = "data/"

if (length(intersect(cpath_args, names(args))) > 0) {
  cpath = unlist(args[intersect(cpath_args, names(args))[1]])
}

if (length(intersect(dpath_args, names(args))) > 0) {
  dpath = unlist(args[intersect(dpath_args, names(args))[1]])
}


if (dpath != "data/exampledata" & !file.exists(dpath)) {
  stop("Directory does not exist. Default is 'data' in current directory", call.=FALSE)
}

dir.create(paste0(getwd(), "/results/"), showWarnings = F, recursive = T)
dir.create(paste0(getwd(), "/plots/"), showWarnings = F, recursive = T)

############################
# Process data and run lasso
############################

#default is to have data saved as csv where column 1 is Y values
# e.g data is [Y, X]

################################################################################
# Set data to use
################################################################################

if (dpath == "data/exampledata") {
  # no data path set, use example data
  df = LassoReg::exampledata
} else {
  df = readr::read_csv(dpath)

}
################################################################################

alphaValues = c(0.5)

# test with 2-way
results = list()
file_id = ""
if (cpath != "") {
  modules = LassoReg::getModules(df, cpath)

  results = LassoReg::LASSO_Grid(modules$modulePA, alphaValues)

  # save these results
  # use regex to get the cluster file name wihtout the extension
  file_id = stringr::str_split_1(cpath,
                                 pattern = regex("\\.[a-zA-Z0-9]+$"))[1]

  # get last element (don't want full file paht)
  file_id = tail(stringr::str_split_1(file_id, "/"), n = 1)

  saveRDS(modules, file = paste0(getwd(), "/results/modules.RDS"))

} else {
  # run normal lasso
  # run lasso/rf/svm
  results = LassoReg::LASSO_Grid(df, alphaValues)

  file_id = stringr::str_flatten(alphaValues, collapse = "_")

}

results$lambdas = alphaValues
chosen_vars_lambda = LassoReg::extractVars(results)

saveRDS(results,
        file = paste0(getwd(), "/results/", file_id, ".RDS"))

# saveRDS(chosen_vars_lambda, file = paste0(getwd(), "/results/chosen_features.RDS"))

# plot results
plotResults(results, df, paste0(getwd(), "/plots"))


# done
