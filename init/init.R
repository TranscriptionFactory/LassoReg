#!/usr/bin/env Rscript
args = R.utils::commandArgs(asValues = TRUE, excludeReserved = TRUE, args = TRUE)

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


# usage: init.R ~/path/to/data where input data is only file

############################
# init.R
# initialize working directory to data location and source
# run scripts
############################
# set directories
dirs = list(main = getwd(),
            scripts = paste0(getwd(), '/src/'),
                             data = paste0(getwd(), "/", dpath),
            cluster_filepath = cpath)

############################
# Code for running arguments from command line
############################
if (!file.exists(dirs$data)) {
  stop("Directory does not exist. Default is 'data' in current directory", call.=FALSE)
}

############################
# Source run scripts
############################
# first source imports
source(paste0(dirs$scripts, 'imports.R'))

# source imports and functions
for (file in list.files(dirs$scripts, full.names = T)) {
  source(file)
}

############################
# Process data and run lasso
############################

#default is to have data saved as csv where column 1 is Y values
# e.g data is [Y, X]
df = read_csv(dpath)

# test with 2-way
results = list()
file_id = ""
if (cpath != "") {
  modules = getModules(df, dirs$cluster_filepath)

  results = LASSO_Grid(modules$modulePA, alphaValues = c(0.75, 1.0, 1.25))

  # save these results
  # use regex to get the cluster file name wihtout the extension
  file_id = str_split_1(dirs$cluster_filepath,
                        pattern = regex("\\.[a-zA-Z0-9]+$"))[1]
  
  # get last element (don't want full file paht)
  file_id = tail(str_split_1(file_id, "/"), n = 1)

} else {
  # run normal lasso
  # run lasso/rf/svm

  alphaValues = c(0.75, 1, 1.25)
  
  results = LASSO_Grid(df, alphaValues)

  file_id = str_flatten(alphaValues, collapse = "_")

  }

saveRDS(results,
        file = paste0("results/", file_id, ".RDS"))

# plot results
plotResults(results, df)


# done
