
processFile = function(cluster_input, df, fullModules = F, filepath = NULL) {
  results = list()

  allModules = list()
  gcount = 0
  # read from filepath
  if (!is.null(filepath)) {
    con = base::file(filepath, "r")
    while ( TRUE ) {
      line = base::readLines(con, n = 1)
      if ( length(line) == 0 ) {
        break
      }
      genes = stringr::str_split(line, "\t|\n")[[1]]

      allModules[length(allModules) + 1] = list(l = genes)
      gcount = gcount + 1
    }
    close(con)

  } else if (!is.null(cluster_input)) {
      # read through clusters dataframe
      for (i in 1:length(cluster_input)) {
        allModules[length(allModules) + 1] = list( stringr::str_c(cluster_input[i, ]))
      }
  } else {
      # somethings wrong
      return()
  }

  # iterate through all modules and create the modules based on our tolerance

  for (i in 1:length(allModules)) {
    genes = allModules[[i]]
    genelist = c()
    for (g in genes) {
      if (g %in% names(df[, -1])) {
        genelist = c(genelist, g)
      }
    }
    if ( length(genelist) >= 3) {
      if (fullModules) {
        # get all genes in the modules
        results[length(results) + 1] = list(l = genes)
      } else {
        # get only genes in data from module
        results[length(results) + 1] = list(l = genelist)
      }
    }
  }

  final_results = list(results = results, allModules = allModules)
  return(final_results)
}

#' @export
getModules = function(dfs, cluster_filepath) {

  modules = processFile(cluster_filepath, dfs)

  modulesInData = modules$results
  tempModules <- modulesInData
  data = dfs

  M <- matrix(0, ncol=nrow(data), nrow = length(modulesInData))
  i <- 1
  Y = data$Group
  moduleNames = list()
  for (i in 1:length(tempModules)) {

    module <- modulesInData[[i]]
    # if(unique(module[module != ""]))
    module <- unique(module[module != ""])

    subDF <- data %>% dplyr::select(all_of(module))
    moduleNames[[length(moduleNames) + 1]] = module

    subDF <- cbind.data.frame(Y, subDF)


    grp_means = matrix(0, nrow = length(unique(Y)), ncol = ncol(subDF))
    grp_sds = matrix(0, nrow = length(unique(Y)), ncol = ncol(subDF))
    for (val in 1:length(unique(Y))) {
      grp_means[val,] = colSums(subset(subDF, Y == unique(Y)[val]))/(length(which(subDF$Y == unique(Y)[val])))

      grp_sds[val,] = apply(base::subset(subDF, Y == unique(Y)[val]), 2, sd)
    }

    grp_means = grp_means[, -1]

    grp_sds[, -1]
    tempSubDF <- subDF[,-1]
    m <- c()
    N <- c()
    SN <- c()
    tempSign <- 1
    for (k in 1:nrow(as.matrix(tempSubDF))) {
      for (j in 1:ncol(as.matrix(tempSubDF))) {
        # group means for col j
        j_means = sum(grp_means[, j]^2)

        # calculate sample vs group scores for each group
        # get value of Y for sample k
        y_k = Y[k]

#         j_score = ((rep(tempSubDF[k, j], times = length(grp_means[,j])) - grp_means[, j]))^2 / grp_means[, j]


        # get abs values

#         m[j] = sum(j_score)
          m[j] = log1p( ((tempSubDF[k, j] - grp_means[which(unique(Y) == y_k), j])^2) / (prod(grp_means[, j]) + 1) )

#         m[j] = log1p( (tempSubDF[k, j]^2) / prod(grp_means[, j]) )

        SN[j] <- sign(m[j])


      }

      m[!is.finite(m)] <- 0.0001

      posSN <- 0
      negSN <- 0
      posSN <- sum(SN > 0, na.rm = T)
      negSN <- sum(SN < 0, na.rm = T)

      m  <- abs(m)

      if (posSN > negSN){
        N[k] <- sum(m)*1
      } else{
        N[k] <- sum(m)*(-1)
      }
    }

    M[i, ] <- N
  }

  modulePA <- t(M)  #%>% data.frame(stringsAsFactors = F)

  # module names correspond to the genes in each module
  names(modulePA) = 1:length(moduleNames)
  modulePA = apply(modulePA, 2, function(x) scale(x, F, T))
  # modulePA = modulePA[, !duplicated(names(modulePA))]
  Group = data$Group
  modulePA <- cbind.data.frame(Group, modulePA)


  retenv = list()
  retenv$modulePA = modulePA
  retenv$modules = modules
  return(retenv)

}
