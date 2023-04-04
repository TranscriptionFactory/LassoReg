#' @export
calculateAUC = function(trueY, permutedY,
                        append_model, append_model_permute,
                        unique_levels = 2) {

    if (length(unique(trueY)) == 2) {

    auc_model = glmnet:::auc(trueY, append_model)

    # auc_model_permute = auc(permutedY, append_model_permute)
    auc_model_permute = glmnet:::auc(permutedY, append_model_permute)

  }
  else {
    # calculate multiclass auc
    auc_model = pROC::auc(pROC::multiclass.roc(trueY, append_model))

    # auc_model_permute = auc(multiclass.roc(permutedY, append_model_permute))
    auc_model_permute = pROC::auc(pROC::multiclass.roc(permutedY, append_model_permute))

  }

  trueY = factor(trueY)
  permutedY = factor(permutedY)
#   cfm = caret::confusionMatrix(table(factor(round(append_model)),
#                                          trueY))

#   cfm_permute = caret::confusionMatrix(table(factor(round(append_model_permute)),
#                                 permutedY))
  # these might break if the classes need to be rounded or round incorrectly
  cfm = caret::confusionMatrix(table(factor(round(append_model), levels = levels(trueY)),
                                         trueY))

  cfm_permute = caret::confusionMatrix(table(factor(round(append_model_permute), levels = levels(permutedY)),
                                permutedY))

  return(list(auc_model = auc_model, auc_model_permute = auc_model_permute,
              cfm = cfm, cfm_permute = cfm_permute))
}


#' @export
LASSO_Grid = function(fulldata, lambdaValues = c(1.0)) {
  retenv = list()
  #############################################################

  # set first column to be Y (called Group here)
  names(fulldata)[1] = "Group"

  lambdaValues = as.numeric(lambdaValues)

  # number of folds to create
  k = 10

  # number of repeats
  numRuns = 10

  # initialize results list
  gridResults = vector("list", numRuns)

  for(runCount in 1:numRuns) {

    gridValues = list()
    for (val in lambdaValues) {
      gridValues[[length(gridValues) + 1]] = list(alpha = val)
    }

    #create kfold associate to a data frame.
    kfolds <- caret::createFolds(y = fulldata$Group, k=k, list = FALSE, returnTrain = FALSE)
    myData <- cbind(fulldata, kfolds)

    # scramble response vector
    permutedY = sample(fulldata$Group)

    unique_levels = length(unique(fulldata$Group))

    ##############################################

    for (fold in 1:k) {
      #Create different subdata by different fold.
      currentFold <- which(kfolds == fold)

      # copy original data
      newData = myData
      variables = list()

      for (entry in 1:length(gridValues)) {

        for (n_run in 1:2) {
          #take out the random row in newData and put into train,
          #and the test data only contain the one which the train data does not have
          train <- newData[-currentFold, ]
          test <- newData[currentFold, ]

          #remove the kfold column which is the last column out of the data to calculate
          train <- train[, -ncol(train)]
          test <- test[, -ncol(test)]

          #Since X data will be constant
          X_train <- as.matrix(train[, -1])

          # not used until later
          X_test <- as.matrix(test[, -1])

          y = train$Group

          lasso_run <- glmnet::cv.glmnet(x=X_train, y=y, alpha=1, nfolds = k,
                                 standardize = T)


          #run lasso on training data with best lambda
          lasso_run_lmin <- glmnet::glmnet(x=X_train, y=y, alpha= 1.0,
                                   lambda = lasso_run$lambda.min * gridValues[[entry]]$alpha)
          c = stats::coef(lasso_run_lmin)
          inds = which(c!=0)


          if (length(inds) == 0) {
            print("No variables selected\n")
            next
          }

            #remove intercept from the list
            variables[[entry]] <- tail(row.names(c)[inds], -1)

            selTrain = train[, (names(train) %in% variables[[entry]])]
            newTrain = cbind.data.frame(train$Group, selTrain)

            selTest = test[, (names(test) %in% variables[[entry]])]
            newTest = cbind.data.frame(test$Group, selTest)

            colnames(newTrain)[1] <- "Y"
            colnames(newTest)[1] <- "Y"

            # for true Y and permuteY, put in separate lists but both named Y
            if (n_run == 1) {
              gridValues[[entry]]$trueY = append(gridValues[[entry]]$trueY, newTest$Y)

            } else {
              gridValues[[entry]]$permuteY = append(gridValues[[entry]]$permuteY, newTest$Y)

            }

            newTrain = droplevels(newTrain)
            newTest = droplevels(newTest)

            svmfit = e1071::svm(y = as.factor(newTrain$Y), x=as.matrix(newTrain[,-1]),
                         kernel="linear", cost=10, scale=T, na.action = na.omit)
            yhat.SVM = predict(svmfit, newdata = as.matrix(newTest[, -1]), type = "response")

              # random forest
            RFfit <- randomForest::randomForest(y = as.factor(newTrain$Y), x=as.matrix(newTrain[,-1]),
                                  importance=TRUE, ntree = 100)

            yhat.RF = predict(RFfit, newdata = as.matrix(newTest[, -1]), type = "response")

            # store yhat in appropriate location depending on mode
            if (n_run == 1) {
              gridValues[[entry]]$chosenFeats = variables[[entry]]

              #svm and rf
              gridValues[[entry]]$append_svm = c(gridValues[[entry]]$append_svm, yhat.SVM)
              gridValues[[entry]]$svm_models = list(fit = svmfit, variables = variables[[entry]])

              gridValues[[entry]]$append_RF = c(gridValues[[entry]]$append_RF, yhat.RF)
              gridValues[[entry]]$rf_models = list(fit = RFfit, variables = variables[[entry]])

            } else {
              gridValues[[entry]]$chosenFeat_perm = variables[[entry]]

              gridValues[[entry]]$append_svm_permute = c(gridValues[[entry]]$append_svm_permute, yhat.SVM)

              gridValues[[entry]]$append_RF_permute = c(gridValues[[entry]]$append_RF_permute, yhat.RF)
          }
        }
      }
    }

    # calculate AUC and confusion matrices for svm and RF models
    for (entry in 1:length(gridValues)) {

      gridValues[[entry]]$svm = calculateAUC(gridValues[[entry]]$trueY,
                                             permutedY,
                                             gridValues[[entry]]$append_svm,
                                             gridValues[[entry]]$append_svm_permute,
                                             unique_levels)

      gridValues[[entry]]$rf = calculateAUC(gridValues[[entry]]$trueY,
                                             permutedY,
                                             gridValues[[entry]]$append_RF,
                                             gridValues[[entry]]$append_RF_permute,
                                            unique_levels)

    }

    gridResults[[runCount]] = gridValues
  }

  retenv$gridResults = gridResults
  retenv$lambdas = lambdaValues
  retenv$lasso_input = fulldata

  # extractVars here
  retenv$vars = LassoReg::extractVars(retenv)
  return(retenv)
}
