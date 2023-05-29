#' @export
calculateAUC = function(trueY, permutedY,
                        append_model, append_model_permute,
                        unique_levels = 2) {

  trueY = as.factor(trueY)
  permutedY = as.factor(permutedY)
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

  cfm = caret::confusionMatrix(data = factor(round(append_model)),
                                                   reference = trueY)

  cfm_permute = caret::confusionMatrix(data = factor(round(append_model_permute)),
                                       reference = permutedY)

  return(list(auc_model = auc_model, auc_model_permute = auc_model_permute,
              cfm = cfm, cfm_permute = cfm_permute))
}

#' @export
LASSO_Grid = function(fulldata, lambdaValues = c(1.0), numFolds = 10, numRuns = 10) {
  # all results to be returned
  retenv = list()
  #############################################################

  # need to remove columns with all zeros and zero standard deviation
  temp = fulldata[, -1]
  temp = temp[, which(apply(temp, 2, sd) != 0)]

  fulldata = cbind.data.frame(fulldata[, 1], temp)

  # set first column to be Y (called Group here)
  names(fulldata)[1] = "Group"

  lambdaValues = as.numeric(lambdaValues)

  # number of repeats
  numRuns = ifelse(numRuns > 0, numRuns, 10)

  # initialize results list
  gridResults = vector("list", numRuns)

  # start replicates of cv
  for(runCount in 1:numRuns) {

    gridValues = list()
    for (val in lambdaValues) {
      gridValues[[length(gridValues) + 1]] = list(lambda = val)
    }

    #create folds for dataframe and attach to data
    kfolds <- caret::createFolds(y = fulldata$Group, k=numFolds, list = FALSE, returnTrain = FALSE)
    myData <- cbind(fulldata, kfolds)


    # true response
    trueY = fulldata$Group
    # scramble response vector
    permutedY = sample(fulldata$Group)

    # copy original data
    newData_true_y = myData

    # create dataframe with permuted Y as group
    newData_permuted_y = myData
    newData_permuted_y$Group = permutedY

    unique_levels = length(unique(fulldata$Group))

    ##############################################
    # begin kfold validation
    for (fold in 1:numFolds) {
      #Create different subdata by different fold.
      currentFold <- which(kfolds == fold)


      # repeat the lasso selection/svm/rf for each lamda value
      for (entry in 1:length(gridValues)) {

        # this function runs lasso, selects variables and fits the svm/rf
        select_lasso_variables_fit_model = function(newData, currentFold, lambda_val, numFolds) {

          lasso_run_results = list()

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

          lasso_run <- glmnet::cv.glmnet(x=X_train, y=y, alpha=1, nfolds = numFolds,
                                         standardize = T)

          #run lasso on training data with best lambda
          lasso_run_lmin <- glmnet::glmnet(x=X_train, y=y, alpha= 1.0,
                                           lambda = lasso_run$lambda.min * lambda_val)

          # take the indices of the coefficients that are nonzero
          # (zero is the intercept value for the regression fit)
          c = stats::coef(lasso_run_lmin)
          inds = c@i[-1]

          if (length(inds) == 0) {
            print("No variables selected\n")
            return(NULL)
          }

          #remove intercept from the list
          lasso_run_results$variables <- colnames(X_train)[inds]

          selTrain = train[, (names(train) %in% lasso_run_results$variables)]
          newTrain = cbind.data.frame(train$Group, selTrain)

          selTest = test[, (names(test) %in% lasso_run_results$variables)]
          newTest = cbind.data.frame(test$Group, selTest)

          colnames(newTrain)[1] <- "Y"
          colnames(newTest)[1] <- "Y"

          lasso_run_results$Y = newTest$Y

          newTrain = droplevels(newTrain)
          newTest = droplevels(newTest)

          svmfit = e1071::svm(y = as.factor(newTrain$Y), x=as.matrix(newTrain[,-1]),
                              kernel="linear", cost=10, scale=T, na.action = na.omit)
          yhat.SVM = predict(svmfit, newdata = as.matrix(newTest[, -1]), type = "response")

          # random forest
          RFfit <- randomForest::randomForest(y = as.factor(newTrain$Y), x=as.matrix(newTrain[,-1]),
                                              importance=TRUE, ntree = 100)

          yhat.RF = predict(RFfit, newdata = as.matrix(newTest[, -1]), type = "response")


          # save model predictions
          lasso_run_results$yhat.SVM = yhat.SVM
          lasso_run_results$yhat.RF = yhat.RF

          # save models
          lasso_run_results$svmfit = svmfit
          lasso_run_results$RFfit = RFfit

          return(lasso_run_results)
        }

        current_lambda_val = gridValues[[entry]]$lambda

        # here is where we call the function above for data with the true Y and
        # data with the permuted Y
        ######################################################
        # DATA WITH TRUE Y
        ######################################################
        # first run lasso/svm/rf for data with the true Y values
        lasso_run_results_trueY = select_lasso_variables_fit_model(newData_true_y, currentFold, current_lambda_val, numFolds = numFolds)

        if (is.null(lasso_run_results_trueY)) {
          print("no variables selected using true Y; skipping fold")
          next
        }
        # now assign these results to our big list
        gridValues[[entry]]$chosenFeats = c(gridValues[[entry]]$chosenFeats, lasso_run_results_trueY$variables)

        #svm and rf
        gridValues[[entry]]$append_svm = c(gridValues[[entry]]$append_svm, lasso_run_results_trueY$yhat.SVM)
        gridValues[[entry]]$svm_models = list(fit = lasso_run_results_trueY$svmfit, variables = lasso_run_results_trueY$variables)

        gridValues[[entry]]$append_RF = c(gridValues[[entry]]$append_RF, lasso_run_results_trueY$yhat.RF)
        gridValues[[entry]]$rf_models = list(fit = lasso_run_results_trueY$RFfit, variables = lasso_run_results_trueY$variables)
        gridValues[[entry]]$trueY = c(gridValues[[entry]]$trueY, lasso_run_results_trueY$Y)

        ######################################################
        # DATA WITH PERMUTED Y
        ######################################################
        lasso_run_results_permutedY = select_lasso_variables_fit_model(newData_permuted_y, currentFold, current_lambda_val, numFolds = numFolds)

        if (is.null(lasso_run_results_permutedY)) {
          print("no variables selected using permuted Y; skipping fold")
          next
        }
        # assign these to the spots for our permuted results
        gridValues[[entry]]$chosenFeat_perm = c(gridValues[[entry]]$chosenFeat_perm, lasso_run_results_permutedY$variables)
        gridValues[[entry]]$append_svm_permute = c(gridValues[[entry]]$append_svm_permute, lasso_run_results_permutedY$yhat.SVM)
        gridValues[[entry]]$append_RF_permute = c(gridValues[[entry]]$append_RF_permute, lasso_run_results_permutedY$yhat.RF)
        gridValues[[entry]]$permutedY = c(gridValues[[entry]]$permutedY, lasso_run_results_permutedY$Y)
      }
    }

    # calculate AUC and confusion matrices for svm and RF models
    for (entry in 1:length(gridValues)) {

      # we evaluate both models (model trained on true Y and model trained on permuted Y)
      # against the true Y
      gridValues[[entry]]$svm = LassoReg::calculateAUC(gridValues[[entry]]$trueY,
                                             gridValues[[entry]]$trueY,
                                             gridValues[[entry]]$append_svm,
                                             gridValues[[entry]]$append_svm_permute,
                                             unique_levels)

      gridValues[[entry]]$rf = LassoReg::calculateAUC(gridValues[[entry]]$trueY,
                                            gridValues[[entry]]$trueY,
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
