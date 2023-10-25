library(R6)
library(xgboost)
library(Matrix)
library(caret)
library(e1071)
library(gtools)
library(tidyverse)
library(randomForest)
library(ROCR)
library(tensorflow)
library(reticulate)
library(keras)
#use_condaenv(condaenv = "tensorflow_cpu", conda = "/Users/vbucci/anaconda3/bin/conda")
#use_condaenv(condaenv = "tensorflow_cpu", conda = "/data/tools/miniconda3/bin/conda")
#use_condaenv(condaenv = "tensorflow_cpu", conda = "/home/vbucci/anaconda3/bin/conda")

library(vita)
library(plyr)
library(dplyr)
library(caret)
library(yardstick)
library(lime)
library(pROC)
library(ranger)
library(rsample)
library(Boruta)
library(MLmetrics)

#library(MEmlm) # this is a custom version of the MEml package that I mad. Will share if anyone wants it!

# get macro f1- score 
f1_score <- function(predicted, expected, positive.class="1") {
  predicted <- factor(as.character(predicted), levels=unique(as.character(expected)))
  expected  <- as.factor(expected)
  cm = as.matrix(table(expected, predicted))
  
  precision <- diag(cm) / colSums(cm)
  recall <- diag(cm) / rowSums(cm)
  f1 <-  ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  
  #Assuming that F1 is zero when it's not possible compute it
  f1[is.na(f1)] <- 0
  
  #Binary F1 or Multi-class macro-averaged F1
  ifelse(nlevels(expected) == 2, f1[positive.class], mean(f1))
  return(f1)
}

'%notin%' <- Negate('%in%')

PFI = function(X, Y, f){
  score = c()
  FI = c()
  
  #score_orig = f %>% evaluate(
  #  X, Y, verbose = 0)
  
  for (i in 1: ncol(X)){
    xi = sample(X[,i], replace = FALSE)
    Xi = cbind(xi, X[,-i])
    
    score[i] = f %>% evaluate(
      Xi, Y, verbose = 0)

    #FI [i] = score[i]/score_orig
    FI [i] = score[i]
  }
  PFI_final = data.frame("Features" = colnames(X), "FI" = FI)
  return(PFI_final[order(PFI_final$Features, decreasing = TRUE),])
}

ml_data <- R6Class("ml_data",
                  public = list(
                    all_frames = NULL,
                    fulldata = NULL,
                    xdata = NULL,
                    ydata = NULL,
                    ptype = NULL,
                    train_data = NULL,
                    test_data = NULL,
                    
                    # boruta
                    boruta_model = NULL,
                    boruta_importance = NULL,
                    
                    #rfr
                    rfr_model = NULL,
                    rfr_importance = NULL,
                    rfr_confusion = NULL, 
                    rfr_roc = NULL, 
                    rfr_rmse = NULL,
                    rfr_lime_explanation = NULL,
                    
                    # rfc
                    rfc_model = NULL,
                    rfc_importance = NULL,
                    rfc_confusion = NULL,
                    rfc_pred = NULL,
                    rfc_roc = NULL, 
                    rfc_f1 = NULL,
                    rfc_lime_explanation = NULL,
                    
                    #xgb
                    xgb_importance = NULL,
                    xgb_confusion = NULL,
                    xgb_roc = NULL,
                    xgb_lime_explanation = NULL, 
                    
                    # nn
                    keras_model = NULL,
                    keras_lime_explanation = NULL,
                    keras_confusion = NULL, 
                    keras_accuracy = NULL,
                    keras_roc = NULL,
                    keras_estimates_tbl = NULL,
                    keras_pfi_importance = NULL,
                    
                    # ordinal RF
                    of_model = NULL,
                    of_importance = NULL,
                    of_confusion = NULL,
                    of_roc = NULL,
                    of_f1= NULL,
                    of_lime_explanation = NULL,
                    
                    # MERF
                    merf_model = NULL,
                    merf_importance = NULL,
                    merf_confusion = NULL,
                    merf_pred = NULL,
                    merf_roc = NULL,
                    merf_f1 = NULL,
                    merf_out = NULL, 
                    merf_rhs = NULL,
                    
                    # merf_lime_explanation = NULL, # havent figured out how to do this yet
                    
                    # rf multiclass classification
                    rfmc_model = NULL,
                    rfmc_importance = NULL,
                    rfmc_confusion = NULL, 
                    rfmc_roc = NULL, 
                    rfmc_f1 = NULL,
                    rfmc_lime_explanation = NULL,
                    
                    initialize = function(all_frames = NA, 
                                          fulldata = NA, xdata = NA, ydata = NA, ptype = NA) {
                      self$all_frames <- all_frames
                      self$fulldata <- fulldata
                      self$xdata <- xdata
                      self$ydata <- ydata
                      self$ptype <- ptype
                    },
                    
                    load_fulldata_folder = function(indir,c_var) {
                      s_files <- mixedsort(list.files(path = paste0(indir,c_var),pattern = ".csv",full.names = T))
                      # Read all the csv files in a given path and create a column for the source information (Sample =  1,2,3...)
                      all_frames <- 1:length(s_files)  %>%
                        map_dfr(function(x) {
                          read.csv(s_files[x],row.names = 1)%>% mutate(Sample = x ) } )
                      print(all_frames)
                      self$all_frames <- all_frames
                    },
                    
                    load_fulldata = function(batch_id) {
                      dem_dt <-  self$all_frames %>% filter(Sample ==  batch_id)
                      print(dem_dt)
                      dem_dt$Sample <- NULL
                      dem_dt$ID <- NULL
                      self$fulldata <- dem_dt
                    },
                    
                    load_fulldata_list = function(df){
                      dem_dt <- df
                      dem_dt$Sample <- NULL
                      dem_dt$ID <- NULL
                      self$fulldata <- dem_dt
                    },
                    
                    factorize_predictors = function(to_be_factorized) {
                      self$fulldata[to_be_factorized] <- lapply(self$fulldata[to_be_factorized] , factor)
                    },
                    
                    divide_train_and_test = function() {
                      train_dt <- self$fulldata[self$fulldata$Set == "train",]
                      train_dt$Set <-  NULL
                      test_dt <-  self$fulldata[self$fulldata$Set == "test",]
                      test_dt$Set <- NULL
                      self$train_data <- train_dt
                      self$test_data <- test_dt
                    },
                
                    run_rf_regression = function(myseed, ntrees, ncores, opt_encoding, is_splitted_data, o_run_pimp) {
                      ## YET TO CONFIRM THIS WORKS
                      # one hot encoding for train set
                      if (opt_encoding == TRUE){
                        # dummify the data (One-hot-encoding)
                        dmy <- dummyVars(" ~ .", data = self$train_data)
                        train_dt_rf <- data.frame(predict(dmy, newdata = self$train_data))
                      }else{
                        train_dt_rf <- self$train_data
                      }
                      
                      if (is_splitted_data==T){
                        # one hot encoding for test set
                        if (opt_encoding == TRUE){
                          # dummify the data (One-hot-encoding)
                          dmy <- dummyVars(" ~ .", data = self$test_data)
                          test_dt_rf <- data.frame(predict(dmy, newdata = self$test_data))
                        }else{
                          test_dt_rf <- self$test_data
                        }
                      }
                      # Random Forest
                      # Run RF regression
                      set.seed(myseed)
                      #rf <- randomForest(class ~.,train_dt_rf, ntree = ntrees, importance = T)
                      rf <- ranger(class~.,data=train_dt_rf,importance = "permutation",
                                   probability = T,
                                   num.trees = 5000,
                                   replace=FALSE)
                      self$rfr_model <- rf
                      
                      print(rf)
                      print("RF regression training done....")
                      
                      if (o_run_pimp == TRUE) {
                        # Run permutated importance 
                        pimp <- PIMP(train_dt_rf[,-c("class")], train_dt_rf$class, rf, parallel=TRUE, ncores = ncores, seed = myseed)
                        p_test <- PimpTest(pimp)
                        pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
                        imp_dt <- pimp_all     
                        imp_dt$Predictors <- rownames(imp_dt)
                        rownames(imp_dt) <- NULL
                        self$rfr_importance<-imp_dt
                      }else {
                        self$rfr_importance <- NULL
                      }
                      
                      if (is_splitted_data == T){
                        p_rf_test <- predict(self$rfr_model, test_dt_rf)
                        print(p_rf_test$predictions)
                        rmse <- rmse(p_rf_test$predictions, test_dt_rf$class)
                        print(rmse)
                        self$rfr_rmse <- rmse
                        
                        # lime explanation
                        train_X <-  train_dt_rf[,!colnames(train_dt_rf) %in% "class"]
                        rf <- as_regressor(rf,labels = NULL)
                        expln <- lime(train_X, model = rf)
                        #print(expln)
                        test_X <-  test_dt_rf[,!colnames(test_dt_rf) %in% "class"]
                        lime_reasons <- lime::explain(x=test_X, explainer=expln,
                                                      n_labels = 1,
                                                      n_features = 4,
                                                      feature_select = "tree")
                        self$rfr_lime_explanation <- lime_reasons
                      }
                    },
                    
                    run_rf_classification = function(myseed, ntrees, ncores, opt_encoding, 
                                                     is_splitted_data, o_run_pimp, is_LOO) {
                      
                      # one hot encoding for train set
                      if (opt_encoding == TRUE){
                        # dummify the data (One-hot-encoding)
                        dmy <- dummyVars(" ~ .", data = self$train_data)
                        train_dt_rf <- data.frame(predict(dmy, newdata = self$train_data))
                        train_dt_rf$class <- factor(train_dt_rf$class, levels = c(0,1))
                      }else{
                        train_dt_rf <- self$train_data
                        train_dt_rf$class <- factor(train_dt_rf$class, levels = c(0,1))
                      }
                      
                      if (is_splitted_data==T){
                        # one hot encoding for test set
                        if (opt_encoding == TRUE){
                          # dummify the data (One-hot-encoding)
                          dmy <- dummyVars(" ~ .", data = self$test_data)
                          test_dt_rf <- data.frame(predict(dmy, newdata = self$test_data))
                          test_dt_rf$class <- factor(test_dt_rf$class,levels = c(0,1))
                        }else{
                          train_dt_rf <- self$train_data
                          train_dt_rf$class <- factor(train_dt_rf$class, levels = c(0,1))
                          test_dt_rf <- self$test_data
                          test_dt_rf$class <- factor(test_dt_rf$class,levels = c(0,1))
                        }
                      }
                      # Random Forest
                      # Run RF classification 
                      set.seed(myseed)
                      #rf <- randomForest(class ~.,train_dt_rf, ntree = ntrees, importance = T)
                      rf <- ranger(class~.,data=train_dt_rf,importance = "permutation",
                                         probability = T,
                                         num.trees = 5000,
                                         replace=FALSE)
                      self$rfc_model <- rf
                      
                      print(rf)
                      print("RF classification training done....")
                      
                      if (o_run_pimp == TRUE) {
                        # Run permutated importance 
                        X <- train_dt_rf %>% select(-class)
                        pimp <- PIMP(X, train_dt_rf$class, rf, parallel=TRUE, ncores = ncores, seed = myseed)
                        p_test <- PimpTest(pimp)
                        pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
                        imp_dt <- pimp_all     
                        imp_dt$Predictors <- rownames(imp_dt)
                        rownames(imp_dt) <- NULL
                        self$rfc_importance<-imp_dt
                      }else {
                        self$rfc_importance <- NULL
                      }
                      
                      if (is_splitted_data == T){
                        
                        if(is_LOO == T){
                          print("Is LOOCV, running prediction and saving results...")
                          p_rf_test <- predict(self$rfc_model, test_dt_rf)
                          self$rfc_pred <- p_rf_test
                          
                          print(p_rf_test$predictions)
                         
                        }else{
                          p_rf_test <- predict(self$rfc_model, test_dt_rf)
                          self$rfc_pred <- p_rf_test
                          
                          print(p_rf_test$predictions)
                          #pred <-  p_rf_test  %>%
                          #  data.frame() %>%
                          #  mutate(label = test_dt_rf$class,max_prob = max.col(.,"last") -1)
                          #res_pred = ifelse(pred[,"X1"] >= 0.5,1,0)
                          res_pred = ifelse(p_rf_test$predictions[,1] >= 0.5,1,0)
                          #print(res_pred)
                          f1_df<-data.frame(pred= res_pred, test= ifelse(test_dt_rf$class == 1,1,0))
                          print(f1_df)
                          self$rfc_f1 <- f1_df
                      
                          ## cant run lime with only one variable so have to add an if statement
                          
                          if(length(train_dt_rf) >= 3){
                          # lime explanation
                            train_X <-  train_dt_rf %>% select(-class)
                            rf <- as_classifier(rf,labels = NULL)
                            expln <- lime(train_X, model = rf)
                            #print(expln)
                            test_X <-  test_dt_rf %>% select(-class)
                            lime_reasons <- lime::explain(x=test_X, explainer=expln,
                                                          n_labels = 1,
                                                          n_features = 4,
                                                          feature_select = "tree")
                            self$rfc_lime_explanation <- lime_reasons
                          }else{
                            self$rfc_lime_explanation <- "Not enough variables maintained to run lime"
                          }
                        }
                      }
                    },
                    
                    do_rsample_split_fulldata = function(v_in, repeats_in, fold_id){
                      folds<- vfold_cv(self$fulldata, v= v_in, repeats= repeats_in) 
                      split<-folds$splits[[fold_id]]
                      train_dt <- analysis(split)
                      test_dt <- assessment(split)
                      self$train_data <- train_dt
                      self$test_data <- test_dt
                    },
                    
                    run_boruta = function(myseed, doTrace, maxRuns, opt_encoding, is_splitted_data){
                      
                      if (opt_encoding == TRUE){
                        # dummify the data (One-hot-encoding)
                        dmy <- dummyVars(" ~ .", data = self$train_data)
                        train_dt_rf <- data.frame(predict(dmy, newdata = self$train_data))
                        train_dt_rf$class <- factor(train_dt_rf$class, levels = c("0","1"))
                      }else{
                        train_dt_rf <- self$train_data
                        train_dt_rf$class <- as.factor(train_dt_rf$class)
                      }
                      
                      if (is_splitted_data == T){
                        # one hot encoding for test set
                        if (opt_encoding == TRUE){
                          # dummify the data (One-hot-encoding)
                          dmy <- dummyVars(" ~ .", data = self$test_data)
                          test_dt_rf <- data.frame(predict(dmy, newdata = self$test_data))
                          test_dt_rf$class <- factor(test_dt_rf$class,levels = c("0","1"))
                        }else{
                          test_dt_rf <- self$test_data
                          test_dt_rf$class <- as.factor(test_dt_rf$class)
                        }
                      }
                      # Run BORUTA
                      set.seed(myseed)
                      boruta_model <- Boruta(class ~., data=train_dt_rf, doTrace=doTrace, mcAdj = T, maxRuns = maxRuns)
                      self$boruta_model <- boruta_model
                      
                      boruta_signif <- names(boruta_model$finalDecision[boruta_model$finalDecision %in% c("Confirmed","Tentative")])  # collect Confirmed and Tentative variables
                      self$boruta_importance <- boruta_signif
                    },
                    
                    run_xgb_classification = function(myseed, nrounds, ncores){
                      
                      # Create matrix: One-hot-Encoding for factor variables
                      trainm <- sparse.model.matrix(class ~ .-1,data = self$train_data)

                      train_label <- self$train_data[,"class"]
                      train_matrix <-  xgb.DMatrix(data =  as.matrix(trainm), label = train_label)
                      
                      # Same encoding for test matrix
                      testm <- sparse.model.matrix(class ~ .-1,data = self$test_data)
                      test_label <- self$test_data[,"class"]
                      test_matrix <-  xgb.DMatrix(data =  as.matrix(testm), label = test_label) 
                      
                      nc <- length(unique(train_label))
                      xgb_params <-  list("objective" = "multi:softprob",
                                          #"eval.metric"= "logloss",
                                          "num_class" = nc,
                                          # Learning parameter. Controls how much information from a new tree will be used in the Boosting
                                          # Prevents overfitting
                                          # Range is 0 to 1
                                          "eta" = 0.3,
                                          
                                          # look at only few variables to grow each new node in a tree
                                          #"colsample_bylevel" = 0.5,
                                          # subsample ratio of columns when constructing each tree. 
                                          # Subsampling occurs once for every tree constructed.
                                          # colsample_bytree
                                          
                                          # Maximum depth of the trees. Deeper trees have more terminal does and fit more data
                                          # default is 6
                                          "max_depth" = 6,
                                          
                                          # Stochastic boosting. Stochastic boosting uses only a fraction of the data to grow each tree
                                          "subsample" =1 ,
                                          
                                          #Minimum reduction in the loss function required to grow a new node in a tree.
                                          # Ranges from 0 to inf
                                          "gamma" = 0,
                                          # Controls the minimum number of observations in a terminal node
                                          "min_child_weight" = 1,
                                          nthread = ncores)
                      
                      watchlist <-  list(train =  train_matrix, test =  test_matrix)
                      
                      # eXtreme Gradient Boosting Model
                      set.seed(myseed)
                      bst_model <-  xgb.train(params = xgb_params,
                                              data = train_matrix,
                                              nrounds = 500,
                                              watchlist = watchlist,
                                              seed= myseed,verbose = F)
                      
                      self$xgb_importance <- xgb.importance(colnames(train_matrix), model = bst_model)
                      
                      # prediction
                      p <- predict(bst_model, newdata = test_matrix)
                      pred <-  matrix(p, nrow = nc, ncol = length(p)/nc)  %>%
                        t() %>%
                        data.frame() %>%
                        mutate(label = test_label,max_prob = max.col(.,"last") -1)
                      cm_xgb  <- confusionMatrix(factor(pred$max_prob), factor(pred$label), positive = "1" )
                      print(cm_xgb)
                      print(paste0("XgB Accuracy:",cm_xgb$overall[1])) 
                      self$xgb_confusion <- cm_xgb
                      
                      # Use ROCR package to plot ROC Curve
                      xgb.pred <- prediction(pred$X2, test_label)
                      xgb.perf <- performance(xgb.pred, "tpr", "fpr")
                      plot(xgb.perf)
                      roc_dt_xgb <- data.frame(FPR = unlist(xgb.perf@x.values),TPR = unlist(xgb.perf@y.values),
                                               thres =  unlist(xgb.perf@alpha.values),Method = "xGb" )
                      self$xgb_roc <-  roc_dt_xgb
                    },
                    
                    run_keras = function(myseed,center_scale_id_start,nc){
                      
                      #https://rstudio-pubs-static.s3.amazonaws.com/452498_2bb5b64288b94710a86982c3f70bb483.html
                      
                      # one hot encoding for RF
                      # dummify the data (One-hot-encoding)
                      dmy <- dummyVars(" ~ .", data = self$train_data)
                      train_dt_rf <- data.frame(predict(dmy, newdata = self$train_data))
                      train_dt_rf$class <- factor(train_dt_rf$class, levels = c("0","1"))
                      
                      # Same thing for test set
                      dmy <- dummyVars(" ~ .", data = self$test_data)
                      test_dt_rf <- data.frame(predict(dmy, newdata = self$test_data))
                      test_dt_rf$class <- factor(test_dt_rf$class,levels = c("0","1"))
                      
                      X_train <- as.matrix(train_dt_rf[,2:ncol(train_dt_rf)])
                      print(X_train)
                      X_train[,center_scale_id_start:ncol(X_train)] =  scale(X_train[,center_scale_id_start:ncol(X_train)],center = T,scale = T)
                      y_train=to_categorical(as.matrix(train_dt_rf$class))
                      y_train=as.matrix(train_dt_rf$class)
                      X_test=as.matrix(test_dt_rf[,2:ncol(test_dt_rf)])
                      X_test[,center_scale_id_start:ncol(X_test)] =  scale(X_test[,center_scale_id_start:ncol(X_test)],center = T,scale = T)
                      y_test=as.matrix(test_dt_rf$class)
                      
                      model <- keras_model_sequential() # initialize the model
                      use_session_with_seed(myseed)
                      
                      # set the models two hidden states -- this part is ad-hoc needs experimentation (heuristic) depending on the dataset
                      # model   %>%
                      #   layer_dense(units=150,
                      #               activation="relu",
                      #               input_shape=ncol(X_train),
                      #               kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = myseed)) 
                      # model %>%  
                      #   layer_dense(units=1,
                      #               activation="softmax",
                      #               kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = myseed))
                      # 
                      # model %>%
                      #   compile(loss="binary_crossentropy",
                      #           optimizer=optimizer_sgd(lr=0.001),
                      #           metrics = "accuracy")
                      
                      model %>% 
                        layer_dense(units=60,activation="relu", ncol(X_train),
                                    kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 100)) %>%
                        layer_dropout(0.5) %>%
                        layer_dense(units=40,activation="relu"
                                    ,kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 100)) %>%
                        layer_dense(units=1,activation="sigmoid",
                                    kernel_initializer=initializer_random_uniform(minval = -0.05, maxval = 0.05, seed = 100))
                      
                      compile(model,loss="binary_crossentropy",
                              optimizer=optimizer_sgd(lr=0.01,clipvalue=0.5,clipnorm=1), metrics = "accuracy")
                      
                      model
                      
                      self$keras_model <- model
                      
                      # fit the model to the training data
                      ev.model  =    fit(model,X_train, y_train, 
                                         epochs = 10, 
                                         batch_size = 10,
                                         verbose=1,
                                         validation_data = list(X_test,y_test))
                      
                      # plot loss and accuracy (modify using ggplot next stage)
                      #plot(ev.model)
                      #plot(ev.model$metrics$loss, main="Model Loss", xlab = "epoch", ylab="loss", col="orange", type="l", ylim = c(0,1))
                      #lines(ev.model$metrics$val_loss, col="skyblue")
                      #legend("topright", c("Training","Testing"), col=c("orange", "skyblue"), lty=c(1,1))
                      
                      # predict the test data (classes) 
                      evaluate(model, X_test, y_test, verbose = 1)
                      prediction_classes <- model %>% predict_classes(x = as.matrix(X_test)) %>% as.vector()
                      print(prediction_classes)
                      
                      # predict the test data (probabilities) 
                      prediction_probs <- 
                        predict_proba (object = model, 
                                       x = as.matrix(X_test)) %>%
                        as.vector()                      
                      
                      estimates_keras_tbl <- tibble(
                       truth      = as.factor(y_test),
                       estimate   = as.factor(prediction_classes),
                       class_prob = prediction_probs )

                      estimates_keras_tbl$estimate <- factor(estimates_keras_tbl$estimate, 
                                                             levels = levels(estimates_keras_tbl$truth))
                      print(estimates_keras_tbl)
                      
                      self$keras_estimates_tbl <- estimates_keras_tbl
                      
                      print(confusionMatrix(estimates_keras_tbl$estimate,
                                            estimates_keras_tbl$truth,
                                            positive="1"))
                      
                      self$keras_confusion <- 
                        confusionMatrix(estimates_keras_tbl$estimate,
                                        estimates_keras_tbl$truth,
                                        positive="1")
                      
                      estimates_keras_tbl_accuracy <- estimates_keras_tbl %>% metrics (truth, estimate)
                      print(estimates_keras_tbl_accuracy)
                      self$keras_accuracy <- estimates_keras_tbl_accuracy
                      
                      estimates_keras_tbl_roc <- estimates_keras_tbl %>% roc_auc(truth, class_prob) 
                      self$keras_roc <- estimates_keras_tbl_roc
                      
                      pfi_importance <-PFI(X_test, y_test, model)
                      self$keras_pfi_importance <- pfi_importance
                      
                      # analysis with lime
                      model_type.keras.engine.sequential.Sequential <- function(x, ...) {
                        "classification"}
                      
                      # Setup lime::predict_model()
                      predict_model.keras.engine.sequential.Sequential <- function (x, newdata, type, ...) {
                        pred <- predict(object = x, x = as.matrix(newdata))
                        data.frame (Positive = pred, Negative = 1 - pred) }
                      
                      # Usando lime() no train set
                      #explainer <- lime(x = as.data.frame(X_train), model= model)
                      
                      # Explainer
                      #system.time(
                      #  explanation <- lime::explain (
                      #    x = as.data.frame(X_test),
                      #    explainer = explainer, 
                      #    n_features = ncol(X_train),
                      #    n_labels = 1))
                      #print(explanation)
                      
                      #self$keras_lime_explanation <- explanation
                      
                    },
                    
                    greet = function() {
                      cat(paste0("Hello, my name is ", self$name, ".\n"))
                    },
                    
                    run_ordinalForest = function(myseed, ntrees, ncores, is_splitted_data, o_run_pimp, perff_type) {
                      
                      if (is_splitted_data==T){
                          test_dt_rf <- self$test_data
                          train_dt_rf <- self$train_data
                          #test_dt_rf$class <- factor(test_dt_rf$class,levels = c("0","1"))
                      }else{
                       train_dt_rf <- self$train_data
                       test_dt_rf <- self$test_data
                      }
                      # ordinal Forest
                      set.seed(myseed)
                      of <- ordfor(depvar="class", data=train_dt_rf, nsets=1000, ntreeperdiv=100, 
                                          ntreefinal=5000, perffunction = perff_type)
                      self$of_model <-of
                      
                      print(of)
                      print("Ordinal Forest training done....")
                      
                      if (o_run_pimp == TRUE) {
                        # Run permutated importance 
                        ix <- which(colnames(train_dt_rf) %in% c("class"))
                        pimp <- PIMP(train_dt_rf[,-ix], train_dt_rf$class, of$forestfinal, parallel=TRUE, ncores = ncores, seed = myseed)
                        p_test <- PimpTest(pimp)
                        pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
                        imp_dt <- pimp_all     
                        imp_dt$Predictors <- rownames(imp_dt)
                        rownames(imp_dt) <- NULL
                        self$of_importance<-imp_dt
                      }else {
                        self$of_importance <- of$varimp
                      }
                      
                      if (is_splitted_data == T){
                        # make confusion matrix
                        library(ordinalForest)
                        
                        if(perff_type == "probability"){ 
                        p_rf_test <- predict(of, test_dt_rf,type = "prob")
                        cm_rf <-  confusionMatrix(p_rf_test[["ypred"]], test_dt_rf$class,positive = "1" )
                        print(paste0("RF Accuracy :",cm_rf$overall[1])) 
                        self$of_confusion <- unlist(cm_rf$byClass)    
                        
                        f1 <- f1_score(p_rf_test$ypred, test_dt_rf$class)
                        print(f1)
                        self$of_f1 <- f1
                        
                        }
                  
                        print('Num classes > 3 - multiROC (pROC)')
                        # USE pROC for 3 or more cases
                        p_rf_test_pROC <- predict(of, test_dt_rf,type = "response")
                        roc.multi <- multiclass.roc(as.numeric(test_dt_rf$class), as.numeric(p_rf_test_pROC$ypred))
                        print(auc(roc.multi))
                        self$of_roc <- roc.multi
                      }
                      
                    },
                    
                    # run_merf_classification = function(resp.vars, rand.vars, id, rhs.vars, myseed, ntrees,
                    #                                    ncores, is_splitted_data, o_run_pimp) {
                    #   ## YET TO CONFIRM THIS WORKS
                    #   
                    #   if (is_splitted_data==T){
                    #     train_dt_rf <- self$train_data
                    #     train_dt_rf <- factor(train_dt_rf$class, levels = c("0", "1"))
                    #     test_dt_rf <- self$test_data
                    #     test_dt_rf$class <- factor(test_dt_rf$class,levels = c("0","1"))
                    #   }else{
                    #     df <- self$fulldata
                    #     df$class <- factor(df$class, levels = c("0", "1"))
                    #     ix <- sample(unique(df$id), floor(length(unique(df$id))*0.6))
                    #     test_dt_rf <- dplyr::filter(df, id %in% ix)
                    #     train_dt_rf <- dplyr::filter(df, id %notin% ix)
                    #   }
                    #   # Mixed Effect Random Forest
                    #   para <- list(
                    #     method = "cv", # internal cross-validation method for parameter tuning. See caret package
                    #     tuneLength=5, # grid size for parameter search 
                    #     number = 5,  # number of internal cross-validation
                    #     n.trees=100,   # number of trees in gbm 
                    #     ntree = 400,   # number of trees in random forest
                    #     interaction.depth=3,
                    #     shrinkage=0.05,
                    #     n.minobsinnode=10,
                    #     opt.para= TRUE, # perform parameter tuning through internal cross-validation 
                    #     include.RE = TRUE,  ## to include estimated random effect as a predictor in the machine learning model
                    #     max.iter = 10, ## maximum number of iterations for the "expectation maximization" like step  
                    #     alpha=0.05, 
                    #     minsize=5,
                    #     maxdepth=30, #normally 30
                    #     family = "binomial", 
                    #     glmer.Control = lmerControl(optimizer = "bobyqa"), #glmerControl 
                    #     likelihoodCheck = TRUE, 
                    #     nAGQ=0, 
                    #     decay = 0.05, 
                    #     K = 3, 
                    #     tol= 1e-5,
                    #     seed = myseed
                    #   )
                    #   
                    #   # Run MERF
                    #   set.seed(myseed)
                    #   merf <- MEml2(method= "MErf", data = dat.trn, id=id,  resp.vars= resp.vars, rhs.vars= rhs.vars,
                    #                 rand.vars=rand.vars, para=para, likelihoodCheck=TRUE, return.model= TRUE, verbose =3) 
                    #   self$merf_out <- merf
                    #   self$merf_model <- merf$rf.fit
                    #   
                    #   print(merf)
                    #   print("MERF training done....")
                    #   
                    #   if (o_run_pimp == TRUE) {
                    #     # Run permutated importance 
                    #     source('~/az_umms_dropbox/Model_Repo/Functions/PimpTest_RRF.R') # Had to modify pimp package to work with RRF
                    #     source('~/az_umms_dropbox/Model_Repo/Functions/PIMP_RRF.R') # LMK if you want the scripts
                    #     pimp <- PIMP_RRF(train_dt_rf[,rhs.vars], train_dt_rf[,resp.vars], merf,
                    #                      parallel = TRUE, ncores = ncores, seed = myseed)
                    #   
                    #     p_test <- PimpTest_RRF(pimp)
                    #     pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
                    #     imp_dt <- pimp_all     
                    #     imp_dt$Predictors <- rownames(imp_dt)
                    #     rownames(imp_dt) <- NULL
                    #     self$merf_importance<-imp_dt
                    #   }else {
                    #     self$merf_importance <- NULL
                    #   }
                    #   
                    #   if (is_splitted_data == T){
                    #     p_rf_test <- predict.MErf(self$merf_model, test_dt_rf, type="prob",  allow.new.levels = TRUE)
                    #     print(p_rf_test$predictions)
                    #     #pred <-  p_rf_test  %>%
                    #     #  data.frame() %>%
                    #     #  mutate(label = test_dt_rf$class,max_prob = max.col(.,"last") -1)
                    #     #res_pred = ifelse(pred[,"X1"] >= 0.5,1,0)
                    #     res_pred = ifelse(p_rf_test$predictions[,"1"] >= 0.5,1,0)
                    #     #print(res_pred)
                    #     f1_df<-data.frame(pred= res_pred, test= ifelse(test_dt_rf$class == "1",1,0))
                    #     print(f1_df)
                    #     self$merf_f1 <- f1_df
                    #  
                    #   }
                    # },
                    
                    run_rf_muliclass_classification = function(myseed, ntrees, ncores, opt_encoding, is_splitted_data, o_run_pimp) {
                      ## STILL IN WORK ##
                      # one hot encoding for train set
                      if (opt_encoding == TRUE){
                        # dummify the data (One-hot-encoding)
                        dmy <- dummyVars(" ~ .", data = self$train_data)
                        train_dt_rf <- data.frame(predict(dmy, newdata = self$train_data))
                        train_dt_rf$class <- factor(train_dt_rf$class)
                      }else{
                        train_dt_rf <- self$train_data
                        train_dt_rf$class <- factor(train_dt_rf$class)
                      }
                      
                      if (is_splitted_data==T){
                        # one hot encoding for test set
                        if (opt_encoding == TRUE){
                          # dummify the data (One-hot-encoding)
                          dmy <- dummyVars(" ~ .", data = self$test_data)
                          test_dt_rf <- data.frame(predict(dmy, newdata = self$test_data))
                          test_dt_rf$class <- factor(test_dt_rf$class)
                        }else{
                          test_dt_rf <- self$test_data
                          test_dt_rf$class <- factor(test_dt_rf$class)
                        }
                      }
                      # Random Forest
                      # Run RF multiclass classification 
                      set.seed(myseed)
                      #rf <- randomForest(class ~.,train_dt_rf, ntree = ntrees, importance = T)
                      rf <- ranger(class~.,data=train_dt_rf,importance = "permutation",
                                   probability = T,
                                   num.trees = 5000,
                                   replace=FALSE)
                      self$rfmc_model <- rf
                      
                      print(rf)
                      print("RF multiclass classification training done....")
                      
                      if (o_run_pimp == TRUE) {
                        # Run permutated importance 
                        pimp <- PIMP(train_dt_rf[,-c("class")], train_dt_rf$class, rf, parallel=TRUE, ncores = ncores, seed = myseed)
                        p_test <- PimpTest(pimp)
                        pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
                        imp_dt <- pimp_all     
                        imp_dt$Predictors <- rownames(imp_dt)
                        rownames(imp_dt) <- NULL
                        self$rfmc_importance<-imp_dt
                      }else {
                        self$rfmc_importance <- NULL
                      }
                      
                      if (is_splitted_data == T){
                        p_rf_test <- predict(self$rfmc_model, test_dt_rf)
                        print(p_rf_test$predictions)
                        res_pred = ifelse(p_rf_test$predictions[,"1"] >= 0.5,1,0)
                        #print(res_pred)
                        f1_df<-data.frame(pred= res_pred, test= ifelse(test_dt_rf$class == "1",1,0))
                        print(f1_df)
                        self$rfmc_f1 <- f1_df

                   
                        p_rf_test <- predict(of, test_dt_rf,type = "prob")
                        cm_rf <-  confusionMatrix(p_rf_test[["ypred"]], test_dt_rf$class,positive = "1" )
                        print(paste0("RF Accuracy :",cm_rf$overall[1])) 
                        self$of_confusion <- unlist(cm_rf$byClass)    
                            
                        f1 <- f1_score(p_rf_test$ypred, test_dt_rf$class)
                        print(f1)
                        self$of_f1 <- f1
                  
                        
                        # lime explanation
                        train_X <-  train_dt_rf[,!colnames(train_dt_rf) %in% "class"]
                        rf <- as_classifier(rf,labels = NULL)
                        expln <- lime(train_X, model = rf)
                        #print(expln)
                        test_X <-  test_dt_rf[,!colnames(test_dt_rf) %in% "class"]
                        lime_reasons <- lime::explain(x=test_X, explainer=expln,
                                                      n_labels = 1,
                                                      n_features = 4,
                                                      feature_select = "tree")
                        self$rfc_lime_explanation <- lime_reasons
                      }
                    }
                      
                  )
)
