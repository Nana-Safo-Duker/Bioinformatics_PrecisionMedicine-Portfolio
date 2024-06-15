# Protein-Coding Potential Predictor
# ===================================
#
# This script implements machine learning models to predict whether genomic regions
# encode proteins (protein-coding) or non-coding RNAs.
#
# Author: Bioinformatics Research Team
# Date: 2025

# Load required libraries
library(dplyr)
library(caret)
library(randomForest)
library(gbm)
library(e1071)
library(pROC)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Set random seed for reproducibility
set.seed(42)

# Feature extraction function
# Extract comprehensive features from DNA sequences.
# Parameters: sequences - vector of DNA sequences
# Returns: feature_matrix - data frame with extracted features
extract_features <- function(sequences) {
  features_list <- list()
  
  for (seq in sequences) {
    seq <- toupper(seq)
    seq_char <- strsplit(seq, "")[[1]]
    seq_length <- nchar(seq)
    
    # Base composition (nucleotide frequencies)
    A_count <- sum(seq_char == "A")
    T_count <- sum(seq_char == "T")
    G_count <- sum(seq_char == "G")
    C_count <- sum(seq_char == "C")
    
    A_freq <- A_count / seq_length
    T_freq <- T_count / seq_length
    G_freq <- G_count / seq_length
    C_freq <- C_count / seq_length
    
    # GC content
    GC_content <- (G_count + C_count) / seq_length
    
    # Dinucleotide frequencies
    dinucleotides <- c("AA", "AT", "AG", "AC", "TA", "TT", "TG", "TC",
                      "GA", "GT", "GG", "GC", "CA", "CT", "CG", "CC")
    dinuc_freqs <- sapply(dinucleotides, function(dinuc) {
      pattern <- gregexpr(dinuc, seq, fixed = TRUE)[[1]]
      if (pattern[1] == -1) return(0)
      length(pattern) / (seq_length - 1)
    })
    names(dinuc_freqs) <- paste0("dinuc_", names(dinuc_freqs))
    
    # Start codon frequency (ATG)
    ATG_count <- length(gregexpr("ATG", seq, fixed = TRUE)[[1]])
    if (ATG_count == 1 && gregexpr("ATG", seq, fixed = TRUE)[[1]][1] == -1) {
      ATG_freq <- 0
    } else {
      ATG_freq <- ATG_count / (seq_length - 2)
    }
    
    # Sequence length
    seq_len <- seq_length
    
    # Entropy (sequence complexity)
    base_counts <- c(A_count, T_count, G_count, C_count)
    base_probs <- base_counts / seq_length
    base_probs <- base_probs[base_probs > 0]
    entropy <- -sum(base_probs * log2(base_probs))
    
    # AT/GC skew
    AT_skew <- ifelse((A_count + T_count) > 0, 
                     (A_count - T_count) / (A_count + T_count), 0)
    GC_skew <- ifelse((G_count + C_count) > 0,
                     (G_count - C_count) / (G_count + C_count), 0)
    
    # Combine all features
    features <- c(
      A_freq = A_freq,
      T_freq = T_freq,
      G_freq = G_freq,
      C_freq = C_freq,
      GC_content = GC_content,
      dinuc_freqs,
      ATG_freq = ATG_freq,
      seq_length = seq_len,
      entropy = entropy,
      AT_skew = AT_skew,
      GC_skew = GC_skew
    )
    
    features_list[[length(features_list) + 1]] <- features
  }
  
  # Convert to data frame
  feature_matrix <- do.call(rbind, features_list)
  feature_matrix <- as.data.frame(feature_matrix)
  
  return(feature_matrix)
}

# Load data
# Load data from CSV file.
# Parameters: file_path - path to CSV file
# Returns: data - data frame with sequences and labels
load_data <- function(file_path) {
  cat("Loading data from", file_path, "...\n")
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  cat("Data loaded:", nrow(data), "sequences\n")
  cat("Class distribution:\n")
  print(table(data$Labels))
  return(data)
}

# Train model function
# Train a machine learning model.
# Parameters: X - feature matrix, y - labels, model_type - type of model, test_size - proportion of test set
# Returns: results - list with model and evaluation metrics
train_model <- function(X, y, model_type = "randomForest", test_size = 0.2) {
  # Convert labels to factor
  y <- as.factor(y)
  levels(y) <- c("NonCoding", "ProteinCoding")
  
  # Split data
  train_index <- createDataPartition(y, p = 1 - test_size, list = FALSE)
  X_train <- X[train_index, ]
  X_test <- X[-train_index, ]
  y_train <- y[train_index]
  y_test <- y[-train_index]
  
  cat("\nTraining set size:", nrow(X_train), "\n")
  cat("Test set size:", nrow(X_test), "\n")
  
  # Train model
  if (model_type == "randomForest") {
    model <- randomForest(
      x = X_train,
      y = y_train,
      ntree = 100,
      mtry = sqrt(ncol(X_train)),
      importance = TRUE
    )
  } else if (model_type == "gbm") {
    model <- gbm(
      y ~ .,
      data = cbind(y = y_train, X_train),
      distribution = "bernoulli",
      n.trees = 100,
      interaction.depth = 5,
      shrinkage = 0.1,
      verbose = FALSE
    )
  } else if (model_type == "svm") {
    model <- svm(
      x = X_train,
      y = y_train,
      kernel = "radial",
      probability = TRUE
    )
  } else {
    stop("Unknown model type: ", model_type)
  }
  
  # Predictions
  if (model_type == "gbm") {
    y_pred_proba <- predict(model, newdata = X_test, n.trees = 100, type = "response")
    y_pred <- ifelse(y_pred_proba > 0.5, "ProteinCoding", "NonCoding")
    y_pred <- factor(y_pred, levels = levels(y_test))
  } else {
    y_pred <- predict(model, newdata = X_test)
    if (model_type == "svm") {
      y_pred_proba <- attr(predict(model, newdata = X_test, probability = TRUE), "probabilities")[, "ProteinCoding"]
    } else {
      y_pred_proba <- predict(model, newdata = X_test, type = "prob")[, "ProteinCoding"]
    }
  }
  
  # Calculate metrics
  cm <- confusionMatrix(y_pred, y_test)
  accuracy <- cm$overall["Accuracy"]
  precision <- cm$byClass["Precision"]
  recall <- cm$byClass["Recall"]
  f1 <- cm$byClass["F1"]
  
  # ROC-AUC
  roc_obj <- roc(as.numeric(y_test == "ProteinCoding"), y_pred_proba)
  roc_auc <- auc(roc_obj)
  
  # Print results
  cat("\n", "=", rep("=", 48), "\n", sep = "")
  cat("MODEL PERFORMANCE METRICS\n")
  cat("=", rep("=", 48), "\n", sep = "")
  cat("Accuracy: ", sprintf("%.4f", accuracy), "\n", sep = "")
  cat("Precision:", sprintf("%.4f", precision), "\n", sep = "")
  cat("Recall:   ", sprintf("%.4f", recall), "\n", sep = "")
  cat("F1 Score: ", sprintf("%.4f", f1), "\n", sep = "")
  cat("ROC-AUC:  ", sprintf("%.4f", roc_auc), "\n", sep = "")
  cat("=", rep("=", 48), "\n", sep = "")
  
  cat("\nConfusion Matrix:\n")
  print(cm$table)
  
  cat("\nClassification Report:\n")
  print(cm)
  
  # Store confusion matrix as matrix for plotting
  cm_matrix <- as.matrix(cm$table)
  
  # Cross-validation
  cat("\nPerforming cross-validation...\n")
  cv_folds <- createFolds(y_train, k = 5, returnTrain = TRUE)
  cv_scores <- numeric(length(cv_folds))
  
  for (i in 1:length(cv_folds)) {
    cv_train_idx <- cv_folds[[i]]
    cv_X_train <- X_train[cv_train_idx, ]
    cv_X_val <- X_train[-cv_train_idx, ]
    cv_y_train <- y_train[cv_train_idx]
    cv_y_val <- y_train[-cv_train_idx]
    
    if (model_type == "randomForest") {
      cv_model <- randomForest(x = cv_X_train, y = cv_y_train, ntree = 50)
      cv_pred <- predict(cv_model, newdata = cv_X_val)
    } else if (model_type == "gbm") {
      cv_model <- gbm(y ~ ., data = cbind(y = cv_y_train, cv_X_train),
                     distribution = "bernoulli", n.trees = 50, verbose = FALSE)
      cv_pred <- predict(cv_model, newdata = cv_X_val, n.trees = 50, type = "response")
      cv_pred <- ifelse(cv_pred > 0.5, "ProteinCoding", "NonCoding")
      cv_pred <- factor(cv_pred, levels = levels(cv_y_val))
    } else {
      cv_model <- svm(x = cv_X_train, y = cv_y_train)
      cv_pred <- predict(cv_model, newdata = cv_X_val)
    }
    
    cv_scores[i] <- mean(cv_pred == cv_y_val)
  }
  
  cat("CV Accuracy:", sprintf("%.4f", mean(cv_scores)), 
      "(+/-", sprintf("%.4f", sd(cv_scores) * 2), ")\n")
  
  return(list(
    model = model,
    model_type = model_type,
    X_train = X_train,
    X_test = X_test,
    y_train = y_train,
    y_test = y_test,
    y_pred = y_pred,
    y_pred_proba = y_pred_proba,
    metrics = list(
      accuracy = as.numeric(accuracy),
      precision = as.numeric(precision),
      recall = as.numeric(recall),
      f1 = as.numeric(f1),
      roc_auc = as.numeric(roc_auc),
      cv_accuracy = mean(cv_scores),
      cv_std = sd(cv_scores)
    ),
    confusion_matrix = cm_matrix,
    roc_obj = roc_obj
  ))
}

# Plot results
# Plot visualization of results.
# Parameters: results - results list from train_model(), save_path - path to save the plot
plot_results <- function(results, save_path = NULL) {
  # Create plots
  # Confusion Matrix
  cm_melted <- reshape2::melt(results$confusion_matrix)
  colnames(cm_melted) <- c("Reference", "Prediction", "Freq")
  p1 <- ggplot(cm_melted, aes(x = Prediction, y = Reference, fill = Freq)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = "Confusion Matrix", x = "Predicted Label", y = "True Label") +
    theme_minimal() +
    geom_text(aes(label = Freq), color = "black", size = 5)
  
  # ROC Curve
  roc_df <- data.frame(
    FPR = 1 - results$roc_obj$specificities,
    TPR = results$roc_obj$sensitivities
  )
  p2 <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_line(color = "darkorange", size = 1.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "navy") +
    labs(title = paste("ROC Curve (AUC =", sprintf("%.2f", results$metrics$roc_auc), ")"),
         x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal() +
    coord_fixed()
  
  # Metrics bar chart
  metrics_df <- data.frame(
    Metric = c("Accuracy", "Precision", "Recall", "F1 Score", "ROC-AUC"),
    Score = c(results$metrics$accuracy, results$metrics$precision,
             results$metrics$recall, results$metrics$f1, results$metrics$roc_auc)
  )
  p3 <- ggplot(metrics_df, aes(x = Metric, y = Score, fill = Metric)) +
    geom_bar(stat = "identity") +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Model Performance Metrics", y = "Score") +
    theme_minimal() +
    theme(legend.position = "none") +
    ylim(0, 1) +
    geom_text(aes(label = sprintf("%.3f", Score)), vjust = -0.5)
  
  # Feature importance (if available)
  if (results$model_type == "randomForest" && "importance" %in% names(results$model)) {
    imp_df <- data.frame(
      Feature = rownames(results$model$importance),
      Importance = results$model$importance[, "MeanDecreaseGini"]
    )
    imp_df <- imp_df[order(imp_df$Importance, decreasing = TRUE), ][1:15, ]
    p4 <- ggplot(imp_df, aes(x = reorder(Feature, Importance), y = Importance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = "Top 15 Feature Importances", x = "Feature", y = "Importance") +
      theme_minimal()
  } else {
    p4 <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "Feature importance\nnot available\nfor this model") +
      theme_void()
  }
  
  # Arrange plots
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  if (!is.null(save_path)) {
    ggsave(save_path, width = 15, height = 12, dpi = 300)
    cat("Plot saved to", save_path, "\n")
  }
}

# Predict function
# Predict protein-coding potential for new sequences.
# Parameters: model - trained model, sequences - vector of DNA sequences, model_type - type of model
# Returns: list with predictions and probabilities
predict_sequences <- function(model, sequences, model_type = "randomForest") {
  X <- extract_features(sequences)
  
  if (model_type == "gbm") {
    pred_proba <- predict(model, newdata = X, n.trees = 100, type = "response")
    pred <- ifelse(pred_proba > 0.5, "ProteinCoding", "NonCoding")
  } else if (model_type == "svm") {
    pred <- predict(model, newdata = X)
    pred_proba <- attr(predict(model, newdata = X, probability = TRUE), "probabilities")[, "ProteinCoding"]
  } else {
    pred <- predict(model, newdata = X)
    pred_proba <- predict(model, newdata = X, type = "prob")[, "ProteinCoding"]
  }
  
  return(list(predictions = pred, probabilities = pred_proba))
}

# Main function
main <- function() {
  # Load data
  data_path <- "../data/genomics_data.csv"
  df <- load_data(data_path)
  
  # Extract features
  cat("\nExtracting features from sequences...\n")
  sequences <- df$Sequences
  labels <- df$Labels
  X <- extract_features(sequences)
  y <- labels
  
  cat("Feature matrix shape:", nrow(X), "x", ncol(X), "\n")
  cat("Number of features:", ncol(X), "\n")
  
  # Train model
  results <- train_model(X, y, model_type = "randomForest")
  
  # Plot results
  plot_results(results, save_path = "../results/performance_plots_R.png")
  
  # Save model
  saveRDS(results$model, "../models/protein_coding_predictor_R.rds")
  cat("\nModel saved to ../models/protein_coding_predictor_R.rds\n")
  
  # Example predictions
  cat("\n", "=", rep("=", 48), "\n", sep = "")
  cat("EXAMPLE PREDICTIONS\n")
  cat("=", rep("=", 48), "\n", sep = "")
  test_sequences <- c(
    "GTCCACGACCGAACTCCCACCTTGACCGCAGAGGTACCACCAGAGCCCTG",
    "GAGTTTATATGGCGCGAGCCTAGTGGTTTTTGTACTTGTTTGTCGCGTCG"
  )
  pred_results <- predict_sequences(results$model, test_sequences, 
                                   model_type = results$model_type)
  
  for (i in 1:length(test_sequences)) {
    cat("\nSequence", i, ":", substr(test_sequences[i], 1, 30), "...\n")
    cat("Prediction:", as.character(pred_results$predictions[i]), "\n")
    cat("Probability:", sprintf("%.4f", pred_results$probabilities[i]), "\n")
  }
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}

