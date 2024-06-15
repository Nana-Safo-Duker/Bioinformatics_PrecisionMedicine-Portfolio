# R Package Requirements
# ======================
# Install required R packages for protein-coding potential prediction
#
# Usage: source("requirements.R") to install all packages

# List of required packages
required_packages <- c(
  "dplyr",      # Data manipulation
  "caret",      # Classification and regression training
  "randomForest", # Random forest implementation
  "gbm",        # Gradient boosting machine
  "e1071",      # SVM and other utilities
  "pROC",       # ROC curve analysis
  "ggplot2",    # Data visualization
  "gridExtra",  # Grid arrangement for plots
  "reshape2"    # Data reshaping for plots
)

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE)
  }
  # Load all packages
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE)) {
      stop(paste("Failed to load package:", pkg))
    }
  }
}

# Install and load packages
cat("Installing and loading required R packages...\n")
install_if_missing(required_packages)
cat("All packages installed and loaded successfully!\n")

