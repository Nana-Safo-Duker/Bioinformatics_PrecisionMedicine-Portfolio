# Setup Guide

This guide will help you set up the Protein-Coding Potential Prediction project.

## Quick Start

### 1. Clone the Repository

```bash
git clone <repository-url>
cd "Protein-Coding Potential and Functional Annotation"
```

### 2. Python Setup

#### Option A: Using Virtual Environment (Recommended)

```bash
# Create virtual environment
python -m venv venv

# Activate virtual environment
# On Windows:
venv\Scripts\activate
# On macOS/Linux:
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

#### Option B: Using Conda

```bash
# Create conda environment
conda create -n protein_coding python=3.8
conda activate protein_coding

# Install dependencies
pip install -r requirements.txt
```

### 3. R Setup

```r
# Install R packages
source("R/requirements.R")

# Or manually install:
install.packages(c("dplyr", "caret", "randomForest", "gbm", 
                   "e1071", "pROC", "ggplot2", "gridExtra", "reshape2"))
```

### 4. Run the Analysis

#### Python

```bash
# Command line
cd src
python protein_coding_predictor.py

# Or use Jupyter Notebook
jupyter notebook notebooks/protein_coding_prediction.ipynb
```

#### R

```r
# Command line
Rscript R/protein_coding_predictor.R

# Or use RStudio/Jupyter
# Open notebooks/protein_coding_prediction_R.ipynb
```

## Directory Structure

```
Protein-Coding Potential and Functional Annotation/
├── data/                          # Data directory
│   └── genomics_data.csv         # Genomic sequence data
├── notebooks/                     # Jupyter notebooks
│   ├── protein_coding_prediction.ipynb      # Python notebook
│   └── protein_coding_prediction_R.ipynb    # R notebook
├── src/                           # Python source code
│   └── protein_coding_predictor.py
├── R/                             # R source code
│   ├── protein_coding_predictor.R
│   └── requirements.R
├── models/                        # Trained models (generated)
├── results/                       # Results and visualizations (generated)
├── requirements.txt               # Python dependencies
├── LICENSE                        # License file
├── README.md                      # Main documentation
└── SETUP.md                       # This file
```

## Troubleshooting

### Python Issues

1. **ImportError**: Make sure all dependencies are installed
   ```bash
   pip install -r requirements.txt
   ```

2. **ModuleNotFoundError**: Check that you're in the correct directory
   ```bash
   cd src
   python protein_coding_predictor.py
   ```

### R Issues

1. **Package installation fails**: Try installing packages one by one
   ```r
   install.packages("dplyr")
   install.packages("caret")
   # ... etc
   ```

2. **Library not found**: Make sure R packages are installed
   ```r
   source("R/requirements.R")
   ```

## Next Steps

1. Explore the data in `data/genomics_data.csv`
2. Run the Python or R implementation
3. Review the results in `results/` directory
4. Modify parameters in the scripts to experiment with different models
5. Check the README.md for detailed documentation

## Support

For issues or questions, please open an issue on GitHub.

