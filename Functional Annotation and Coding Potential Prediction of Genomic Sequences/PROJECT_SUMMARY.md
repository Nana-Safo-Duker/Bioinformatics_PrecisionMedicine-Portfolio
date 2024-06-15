# Project Summary

## Overview
This project implements machine learning models to predict protein-coding potential of genomic sequences. It distinguishes between protein-coding regions and non-coding RNA regions using sequence features and various machine learning algorithms.

## Project Structure

```
Protein-Coding Potential and Functional Annotation/
├── data/                          # Dataset directory
│   └── genomics_data.csv         # 2,001 genomic sequences with labels
├── notebooks/                     # Jupyter notebooks
│   ├── protein_coding_prediction.ipynb      # Python analysis notebook
│   └── protein_coding_prediction_R.ipynb    # R analysis notebook
├── src/                           # Python source code
│   └── protein_coding_predictor.py         # Main Python implementation
├── R/                             # R source code
│   ├── protein_coding_predictor.R          # Main R implementation
│   └── requirements.R                      # R package installer
├── models/                        # Trained models (generated after training)
├── results/                       # Results and visualizations (generated)
├── requirements.txt               # Python dependencies
├── LICENSE                        # MIT License
├── README.md                      # Comprehensive documentation
├── SETUP.md                       # Setup instructions
└── .gitignore                     # Git ignore rules
```

## Key Features

### 1. Feature Extraction
- Base composition (A, T, G, C frequencies)
- GC content
- Dinucleotide frequencies (16 combinations)
- Start codon (ATG) frequency
- Sequence length
- Entropy (sequence complexity)
- AT/GC skew
- Repetitive pattern analysis

### 2. Machine Learning Models
- **Random Forest**: Ensemble learning with 100 trees
- **Gradient Boosting**: Sequential ensemble learning
- **Support Vector Machine (SVM)**: Kernel-based classification

### 3. Evaluation Metrics
- Accuracy
- Precision
- Recall
- F1 Score
- ROC-AUC
- Cross-validation scores
- Confusion matrix

### 4. Visualizations
- Confusion matrix heatmap
- ROC curve
- Feature importance plots
- Performance metrics bar charts
- Model comparison plots

## Implementation Details

### Python Implementation
- **Main Script**: `src/protein_coding_predictor.py`
- **Notebook**: `notebooks/protein_coding_prediction.ipynb`
- **Dependencies**: pandas, numpy, scikit-learn, matplotlib, seaborn, joblib
- **Features**: 
  - Object-oriented design
  - Comprehensive feature extraction
  - Multiple model support
  - Model saving/loading
  - Visualization tools

### R Implementation
- **Main Script**: `R/protein_coding_predictor.R`
- **Notebook**: `notebooks/protein_coding_prediction_R.ipynb`
- **Dependencies**: dplyr, caret, randomForest, gbm, e1071, pROC, ggplot2, gridExtra, reshape2
- **Features**:
  - Functional programming approach
  - Comprehensive feature extraction
  - Multiple model support
  - Model saving/loading
  - Visualization tools

## Dataset Information

- **Total Sequences**: 2,001
- **Sequence Length**: 50 nucleotides
- **Classes**: 
  - 0: Non-coding RNA
  - 1: Protein-coding
- **Format**: CSV with 'Sequences' and 'Labels' columns

## Usage

### Quick Start (Python)
```bash
cd src
python protein_coding_predictor.py
```

### Quick Start (R)
```r
source("R/protein_coding_predictor.R")
main()
```

### Jupyter Notebooks
```bash
# Python
jupyter notebook notebooks/protein_coding_prediction.ipynb

# R (requires IRkernel)
jupyter notebook notebooks/protein_coding_prediction_R.ipynb
```

## Expected Results

Based on the implemented models:
- **Accuracy**: > 0.85
- **Precision**: > 0.80
- **Recall**: > 0.80
- **F1 Score**: > 0.80
- **ROC-AUC**: > 0.90

*Note: Actual performance may vary based on dataset and model configuration.*

## License

- **Code**: MIT License (see LICENSE file)
- **Dataset**: Please refer to original dataset source for dataset-specific licensing information

## Next Steps

1. Run the analysis using Python or R implementations
2. Experiment with different model parameters
3. Try additional feature engineering
4. Implement deep learning models (CNN, LSTM)
5. Integrate with biological databases
6. Develop web interface for predictions

## Contributing

Contributions are welcome! Please follow the guidelines in README.md.

## Contact

For questions or issues, please open an issue on GitHub.

---

**Project Status**: ✅ Complete and Ready for Use
**Last Updated**: January 2025

