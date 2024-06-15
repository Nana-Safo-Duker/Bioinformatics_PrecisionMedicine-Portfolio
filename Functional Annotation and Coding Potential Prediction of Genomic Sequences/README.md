# Protein-Coding Potential and Functional Annotation

## Overview

This project implements machine learning models to predict whether genomic regions encode proteins (protein-coding) or non-coding RNAs. The project provides comprehensive implementations in both Python and R, making it accessible to researchers using either language.

## Project Structure

```
Protein-Coding Potential and Functional Annotation/
│
├── data/                          # Data directory
│   └── genomics_data.csv         # Genomic sequence data
│
├── notebooks/                     # Jupyter notebooks
│   └── protein_coding_prediction.ipynb  # Python analysis notebook
│
├── src/                           # Python source code
│   └── protein_coding_predictor.py      # Main Python implementation
│
├── R/                             # R source code
│   ├── protein_coding_predictor.R       # Main R implementation
│   └── requirements.R                   # R package requirements
│
├── models/                        # Trained models (generated)
│   ├── protein_coding_predictor.pkl     # Python model
│   └── protein_coding_predictor_R.rds   # R model
│
├── results/                       # Results and visualizations (generated)
│   ├── performance_plots.png      # Python performance plots
│   └── performance_plots_R.png    # R performance plots
│
├── requirements.txt               # Python package requirements
├── .gitignore                     # Git ignore file
└── README.md                      # This file
```

## Dataset

### Description
The dataset contains genomic sequences with binary labels indicating whether each sequence is protein-coding (1) or non-coding RNA (0).

### Dataset License
**Important:** This dataset is provided for research and educational purposes. The original dataset license should be referenced here. If you are using a publicly available dataset, please include:

- Dataset source and citation
- License type (e.g., CC0, CC-BY, MIT, etc.)
- Any restrictions on use
- Attribution requirements

**Example citation format:**
```
If using a public dataset, please cite:
[Dataset Name]. [Dataset Source/Repository]. [License]. [URL]
```

**Note:** Please verify and update the license information based on your specific dataset source. Common genomic dataset licenses include:
- **CC0 1.0** (Public Domain Dedication)
- **CC-BY 4.0** (Creative Commons Attribution)
- **Open Data Commons Open Database License (ODbL)**
- **MIT License**

### Data Format
- **Sequences**: DNA sequences (50 nucleotides each)
- **Labels**: Binary classification (0 = non-coding RNA, 1 = protein-coding)

### Data Statistics
- Total sequences: 2,001
- Sequence length: 50 nucleotides
- Classes: 2 (non-coding RNA, protein-coding)

## Features

### Sequence Features Extracted
1. **Base Composition**: Nucleotide frequencies (A, T, G, C)
2. **GC Content**: Proportion of G and C nucleotides
3. **Dinucleotide Frequencies**: Frequency of all 16 dinucleotide combinations
4. **Start Codon Frequency**: Frequency of ATG (start codon) patterns
5. **Sequence Length**: Length of the sequence
6. **Entropy**: Sequence complexity measure
7. **AT/GC Skew**: Asymmetry in nucleotide distribution
8. **Repetitive Patterns**: Maximum repeat length

### Machine Learning Models
- **Random Forest**: Ensemble learning with decision trees
- **Gradient Boosting**: Sequential ensemble learning
- **Support Vector Machine (SVM)**: Kernel-based classification

## Installation

### Python Setup

1. **Clone the repository:**
```bash
git clone <repository-url>
cd "Protein-Coding Potential and Functional Annotation"
```

2. **Create a virtual environment (recommended):**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. **Install Python dependencies:**
```bash
pip install -r requirements.txt
```

### R Setup

1. **Install R packages:**
```r
source("R/requirements.R")
```

Or manually install:
```r
install.packages(c("dplyr", "caret", "randomForest", "gbm", 
                   "e1071", "pROC", "ggplot2", "gridExtra"))
```

2. **Load the script:**
```r
source("R/protein_coding_predictor.R")
```

## Usage

### Python

#### Command Line
```bash
cd src
python protein_coding_predictor.py
```

#### Jupyter Notebook
```bash
jupyter notebook notebooks/protein_coding_prediction.ipynb
```

#### Programmatic Usage
```python
from src.protein_coding_predictor import ProteinCodingPredictor

# Initialize predictor
predictor = ProteinCodingPredictor(model_type='random_forest')

# Load data
predictor.load_data('data/genomics_data.csv')

# Prepare features
X, y = predictor.prepare_data()

# Train model
results = predictor.train(X, y)

# Make predictions
sequences = ['GTCCACGACCGAACTCCCACCTTGACCGCAGAGGTACCACCAGAGCCCTG']
predictions, probabilities = predictor.predict(sequences)

# Save model
predictor.save_model('models/protein_coding_predictor.pkl')
```

### R

#### Command Line
```r
Rscript R/protein_coding_predictor.R
```

#### R Console
```r
source("R/protein_coding_predictor.R")
main()
```

#### Programmatic Usage
```r
# Load functions
source("R/protein_coding_predictor.R")

# Load data
df <- load_data("../data/genomics_data.csv")

# Extract features
X <- extract_features(df$Sequences)
y <- df$Labels

# Train model
results <- train_model(X, y, model_type = "randomForest")

# Make predictions
test_sequences <- c("GTCCACGACCGAACTCCCACCTTGACCGCAGAGGTACCACCAGAGCCCTG")
pred_results <- predict_sequences(results$model, test_sequences, 
                                 model_type = results$model_type)

# Save model
saveRDS(results$model, "models/protein_coding_predictor_R.rds")
```

## Model Performance

### Evaluation Metrics
- **Accuracy**: Overall correctness of predictions
- **Precision**: Proportion of positive predictions that are correct
- **Recall**: Proportion of actual positives that are correctly identified
- **F1 Score**: Harmonic mean of precision and recall
- **ROC-AUC**: Area under the receiver operating characteristic curve

### Expected Performance
Based on the implemented models, typical performance metrics:
- Accuracy: > 0.85
- Precision: > 0.80
- Recall: > 0.80
- F1 Score: > 0.80
- ROC-AUC: > 0.90

*Note: Actual performance may vary based on the specific dataset and model configuration.*

## Results Visualization

The project generates comprehensive visualizations including:
1. **Confusion Matrix**: Classification performance breakdown
2. **ROC Curve**: Receiver operating characteristic analysis
3. **Feature Importance**: Most important features for prediction
4. **Performance Metrics**: Bar chart of evaluation metrics
5. **Model Comparison**: Comparison across different models

Results are saved in the `results/` directory.

## Contributing

Contributions are welcome! Please follow these guidelines:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Citation

If you use this code in your research, please cite:

```bibtex
@software{protein_coding_predictor,
  title = {Protein-Coding Potential and Functional Annotation},
  author = {Your Name},
  year = {2025},
  url = {https://github.com/yourusername/protein-coding-predictor}
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

**Dataset License:** Please refer to the original dataset source for dataset-specific licensing information. Ensure compliance with the dataset license when using this code.

## Acknowledgments

- Thanks to the open-source community for excellent machine learning libraries
- Dataset providers and curators
- Contributors and reviewers

## Contact

For questions or issues, please open an issue on GitHub or contact:
- Email: [your-email@example.com]
- GitHub: [your-username]

## References

1. Breiman, L. (2001). Random forests. Machine learning, 45(1), 5-32.
2. Friedman, J. H. (2001). Greedy function approximation: a gradient boosting machine. Annals of statistics, 1189-1232.
3. Cortes, C., & Vapnik, V. (1995). Support-vector networks. Machine learning, 20(3), 273-297.
4. Fickett, J. W. (1982). Recognition of protein coding regions in DNA sequences. Nucleic acids research, 10(17), 5303-5318.

## Changelog

### Version 1.0.0 (2025-01-11)
- Initial release
- Python and R implementations
- Random Forest, Gradient Boosting, and SVM models
- Comprehensive feature extraction
- Visualization and evaluation tools

## Future Work

- [ ] Deep learning models (CNN, LSTM)
- [ ] Integration with biological databases
- [ ] Web interface for predictions
- [ ] Additional feature engineering
- [ ] Model interpretability analysis
- [ ] Real-time prediction API

---

**Last Updated:** January 2025

