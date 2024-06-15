"""
Protein-Coding Potential Predictor
===================================

This script implements machine learning models to predict whether genomic regions
encode proteins (protein-coding) or non-coding RNAs.

Author: Bioinformatics Research Team
Date: 2025
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.metrics import (accuracy_score, precision_score, recall_score, 
                             f1_score, roc_auc_score, confusion_matrix, 
                             classification_report, roc_curve, auc)
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
import os
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Set random seed for reproducibility
np.random.seed(42)

class SequenceFeatureExtractor:
    """Extract features from DNA sequences for machine learning."""
    
    def __init__(self):
        self.bases = ['A', 'T', 'G', 'C']
        
    def extract_features(self, sequences):
        """
        Extract comprehensive features from DNA sequences.
        
        Parameters:
        -----------
        sequences : list or array-like
            List of DNA sequences
            
        Returns:
        --------
        features : numpy array
            Feature matrix
        """
        features = []
        
        for seq in sequences:
            seq = seq.upper()
            seq_features = []
            
            # 1. Base composition (nucleotide frequencies)
            total_length = len(seq)
            for base in self.bases:
                seq_features.append(seq.count(base) / total_length)
            
            # 2. GC content
            gc_count = seq.count('G') + seq.count('C')
            seq_features.append(gc_count / total_length)
            
            # 3. Dinucleotide frequencies
            dinucleotides = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC',
                           'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']
            for dinuc in dinucleotides:
                seq_features.append(seq.count(dinuc) / (total_length - 1))
            
            # 4. Codon usage patterns (for protein-coding regions)
            # Look for start codons (ATG)
            seq_features.append(seq.count('ATG') / (total_length - 2))
            
            # 5. Sequence length
            seq_features.append(total_length)
            
            # 6. Entropy (sequence complexity)
            base_counts = [seq.count(base) for base in self.bases]
            base_probs = [count / total_length if total_length > 0 else 0 
                         for count in base_counts]
            entropy = -sum(p * np.log2(p) if p > 0 else 0 for p in base_probs)
            seq_features.append(entropy)
            
            # 7. Repetitive sequences (simple repeats)
            max_repeat = self._max_repeat_length(seq)
            seq_features.append(max_repeat)
            
            # 8. AT/GC skew
            a_count = seq.count('A')
            t_count = seq.count('T')
            g_count = seq.count('G')
            c_count = seq.count('C')
            
            at_skew = (a_count - t_count) / (a_count + t_count) if (a_count + t_count) > 0 else 0
            gc_skew = (g_count - c_count) / (g_count + c_count) if (g_count + c_count) > 0 else 0
            seq_features.extend([at_skew, gc_skew])
            
            features.append(seq_features)
        
        return np.array(features)
    
    def _max_repeat_length(self, seq):
        """Calculate maximum repeat length in sequence."""
        max_len = 1
        for i in range(len(seq)):
            for j in range(i + 1, len(seq)):
                k = 0
                while j + k < len(seq) and seq[i + k] == seq[j + k]:
                    k += 1
                max_len = max(max_len, k)
        return max_len


class ProteinCodingPredictor:
    """Main class for protein-coding potential prediction."""
    
    def __init__(self, model_type='random_forest'):
        """
        Initialize the predictor.
        
        Parameters:
        -----------
        model_type : str
            Type of model to use ('random_forest', 'gradient_boosting', 'svm')
        """
        self.model_type = model_type
        self.feature_extractor = SequenceFeatureExtractor()
        self.model = None
        self.feature_names = None
        
        if model_type == 'random_forest':
            self.model = RandomForestClassifier(
                n_estimators=100,
                max_depth=20,
                min_samples_split=5,
                min_samples_leaf=2,
                random_state=42,
                n_jobs=-1
            )
        elif model_type == 'gradient_boosting':
            self.model = GradientBoostingClassifier(
                n_estimators=100,
                learning_rate=0.1,
                max_depth=5,
                random_state=42
            )
        elif model_type == 'svm':
            self.model = SVC(
                kernel='rbf',
                probability=True,
                random_state=42,
                C=1.0,
                gamma='scale'
            )
        else:
            raise ValueError(f"Unknown model type: {model_type}")
    
    def load_data(self, file_path):
        """
        Load data from CSV file.
        
        Parameters:
        -----------
        file_path : str
            Path to the CSV file
        """
        print(f"Loading data from {file_path}...")
        self.df = pd.read_csv(file_path)
        print(f"Data loaded: {len(self.df)} sequences")
        print(f"Class distribution:\n{self.df['Labels'].value_counts()}")
        return self.df
    
    def prepare_data(self):
        """Prepare features and labels from the dataset."""
        print("Extracting features from sequences...")
        sequences = self.df['Sequences'].values
        labels = self.df['Labels'].values
        
        # Extract features
        X = self.feature_extractor.extract_features(sequences)
        y = labels
        
        print(f"Feature matrix shape: {X.shape}")
        print(f"Number of features: {X.shape[1]}")
        
        return X, y
    
    def train(self, X, y, test_size=0.2, random_state=42):
        """
        Train the model.
        
        Parameters:
        -----------
        X : numpy array
            Feature matrix
        y : numpy array
            Labels
        test_size : float
            Proportion of test set
        random_state : int
            Random seed
        """
        print(f"\nTraining {self.model_type} model...")
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=test_size, random_state=random_state, stratify=y
        )
        
        print(f"Training set size: {X_train.shape[0]}")
        print(f"Test set size: {X_test.shape[0]}")
        
        # Train model
        self.model.fit(X_train, y_train)
        
        # Evaluate on test set
        y_pred = self.model.predict(X_test)
        y_pred_proba = self.model.predict_proba(X_test)[:, 1]
        
        # Calculate metrics
        accuracy = accuracy_score(y_test, y_pred)
        precision = precision_score(y_test, y_pred)
        recall = recall_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred)
        roc_auc = roc_auc_score(y_test, y_pred_proba)
        
        print("\n" + "="*50)
        print("MODEL PERFORMANCE METRICS")
        print("="*50)
        print(f"Accuracy:  {accuracy:.4f}")
        print(f"Precision: {precision:.4f}")
        print(f"Recall:    {recall:.4f}")
        print(f"F1 Score:  {f1:.4f}")
        print(f"ROC-AUC:   {roc_auc:.4f}")
        print("="*50)
        
        # Confusion matrix
        cm = confusion_matrix(y_test, y_pred)
        print("\nConfusion Matrix:")
        print(cm)
        
        # Classification report
        print("\nClassification Report:")
        print(classification_report(y_test, y_pred, 
                                   target_names=['Non-coding', 'Protein-coding']))
        
        # Cross-validation
        print("\nPerforming cross-validation...")
        cv_scores = cross_val_score(self.model, X_train, y_train, 
                                    cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
                                    scoring='accuracy')
        print(f"CV Accuracy: {cv_scores.mean():.4f} (+/- {cv_scores.std() * 2:.4f})")
        
        return {
            'X_train': X_train,
            'X_test': X_test,
            'y_train': y_train,
            'y_test': y_test,
            'y_pred': y_pred,
            'y_pred_proba': y_pred_proba,
            'metrics': {
                'accuracy': accuracy,
                'precision': precision,
                'recall': recall,
                'f1': f1,
                'roc_auc': roc_auc,
                'cv_accuracy': cv_scores.mean(),
                'cv_std': cv_scores.std()
            },
            'confusion_matrix': cm
        }
    
    def predict(self, sequences):
        """
        Predict protein-coding potential for new sequences.
        
        Parameters:
        -----------
        sequences : list or array-like
            List of DNA sequences
            
        Returns:
        --------
        predictions : numpy array
            Binary predictions (0 = non-coding, 1 = protein-coding)
        probabilities : numpy array
            Prediction probabilities
        """
        X = self.feature_extractor.extract_features(sequences)
        predictions = self.model.predict(X)
        probabilities = self.model.predict_proba(X)[:, 1]
        return predictions, probabilities
    
    def save_model(self, file_path):
        """Save the trained model to disk."""
        os.makedirs(os.path.dirname(file_path) if os.path.dirname(file_path) else '.', exist_ok=True)
        joblib.dump(self.model, file_path)
        print(f"Model saved to {file_path}")
    
    def load_model(self, file_path):
        """Load a trained model from disk."""
        self.model = joblib.load(file_path)
        print(f"Model loaded from {file_path}")
    
    def plot_results(self, results, save_path=None):
        """
        Plot visualization of results.
        
        Parameters:
        -----------
        results : dict
            Results dictionary from train() method
        save_path : str
            Path to save the plot
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Confusion Matrix
        cm = results['confusion_matrix']
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[0, 0])
        axes[0, 0].set_title('Confusion Matrix')
        axes[0, 0].set_ylabel('True Label')
        axes[0, 0].set_xlabel('Predicted Label')
        axes[0, 0].set_xticklabels(['Non-coding', 'Protein-coding'])
        axes[0, 0].set_yticklabels(['Non-coding', 'Protein-coding'])
        
        # ROC Curve
        fpr, tpr, _ = roc_curve(results['y_test'], results['y_pred_proba'])
        roc_auc = auc(fpr, tpr)
        axes[0, 1].plot(fpr, tpr, color='darkorange', lw=2, 
                       label=f'ROC curve (AUC = {roc_auc:.2f})')
        axes[0, 1].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        axes[0, 1].set_xlim([0.0, 1.0])
        axes[0, 1].set_ylim([0.0, 1.05])
        axes[0, 1].set_xlabel('False Positive Rate')
        axes[0, 1].set_ylabel('True Positive Rate')
        axes[0, 1].set_title('ROC Curve')
        axes[0, 1].legend(loc="lower right")
        axes[0, 1].grid(True)
        
        # Feature Importance (if available)
        if hasattr(self.model, 'feature_importances_'):
            importances = self.model.feature_importances_
            indices = np.argsort(importances)[::-1][:15]  # Top 15 features
            axes[1, 0].barh(range(len(indices)), importances[indices])
            axes[1, 0].set_yticks(range(len(indices)))
            axes[1, 0].set_yticklabels([f'Feature {i}' for i in indices])
            axes[1, 0].set_xlabel('Importance')
            axes[1, 0].set_title('Top 15 Feature Importances')
            axes[1, 0].invert_yaxis()
        else:
            axes[1, 0].text(0.5, 0.5, 'Feature importance\nnot available\nfor this model',
                           ha='center', va='center', transform=axes[1, 0].transAxes)
            axes[1, 0].set_title('Feature Importance')
        
        # Metrics Bar Chart
        metrics = results['metrics']
        metric_names = ['Accuracy', 'Precision', 'Recall', 'F1 Score', 'ROC-AUC']
        metric_values = [metrics['accuracy'], metrics['precision'], 
                        metrics['recall'], metrics['f1'], metrics['roc_auc']]
        axes[1, 1].bar(metric_names, metric_values, color=['#1f77b4', '#ff7f0e', 
                                                          '#2ca02c', '#d62728', '#9467bd'])
        axes[1, 1].set_ylabel('Score')
        axes[1, 1].set_title('Model Performance Metrics')
        axes[1, 1].set_ylim([0, 1])
        axes[1, 1].grid(True, axis='y', alpha=0.3)
        for i, v in enumerate(metric_values):
            axes[1, 1].text(i, v + 0.02, f'{v:.3f}', ha='center', va='bottom')
        
        plt.tight_layout()
        
        if save_path:
            os.makedirs(os.path.dirname(save_path) if os.path.dirname(save_path) else '.', exist_ok=True)
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Plot saved to {save_path}")
        
        plt.show()


def main():
    """Main function to run the prediction pipeline."""
    # Initialize predictor
    predictor = ProteinCodingPredictor(model_type='random_forest')
    
    # Load data
    data_path = '../data/genomics_data.csv'
    predictor.load_data(data_path)
    
    # Prepare data
    X, y = predictor.prepare_data()
    
    # Train model
    results = predictor.train(X, y)
    
    # Plot results
    predictor.plot_results(results, save_path='../results/performance_plots.png')
    
    # Save model
    predictor.save_model('../models/protein_coding_predictor.pkl')
    
    # Example prediction on new sequences
    print("\n" + "="*50)
    print("EXAMPLE PREDICTIONS")
    print("="*50)
    test_sequences = [
        'GTCCACGACCGAACTCCCACCTTGACCGCAGAGGTACCACCAGAGCCCTG',  # Known protein-coding
        'GAGTTTATATGGCGCGAGCCTAGTGGTTTTTGTACTTGTTTGTCGCGTCG',  # Known non-coding
    ]
    predictions, probabilities = predictor.predict(test_sequences)
    
    for i, (seq, pred, prob) in enumerate(zip(test_sequences, predictions, probabilities)):
        label = 'Protein-coding' if pred == 1 else 'Non-coding'
        print(f"\nSequence {i+1}: {seq[:30]}...")
        print(f"Prediction: {label}")
        print(f"Probability: {prob:.4f}")


if __name__ == '__main__':
    main()

