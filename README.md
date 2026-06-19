# AI/ML Bioinformatics & Precision Medicine

Multi-project laboratory of clinical-grade genomics analytics created by Nana Safo-Duker. Each folder is a self-contained workflow that combines statistical genetics, explainable ML, and deep learning to solve a different precision-medicine task from mutation pathogenicity classification to transcription-factor binding prediction. This README provides the cross-project narrative by documenting structure, shared tooling, and deliverables.

---

## Table of Contents

1. [About](#about)
2. [Portfolio Overview](#portfolio-overview)
3. [Repository Layout](#repository-layout)
4. [Shared Setup Workflow](#shared-setup-workflow)
5. [Project Capsules](#project-capsules)
   - [Mutation Impact & Pathogenicity](#1-ai-driven-mutation-impact-and-pathogenicity-prediction-using-genomic-variants)
   - [DNA Biomarker Discovery](#2-dna-based-biomarker-discovery-and-computational-drug-target-identification)
   - [Genome-Wide SNP Analysis](#3-genome-wide-snp-analysis-for-genetic-variant-and-disease-association-mapping)
   - [Genomic Data Integration](#4-genomic-data-integration-and-analysis-pipeline)
   - [Regulatory Element Identification](#5-identification-and-characterization-of-regulatory-elements-in-the-human-genome)
   - [Gene Expression Modeling](#6-integrative-gene-expression-modeling-from-dna-sequence-features)
   - [TF Binding Site Prediction](#7-predictive-modeling-of-transcription-factor-binding-sites-using-deep-learning)
6. [Visualization Gallery](#visualization-gallery)
7. [Data Sources & Governance](#data-sources--governance)
8. [Contributing](#contributing)
9. [Roadmap](#roadmap)
10. [Contact](#contact)

---

## About

- **Description:** AI/ML bioinformatics portfolio covering mutation impact, biomarker discovery, SNP mapping, regulatory genomics, expression modeling, and transcription-factor binding analysis.
- **Website:** [https://nana-safo-duker.github.io/](https://nana-safo-duker.github.io/)
- **Topics:** bioinformatics · genomics · precision medicine · computational biology · biomarker discovery · gene expression · transcription factors · variant analysis · health analytics · AI in healthcare.

---

## Portfolio Overview

- **Disciplines represented:** variant pathogenicity scoring, biomarker discovery, SNP association mapping, integrative genomics, regulatory genomics, gene expression modeling, transcription factor binding analysis.
- **Languages & runtimes:** Python 3.10+, R 4.x, Jupyter Notebooks, R Markdown, Conda environments, pip, IRkernel.
- **Learning paradigms:** gradient boosting, ensemble methods, random forests, SVM, logistic regression, neural networks, CNNs/LSTMs, probabilistic analyses, statistical hypothesis testing.
- **Deliverables:** reproducible notebooks, reusable scripts/modules, serialized models, diagnostic plots, quick-start guides, and documentation files tailored to each domain problem.
- **Operational footprint:** 50+ commits across seven projects, mirrored Python/R implementations, serialized checkpoints under `models/`, and ready-to-run sample data for every workflow ([GitHub Repo](https://github.com/Nana-Safo-Duker/Bioinformatics_PrecisionMedicine-Portfolio)).

---

## Repository Layout

```
Bioinformatics_PrecisionMedicine-Portfolio/
├── <Project Folder 1>/          # End-to-end workflow (code, notebooks, docs, data placeholder)
├── <Project Folder 2>/
├── ... (7 total projects)
└── README.md                    # This portfolio guide
```

Every project follows a consistent pattern:

- `data/`: placeholder CSV with synthetic sequences; replace before production use.
- `notebooks/`: Python (`.ipynb`) and/or R (`.Rmd`) exploratory analyses.
- `scripts/` or `src/`: modularized pipelines suitable for CLI execution.
- `results/` and `models/`: generated analytics outputs and serialized estimators.
- `docs/`, `PROJECT_SUMMARY.md`, or `QUICK_START.md`: domain-specific instructions.

---

## Technology Stack & Tooling Matrix

| Layer | Tooling | Where Used |
| --- | --- | --- |
| **Languages** | Python 3.10+, R 4.x | All projects provide parallel implementations for parity and reproducibility. |
| **Environments** | `requirements.txt`, `environment.yml`, `R_requirements.R`, `install_R_packages.R` | Each folder pins its dependencies for deterministic setup. |
| **Notebooks** | Jupyter (`.ipynb`), R Markdown (`.Rmd`) | Exploratory analysis, visualization, and step-by-step walkthroughs. |
| **ML/DL Libraries** | scikit-learn, XGBoost, LightGBM, TensorFlow/Keras, caret, randomForest | Mutation prediction, biomarker discovery, TF binding models, gene expression pipelines. |
| **Visualization** | matplotlib, seaborn, ggplot2, corrplot, custom PNG exports | Result dashboards (`results/figures/`, ROC curves, correlation heatmaps, confusion matrices). |
| **Serialization** | `.pkl`, `.h5`, `.rds` models | Mutations, TF binding, biomarker discovery, and expression projects all ship trained checkpoints. |
| **Testing & QA** | Python `tests/`, R scripts, manual validation notebooks | Mutation project includes `tests/test_data_loader.py`; others rely on reproducible notebooks. |

---

## Shared Setup Workflow

1. **Clone once** to access all projects:
   ```bash
   git clone https://github.com/Nana-Safo-Duker/Bioinformatics_PrecisionMedicine-Portfolio.git
   cd Bioinformatics_PrecisionMedicine-Portfolio
   ```
2. **Select a project** from the capsules below and read its local `README.md` plus any `docs/` supplements.
3. **Create an environment**:
   - Python: `python -m venv venv` or `conda env create -f environment.yml`
   - R: run the provided `R_requirements.R`, `R_requirements.txt`, or `install_R_packages.R`
4. **Run notebooks** with Jupyter (`jupyter notebook`) or RStudio. Many R projects ship `.Rmd` notebooks—install Pandoc if you plan to knit them.
5. **Execute scripts** via CLI (`python src/main.py`, `Rscript scripts/dna_ml_pipeline.R`, etc.) for batch workflows, automated training, or integration tests.

> **Note:** All `data/genomics_data.csv` files are placeholders. Swap in IRB-approved or de-identified datasets stored securely before training or inference.

---

## Workflow Blueprint

1. **Discovery:** browse the project capsules to match a workflow to your biological question (variant scoring, biomarker discovery, TF binding, etc.).
2. **Environment provisioning:** create Python and/or R environments using the provided lock files to ensure version parity.
3. **Notebook rehearsal:** execute the notebooks end-to-end to understand feature engineering choices, metrics, and visualization outputs.
4. **Script automation:** move to the CLI scripts (`src/`, `scripts/`) for batch training, hyperparameter sweeps, or CI/CD integration.
5. **Model management:** persist trained artifacts into the `models/` directory and capture diagnostics inside `results/`.
6. **Deployment prep:** export metrics, confusion matrices, and ROC plots for downstream reporting or regulatory submissions.

---

## Project Capsules

### 1. AI-Driven Mutation Impact and Pathogenicity Prediction Using Genomic Variants
- **Problem statement:** classify missense, nonsense, and regulatory mutations as benign vs. pathogenic.
- **Highlights:** one-hot, k-mer, and numeric encoders; models include RF, SVM, logistic regression, gradient boosting, neural networks, XGBoost, and ensemble voting; Python + R parity; unit tests in `tests/`.
- **Key assets:** `src/data_loader.py|R`, `src/models.py|R`, `notebooks/mutation_prediction_python.ipynb`, `examples/quick_start.py`, pre-trained artifacts under `models/`.
- **Use cases:** variant triage, functional genomics pipelines, pathology decision support.

### 2. DNA-Based Biomarker Discovery and Computational Drug Target Identification
- **Problem statement:** engineer interpretable DNA features to pinpoint biomarker candidates and drug targets.
- **Highlights:** exhaustive feature engineering (GC/AT skew, dinucleotide/trinucleotide spectra, Shannon entropy, homopolymer runs), mirrored Python/R scripts, cross-validation, ROC/feature-importance plots, `feature_importance.csv`.
- **Key assets:** `scripts/dna_feature_extraction.*`, `scripts/dna_ml_pipeline.*`, `results/figures/`, notebooks in both languages, Conda and pip lock files.
- **Use cases:** early drug discovery, biomarker prioritization, companion diagnostic research.

### 3. Genome-Wide SNP Analysis for Genetic Variant and Disease Association Mapping
- **Problem statement:** detect SNPs/indels/structural variants and correlate them with disease phenotypes.
- **Highlights:** modular Python/R pipelines (`src/variant_identification.*`, `src/variant_analysis.*`), notebooks, correlation heatmaps, case-control statistics, troubleshooting docs for R/Pandoc upgrades.
- **Key assets:** rich `docs/` folder (setup, troubleshooting), `results/` exports (plots + CSV summaries), tests scaffold.
- **Use cases:** GWAS-style studies, cohort analytics, variant annotation dashboards.

### 4. Genomic Data Integration and Analysis Pipeline
- **Problem statement:** provide a staged workflow covering statistics, EDA, and ML for sequence classification.
- **Highlights:** four-step notebook + script suites (univariate/bivariate/multivariate, descriptive vs inferential, comprehensive EDA, ML benchmarking) in both Python and R; `environment.yml` for reproducibility.
- **Key assets:** `notebooks/python|r/1_*` through `4_*`, `scripts/python|r/`, `docs/QUICK_START.md`.
- **Use cases:** curriculum material, analytics playbooks, comparative method evaluation.

### 5. Identification and Characterization of Regulatory Elements in the Human Genome
- **Problem statement:** detect promoters, enhancers, silencers, and insulators using sequence features and motif signals.
- **Highlights:** motif-aware feature extraction (TATA box, CAAT, GC box, motif variants), ML models (RF, GBM, SVM, MLP), comprehensive visualization suite (ROC, confusion matrices, feature importance), dual-language implementations.
- **Key assets:** `src/regulatory_element_identification.py`, `scripts/regulatory_element_identification.R|Rmd`, `notebooks/regulatory_element_identification.ipynb`, `results/*.png`.
- **Use cases:** epigenomic annotation, regulatory genomics research, candidate enhancer prioritization.

### 6. Integrative Gene Expression Modeling from DNA Sequence Features
- **Problem statement:** predict high vs low gene expression directly from underlying DNA sequence context.
- **Highlights:** k-mer encoders, ML stack comprising RF, GBM, XGBoost, LightGBM, SVM; `models/` and `results/` scaffolding for exports; dual-language scripts and notebooks.
- **Key assets:** `scripts/gene_expression_prediction.py|R`, `notebooks/gene_expression_prediction.ipynb`, `docs/` references.
- **Use cases:** promoter design, synthetic biology, expression regulation studies.

### 7. Predictive Modeling of Transcription Factor Binding Sites Using Deep Learning
- **Problem statement:** locate TF binding motifs using hybrid ML/DL approaches.
- **Highlights:** CNN and LSTM architectures alongside RF, XGBoost, SVM, and R-side MLP/logistic regression; one-hot encoding pipeline; saved models (`models/*.h5`, `.pkl`); evaluation plots (training history, ROC, confusion matrices).
- **Key assets:** `scripts/python/tf_binding_prediction.py`, `scripts/r/tf_binding_prediction.R`, `notebooks/python/tf_binding_prediction.ipynb`, `results/`.
- **Use cases:** enhancer mapping, TF motif discovery, regulatory circuit modeling.

---

## Visualization Gallery

Each project ships showcase figures under `docs/images/` (committed for GitHub rendering). Full run outputs also land in that project's `results/` folder.

<table>
<tr>
<td width="50%" valign="top">

### 1. Mutation Impact & Pathogenicity

[![Mutation impact dashboard](AI-Driven%20Mutation%20Impact%20and%20Pathogenicity%20Prediction%20Using%20Genomic%20Variants/docs/images/model_comparison.png)](AI-Driven%20Mutation%20Impact%20and%20Pathogenicity%20Prediction%20Using%20Genomic%20Variants/README.md#visualizations)

</td>
<td width="50%" valign="top">

### 2. DNA Biomarker Discovery

[![DNA biomarker ROC](DNA-Based%20Biomarker%20Discovery%20and%20Computational%20Drug%20Target%20Identification/docs/images/roc_curves.png)](DNA-Based%20Biomarker%20Discovery%20and%20Computational%20Drug%20Target%20Identification/README.md#visualizations)

</td>
</tr>
<tr>
<td width="50%" valign="top">

### 3. Genome-Wide SNP Analysis

[![Variant distributions](Genome-Wide%20SNP%20Analysis%20for%20Genetic%20Variant%20and%20Disease%20Association%20Mapping/docs/images/variant_distributions.png)](Genome-Wide%20SNP%20Analysis%20for%20Genetic%20Variant%20and%20Disease%20Association%20Mapping/README.md#visualizations)

</td>
<td width="50%" valign="top">

### 4. Genomic Data Integration

[![ML model comparison](Genomic%20Data%20Integration%20and%20Analysis%20Pipeline/docs/images/ml_model_comparison.png)](Genomic%20Data%20Integration%20and%20Analysis%20Pipeline/README.md#visualizations)

</td>
</tr>
<tr>
<td width="50%" valign="top">

### 5. Regulatory Element Identification

[![Regulatory ROC](Identification%20and%20Characterization%20of%20Regulatory%20Elements%20in%20the%20Human%20Genome/docs/images/roc_curves.png)](Identification%20and%20Characterization%20of%20Regulatory%20Elements%20in%20the%20Human%20Genome/README.md#visualizations)

</td>
<td width="50%" valign="top">

### 6. Gene Expression Modeling

[![Gene expression comparison](Integrative%20Gene%20Expression%20Modeling%20from%20DNA%20Sequence%20Features/docs/images/model_comparison.png)](Integrative%20Gene%20Expression%20Modeling%20from%20DNA%20Sequence%20Features/README.md#visualizations)

</td>
</tr>
<tr>
<td width="50%" valign="top">

### 7. TF Binding Site Prediction

[![TF binding ROC](Predictive%20Modeling%20of%20Transcription%20Factor%20Binding%20Sites%20Using%20Deep%20Learning/docs/images/roc_curves.png)](Predictive%20Modeling%20of%20Transcription%20Factor%20Binding%20Sites%20Using%20Deep%20Learning/README.md#visualizations)

</td>
<td width="50%" valign="top">

Open any project README for the full set of figures (EDA panels, confusion matrices, feature importance, motif enrichment, and more).

</td>
</tr>
</table>

---

## Data Sources & Governance

- All datasets in `data/` are placeholders. Replace them with institutional data or public datasets whose licensing terms (CC-BY, CC0, GPL, etc.) allow reuse.
- Each project README highlights licensing placeholders—update those sections when integrating real data.
- For PHI/PII or clinical genomics data:
  - Store raw data outside the repo (encrypted storage, secure buckets).
  - Use `.env` or config files (added to `.gitignore`) for paths, secrets, and API keys.
  - Follow institutional review board (IRB) and GDPR/HIPAA compliance guidelines.

---

## Testing & Validation Hooks

- **Mutation Project:** `tests/test_data_loader.py` plus serialized models for regression testing.
- **Biomarker & TF Binding Pipelines:** saved confusion matrices, ROC curves, and training-history plots to compare future runs against baseline metrics.
- **SNP Analysis:** dual-language notebooks for cross-validation between Python and R outputs; correlation heatmaps provide visual regression tests.
- **Gene Expression & Regulatory Projects:** feature-importance exports and top-feature CSVs allow drift monitoring when retraining.

Integrate these hooks into CI pipelines (e.g., GitHub Actions) by executing notebooks via `papermill` or running the Python scripts with sample data.

---

## Extensibility Playbook

1. **Add a new modality:** copy the structure of an existing folder (data, notebooks, scripts, results, docs) to onboard additional omics layers (proteomics, metabolomics, single-cell RNA-seq).
2. **Share utilities:** factor reusable encoders and evaluation helpers into a `common/` package; import them via relative paths or pip-installable package for consistency.
3. **Automate metadata:** maintain a small PowerShell/Bash helper script (added at the repo root) to keep GitHub descriptions and websites in sync across sibling repositories.
4. **Document reproducibility:** include `PROJECT_SUMMARY.md` and `SETUP_GUIDE.md` files similar to the mutation and biomarker projects to describe datasets, hyperparameters, and validation cohorts.
5. **Benchmark additions:** whenever introducing a new model, add evaluation plots to `results/` and update this README’s relevant capsule to maintain reader clarity.

---

## Contributing

1. **Create a feature branch**: `git checkout -b feature/<feature-name>`.
2. **Work within a project folder** to avoid cross-contamination; add missing documentation, notebooks, or scripts there.
3. **Document everything** update both the project-level README and this portfolio README when scope changes.
4. **Validate** by executing relevant notebooks/scripts or running any included unit tests.
5. **Submit a Pull Request** detailing datasets, methods, evaluation metrics, and any compliance considerations.

---

## Roadmap

- [ ] Add a new functional annotation / coding-potential project (e.g., CPC2-style ORF assessment, RNA coding potential) with preprocessing scripts, modeling notebooks, and evaluation reports.
- [ ] Harmonize dataset license declarations and add explicit source citations across all subprojects.
- [ ] Factor shared utilities (sequence encoders, visualization helpers, model evaluators) into a reusable `common/` package.
- [ ] Add automated testing workflows (GitHub Actions) for linting, notebook smoke tests, and environment validation.
- [ ] Publish release tags bundling related analytics (e.g., mutation-centric bundle, regulatory genomics bundle).

---

## Contact

Questions, collaboration requests, or demo inquiries? Open an issue or reach out via the portfolio site: [https://nana-safo-duker.github.io/](https://nana-safo-duker.github.io/).

> Also see: [Clinical-Data-Science-Health-Analytics](https://github.com/Nana-Safo-Duker/Clinical-Data-Science-Health-Analytics) for companion clinical analytics projects.

**Last Updated**: May 2024
