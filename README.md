# AMR Pangenome Machine Learning

Machine-learning analysis of antimicrobial resistance burden in *Escherichia coli* using pangenome gene presence/absence data.

## Project overview

This project investigates whether non-AMR pangenome features can predict high antimicrobial resistance burden in *E. coli*.

The dataset contains 2,262 isolates selected from 17,478 publicly available isolates based on core- and accessory-genome distances estimated using k-mer analysis with PopPUNK, along with approximately 45,000 pangenome gene clusters generated using Panaroo. AMR determinants used to define the outcome were removed from the predictor matrix to reduce target leakage.

## Models

Three supervised learning approaches were compared:

- LASSO logistic regression
- Random Forest
- Neural network

## Evaluation

Models were evaluated using:

- train/test splitting
- cross-validation
- ROC-AUC
- PR-AUC
- sensitivity
- specificity
- precision
- feature importance

## Biological interpretation

Top predictors shared across models were investigated further using annotation and sequence similarity searches.

Recurring signals included genes associated with:

- mobile genetic elements
- transposition and recombination
- multidrug efflux
- co-selection
- accessory-genome regulation
