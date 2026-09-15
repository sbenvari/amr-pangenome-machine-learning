# AMR Pangenome Machine Learning

Machine-learning analysis of antimicrobial resistance (AMR) burden in *Escherichia coli* using pangenome gene presence/absence data.

## Project overview

This project investigates whether non-AMR pangenome features can predict high AMR burden in *E. coli*.

The dataset contains a representative subset of publicly available E. coli isolates, with approximately 45,000 pangenome gene clusters generated using Panaroo. AMR determinants used to define the outcome were removed from the predictor matrix to reduce target leakage.

## Outcome definition

AMR genes were identified for each isolate using AMRFinderPlus.

For each isolate, the total number of detected AMR genes was calculated to obtain an AMR burden score.

The distribution of AMR gene counts across the isolates was then examined. Isolates carrying 8 or more AMR genes were classified as having a high AMR burden, while isolates with fewer than 8 AMR genes were classified as the lower-burden group.

This produced a binary outcome used for supervised classification:

- `0` = fewer than 8 AMR genes
- `1` = 8 or more AMR genes

AMR genes used to define this outcome were excluded from the pangenome predictor matrix to reduce target leakage.

## Models

Three supervised learning approaches were compared:

- LASSO logistic regression
- Random Forest
- Neural network

Two neural-network variants were implemented—an original model and a regularized model—resulting in four fitted models in total.

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
