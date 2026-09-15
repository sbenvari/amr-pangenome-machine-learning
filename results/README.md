## Results

In this analysis, we aimed to identify genes that could distinguish isolates
carrying a high number of AMR genes (≥8 AMR genes).
High-AMR-burden isolates represented approximately 15.3% of the dataset, while the remaining 1,916 isolates were classified
as lower burden.

Four models were evaluated for this purpose:

- LASSO logistic regression
- Random Forest
- Feed-forward neural network
- Regularized feed-forward neural network

Model performance was assessed on an independent held-out test set using
ROC-AUC and precision-recall AUC (PR-AUC), together with sensitivity,
specificity, precision, and accuracy at thresholds selected using the
validation set.

The strongest predictive performance was obtained with LASSO, while
Random Forest and neural-network models provided complementary nonlinear
approaches for comparison.

Despite differences between model architectures, several genomic features
were repeatedly ranked among the most informative predictors across multiple
models.

These features suggest that high AMR burden is associated with broader genomic
backgrounds enriched in mobile genetic elements, transposition and
recombination, multidrug efflux, and accessory-genome regulation.

For example, recurrent predictors included an IS256-family transposase,
a Tn3-family resolvase, an EmrE/QacC-like efflux-associated cluster, and a
TetR/AcrR-family regulatory protein.

These genes should not currently be interpreted as causing high AMR burden.
Instead, they may act as markers of plasmids, transposons, genomic islands,
or other mobile genomic environments that facilitate the acquisition,
maintenance, or co-selection of multiple resistance determinants.


## Future work

Several extensions could strengthen the biological and predictive conclusions
of this analysis.

1. **Lineage-aware validation**

   Bacterial population structure can inflate predictive performance when
   closely related isolates occur in both training and test sets. Future
   analyses will therefore evaluate model performance using lineage-aware
   splitting based on MLST, PopPUNK clusters, or another population-structure
   definition.

2. **Control for population structure**

   Lineage information could also be incorporated as a covariate to determine
   whether accessory-genome predictors provide information beyond clonal
   background.

3. **Genomic-context analysis**

   The genomic neighbourhoods of recurrent predictors will be examined to
   determine whether they occur on plasmids, transposons, integrons, or other
   mobile genetic elements and whether they are physically associated with
   AMR determinants.

4. **Functional characterization of hypothetical proteins**

   Remaining uncharacterized Panaroo clusters will be investigated using
   sequence similarity searches, conserved-domain analysis, and genomic
   neighbourhood information.

5. **Outcome sensitivity analysis**

   Alternative definitions of high AMR burden and models using AMR gene count
   directly could be evaluated to determine whether the identified genomic
   signals are robust to the choice of outcome threshold.

6. **External validation**

   The final models could be evaluated on an independent *E. coli* collection
   from a different geographic, temporal, or epidemiological source.
