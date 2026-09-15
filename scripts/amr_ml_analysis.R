
# ============================================================
# AMR PREDICTION FROM PAN-GENOME GENE PRESENCE/ABSENCE
# ============================================================
#
# This script compares supervised machine-learning approaches
# for predicting high antimicrobial resistance (AMR) burden
# in Escherichia coli using non-AMR pangenome gene
# presence/absence features.
#
# Models:
#   1. LASSO logistic regression
#   2. Random Forest
#   3. Feed-forward neural network
#   4. Regularized feed-forward neural network
#
# Evaluation:
#   - ROC-AUC
#   - PR-AUC
#   - Sensitivity
#   - Specificity
#   - Precision
#   - Accuracy
#   - Feature importance
#
# Reproducibility:
#   Random seeds are fixed where applicable.
#
# ============================================================


# ============================================================
# 0. LOAD PACKAGES
# ============================================================

library(readxl)
library(data.table)
library(glmnet)
library(Matrix)
library(pROC)
library(PRROC)
library(ranger)
library(torch)

set.seed(42)
torch_manual_seed(42)


# ============================================================
# 1. DEFINE PATHS
# ============================================================

# Expected project structure:
#
# amr-pangenome-machine-learning/
# ├── data/
# │   ├── amr.xlsx
# │   └── gene_presence_absence_no_AMR.csv
# ├── results/
# ├── scripts/
# │   └── amr_ml_analysis.R
# └── README.md

amr_file <- "data/amr.xlsx"

panaroo_file <- "data/gene_presence_absence_no_AMR.csv"

results_dir <- "results"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


# ============================================================
# 2. LOAD DATA
# ============================================================

amr <- read_excel(amr_file)

pan <- fread(
  panaroo_file
)


# ============================================================
# 3. REMOVE KNOWN AMR-ASSOCIATED FEATURES
# ============================================================

# AMR determinants used to define the outcome were removed
# upstream. This additional filter removes obvious AMR-related
# annotations that remained in the Panaroo table.

remove_amr <-
  grepl(
    "^(bla|ble|aad|ant|tet)",
    pan$Gene,
    ignore.case = TRUE
  ) |
  grepl(
    "^(bla|ble|aad|ant|tet)",
    pan$`Non-unique Gene name`,
    ignore.case = TRUE
  ) |
  grepl(
    paste(
      "aminoglycoside.*phosphotransferase",
      "aminoglycoside.*nucleotidyltransferase",
      "bleomycin resistance",
      "beta-lactamase",
      "tetracycline repressor",
      "transposon Tn10 Tet",
      sep = "|"
    ),
    pan$Annotation,
    ignore.case = TRUE
  )


# Inspect number of removed features

table(remove_amr)


# Inspect removed features

pan[
  remove_amr,
  .(
    Gene,
    `Non-unique Gene name`,
    Annotation
  )
]


# Remove them

pan <- pan[!remove_amr]


# ============================================================
# 4. CREATE GENE PRESENCE/ABSENCE MATRIX
# ============================================================

# Panaroo isolate columns contain locus tags when a gene is
# present and empty cells when absent.
#
# Convert these values into binary presence/absence:
#
# present = 1
# absent  = 0

sample_cols <- names(pan)[4:ncol(pan)]


X <- as.data.frame(
  t(
    !is.na(pan[, ..sample_cols]) &
      pan[, ..sample_cols] != ""
  )
)


# Use Panaroo cluster names as feature names

colnames(X) <- make.unique(pan$Gene)


# Use isolate names as row names

rownames(X) <- sample_cols


# Convert TRUE/FALSE into 1/0

X[] <- lapply(
  X,
  as.integer
)


# Check dimensions

dim(X)


# ============================================================
# 5. MATCH AMR OUTCOME TO ISOLATES
# ============================================================

amr2 <- amr[
  match(
    rownames(X),
    amr$Sample
  ),
]


# Confirm exact sample alignment

stopifnot(
  all(
    amr2$Sample == rownames(X)
  )
)


# Binary outcome:
#
# 0 = lower AMR burden
# 1 = high AMR burden

y <- amr2$binary


# Check class distribution

table(y)


# ============================================================
# 6. TRAIN / VALIDATION / TEST SPLIT
# ============================================================
#
# Stratified split:
#
# Training   = 70%
# Validation = 15%
# Test       = 15%
#
# Validation data are used for threshold selection.
# Test data are reserved for final evaluation.
#
# ============================================================

set.seed(42)


idx0 <- which(y == 0)
idx1 <- which(y == 1)


split_class <- function(idx) {

  idx <- sample(idx)

  n_train <- floor(
    0.70 * length(idx)
  )

  n_val <- floor(
    0.15 * length(idx)
  )

  list(

    train = idx[
      1:n_train
    ],

    val = idx[
      (n_train + 1):
        (n_train + n_val)
    ],

    test = idx[
      (n_train + n_val + 1):
        length(idx)
    ]
  )
}


split0 <- split_class(idx0)
split1 <- split_class(idx1)


train_idx <- c(
  split0$train,
  split1$train
)

val_idx <- c(
  split0$val,
  split1$val
)

test_idx <- c(
  split0$test,
  split1$test
)


y_train <- y[train_idx]
y_val   <- y[val_idx]
y_test  <- y[test_idx]


# Check outcome distribution

table(y_train)
table(y_val)
table(y_test)


# ============================================================
# 7. CREATE SPARSE MATRICES
# ============================================================

# Gene presence/absence matrices are highly sparse.
# Sparse representation reduces memory requirements.

X_sparse <- Matrix(
  as.matrix(X),
  sparse = TRUE
)


X_train <- X_sparse[
  train_idx,
]

X_val <- X_sparse[
  val_idx,
]

X_test <- X_sparse[
  test_idx,
]


# ============================================================
# 8. HELPER FUNCTIONS
# ============================================================


# ------------------------------------------------------------
# ROC-AUC and PR-AUC
# ------------------------------------------------------------

evaluate_probabilities <- function(
    y_true,
    probabilities
) {

  roc_obj <- roc(
    y_true,
    probabilities,
    quiet = TRUE
  )


  pr_obj <- pr.curve(

    scores.class0 =
      probabilities[
        y_true == 1
      ],

    scores.class1 =
      probabilities[
        y_true == 0
      ]
  )


  list(

    roc = roc_obj,

    pr = pr_obj,

    roc_auc =
      as.numeric(
        auc(roc_obj)
      ),

    pr_auc =
      pr_obj$auc.integral
  )
}



# ------------------------------------------------------------
# Performance across thresholds
# ------------------------------------------------------------

threshold_table <- function(
    y_true,
    probabilities
) {

  thresholds <- seq(
    0.05,
    0.95,
    by = 0.05
  )


  out <- lapply(
    thresholds,
    function(thr) {

      pred <- ifelse(
        probabilities >= thr,
        1,
        0
      )


      TP <- sum(
        pred == 1 &
          y_true == 1
      )

      TN <- sum(
        pred == 0 &
          y_true == 0
      )

      FP <- sum(
        pred == 1 &
          y_true == 0
      )

      FN <- sum(
        pred == 0 &
          y_true == 1
      )


      data.frame(

        threshold = thr,

        sensitivity =
          TP / (TP + FN),

        specificity =
          TN / (TN + FP),

        precision =
          ifelse(
            TP + FP == 0,
            NA,
            TP / (TP + FP)
          )
      )
    }
  )


  do.call(
    rbind,
    out
  )
}



# ------------------------------------------------------------
# Performance at one selected threshold
# ------------------------------------------------------------

evaluate_at_threshold <- function(
    y_true,
    probabilities,
    threshold
) {

  pred <- ifelse(
    probabilities >= threshold,
    1,
    0
  )


  TP <- sum(
    pred == 1 &
      y_true == 1
  )

  TN <- sum(
    pred == 0 &
      y_true == 0
  )

  FP <- sum(
    pred == 1 &
      y_true == 0
  )

  FN <- sum(
    pred == 0 &
      y_true == 1
  )


  data.frame(

    threshold = threshold,

    sensitivity =
      TP / (TP + FN),

    specificity =
      TN / (TN + FP),

    precision =
      ifelse(
        TP + FP == 0,
        NA,
        TP / (TP + FP)
      ),

    accuracy =
      (TP + TN) /
        length(y_true)
  )
}


# ============================================================
# 9. LASSO LOGISTIC REGRESSION
# ============================================================

set.seed(42)


cvfit <- cv.glmnet(

  x = X_train,

  y = y_train,

  family = "binomial",

  alpha = 1,

  nfolds = 5,

  type.measure = "deviance"
)


# Cross-validation plot

plot(cvfit)


# Optimal lambda values

cvfit$lambda.min

cvfit$lambda.1se


# ------------------------------------------------------------
# Validation predictions
# ------------------------------------------------------------

prob_lasso_val <- as.numeric(

  predict(

    cvfit,

    newx = X_val,

    s = "lambda.1se",

    type = "response"
  )
)


# Examine threshold performance

lasso_thresholds <- threshold_table(

  y_val,

  prob_lasso_val
)


print(
  lasso_thresholds
)


# ============================================================
# 10. LASSO IMPORTANT GENES
# ============================================================

coef_1se <- coef(
  cvfit,
  s = "lambda.1se"
)


selected <- which(
  coef_1se[, 1] != 0
)


selected_genes <- data.frame(

  gene =
    rownames(coef_1se)[selected],

  coefficient =
    coef_1se[selected, 1]
)


# Remove intercept

selected_genes <-
  selected_genes[
    selected_genes$gene != "(Intercept)",
  ]


# Sort by absolute coefficient magnitude

selected_genes <-
  selected_genes[
    order(
      abs(
        selected_genes$coefficient
      ),
      decreasing = TRUE
    ),
  ]


# Add Panaroo annotation

selected_annotated <- merge(

  selected_genes,

  pan[
    ,
    .(
      Gene,
      `Non-unique Gene name`,
      Annotation
    )
  ],

  by.x = "gene",

  by.y = "Gene",

  all.x = TRUE
)


selected_annotated <-
  selected_annotated[
    order(
      abs(
        selected_annotated$coefficient
      ),
      decreasing = TRUE
    ),
  ]


head(
  selected_annotated,
  30
)


write.csv(

  selected_annotated,

  file.path(
    results_dir,
    "lasso_selected_genes_annotated.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 11. RANDOM FOREST
# ============================================================

# Random Forest is more memory-intensive than LASSO.
#
# Remove extremely rare and nearly universal features:
#
# prevalence < 1%
# prevalence > 99%

freq <- colMeans(X)


keep <-
  freq >= 0.01 &
  freq <= 0.99


X_rf <- X[
  ,
  keep
]


X_train_rf <- X_rf[
  train_idx,
]

X_val_rf <- X_rf[
  val_idx,
]

X_test_rf <- X_rf[
  test_idx,
]


# Positive-class weighting to account for class imbalance

class_weight_positive <-

  sum(
    y_train == 0
  ) /

  sum(
    y_train == 1
  )


set.seed(42)


rf <- ranger(

  x = X_train_rf,

  y = factor(
    y_train,
    levels = c(0, 1)
  ),

  num.trees = 500,

  probability = TRUE,

  importance = "permutation",

  class.weights = c(

    "0" = 1,

    "1" =
      class_weight_positive
  ),

  seed = 42,

  num.threads =
    parallel::detectCores()
)


# ------------------------------------------------------------
# Validation predictions
# ------------------------------------------------------------

prob_rf_val <- predict(

  rf,

  data = X_val_rf

)$predictions[, "1"]


rf_thresholds <- threshold_table(

  y_val,

  prob_rf_val
)


print(
  rf_thresholds
)


# ============================================================
# 12. RANDOM FOREST IMPORTANT GENES
# ============================================================

rf_importance <- data.frame(

  gene =
    names(
      rf$variable.importance
    ),

  importance =
    rf$variable.importance
)


rf_importance <-
  rf_importance[
    order(
      rf_importance$importance,
      decreasing = TRUE
    ),
  ]


rf_importance_annotated <- merge(

  rf_importance,

  pan[
    ,
    .(
      Gene,
      `Non-unique Gene name`,
      Annotation
    )
  ],

  by.x = "gene",

  by.y = "Gene",

  all.x = TRUE
)


rf_importance_annotated <-
  rf_importance_annotated[
    order(
      rf_importance_annotated$importance,
      decreasing = TRUE
    ),
  ]


head(
  rf_importance_annotated,
  30
)


write.csv(

  rf_importance_annotated,

  file.path(
    results_dir,
    "random_forest_feature_importance.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 13. PREPARE TORCH DATA
# ============================================================

# Neural networks require dense numeric matrices.
#
# Note:
# Converting very large sparse matrices to dense matrices can
# require substantial RAM.

X_train_t <- torch_tensor(

  as.matrix(X_train),

  dtype = torch_float()
)


X_val_t <- torch_tensor(

  as.matrix(X_val),

  dtype = torch_float()
)


X_test_t <- torch_tensor(

  as.matrix(X_test),

  dtype = torch_float()
)


y_train_t <- torch_tensor(

  y_train,

  dtype = torch_float()

)$unsqueeze(2)


y_val_t <- torch_tensor(

  y_val,

  dtype = torch_float()

)$unsqueeze(2)


y_test_t <- torch_tensor(

  y_test,

  dtype = torch_float()

)$unsqueeze(2)


# ============================================================
# 14. ORIGINAL NEURAL NETWORK
# ============================================================

torch_manual_seed(42)


model_nn <- nn_sequential(

  nn_linear(
    ncol(X_train),
    128
  ),

  nn_relu(),

  nn_linear(
    128,
    32
  ),

  nn_relu(),

  nn_linear(
    32,
    1
  )
)


loss_fn <-
  nn_bce_with_logits_loss()


optimizer <- optim_adam(

  model_nn$parameters,

  lr = 0.001
)


epochs <- 100


for (epoch in 1:epochs) {

  model_nn$train()


  logits <-
    model_nn(
      X_train_t
    )


  loss <-
    loss_fn(
      logits,
      y_train_t
    )


  optimizer$zero_grad()

  loss$backward()

  optimizer$step()


  if (
    epoch %% 10 == 0
  ) {

    cat(

      "Original NN - Epoch:",

      epoch,

      "Loss:",

      loss$item(),

      "\n"
    )
  }
}


# ------------------------------------------------------------
# Validation predictions
# ------------------------------------------------------------

model_nn$eval()


prob_nn_val <- torch_sigmoid(

  model_nn(
    X_val_t
  )
)


prob_nn_val <- as.numeric(

  as_array(
    prob_nn_val
  )
)


nn_thresholds <- threshold_table(

  y_val,

  prob_nn_val
)


print(
  nn_thresholds
)


# ============================================================
# 15. REGULARIZED NEURAL NETWORK
# ============================================================

# Regularization:
#
# - smaller hidden layers
# - dropout
# - L2 weight decay

torch_manual_seed(42)


model_nn_reg <- nn_sequential(

  nn_linear(
    ncol(X_train),
    64
  ),

  nn_relu(),

  nn_dropout(
    p = 0.3
  ),

  nn_linear(
    64,
    16
  ),

  nn_relu(),

  nn_dropout(
    p = 0.3
  ),

  nn_linear(
    16,
    1
  )
)


loss_fn_reg <-
  nn_bce_with_logits_loss()


optimizer_reg <- optim_adam(

  model_nn_reg$parameters,

  lr = 0.001,

  weight_decay = 1e-4
)


epochs <- 100


for (epoch in 1:epochs) {

  model_nn_reg$train()


  logits <-
    model_nn_reg(
      X_train_t
    )


  loss <-
    loss_fn_reg(
      logits,
      y_train_t
    )


  optimizer_reg$zero_grad()

  loss$backward()

  optimizer_reg$step()


  if (
    epoch %% 10 == 0
  ) {

    cat(

      "Regularized NN - Epoch:",

      epoch,

      "Loss:",

      loss$item(),

      "\n"
    )
  }
}


# ------------------------------------------------------------
# Validation predictions
# ------------------------------------------------------------

model_nn_reg$eval()


prob_nn_reg_val <- torch_sigmoid(

  model_nn_reg(
    X_val_t
  )
)


prob_nn_reg_val <- as.numeric(

  as_array(
    prob_nn_reg_val
  )
)


nn_reg_thresholds <- threshold_table(

  y_val,

  prob_nn_reg_val
)


print(
  nn_reg_thresholds
)


# ============================================================
# 16. COMPARE VALIDATION THRESHOLDS
# ============================================================

cat(
  "\nLASSO\n"
)

print(
  lasso_thresholds
)


cat(
  "\nRANDOM FOREST\n"
)

print(
  rf_thresholds
)


cat(
  "\nORIGINAL NN\n"
)

print(
  nn_thresholds
)


cat(
  "\nREGULARIZED NN\n"
)

print(
  nn_reg_thresholds
)


# ============================================================
# 17. SELECT CLASSIFICATION THRESHOLDS
# ============================================================
#
# These thresholds were selected after inspecting performance
# on the validation set.
#
# Thresholds should NOT be chosen using the final test set.
#
# ============================================================

lasso_thr <- 0.20

rf_thr <- 0.25

nn_thr <- 0.05

nn_reg_thr <- 0.05


# ============================================================
# 18. FINAL TEST PREDICTIONS
# ============================================================


# ------------------------------------------------------------
# LASSO
# ------------------------------------------------------------

prob_lasso_test <- as.numeric(

  predict(

    cvfit,

    newx = X_test,

    s = "lambda.1se",

    type = "response"
  )
)


# ------------------------------------------------------------
# Random Forest
# ------------------------------------------------------------

prob_rf_test <- predict(

  rf,

  data = X_test_rf

)$predictions[, "1"]


# ------------------------------------------------------------
# Original neural network
# ------------------------------------------------------------

model_nn$eval()


prob_nn_test <- torch_sigmoid(

  model_nn(
    X_test_t
  )
)


prob_nn_test <- as.numeric(

  as_array(
    prob_nn_test
  )
)


# ------------------------------------------------------------
# Regularized neural network
# ------------------------------------------------------------

model_nn_reg$eval()


prob_nn_reg_test <- torch_sigmoid(

  model_nn_reg(
    X_test_t
  )
)


prob_nn_reg_test <- as.numeric(

  as_array(
    prob_nn_reg_test
  )
)


# ============================================================
# 19. FINAL ROC-AUC + PR-AUC
# ============================================================

lasso_final <- evaluate_probabilities(

  y_test,

  prob_lasso_test
)


rf_final <- evaluate_probabilities(

  y_test,

  prob_rf_test
)


nn_final <- evaluate_probabilities(

  y_test,

  prob_nn_test
)


nn_reg_final <- evaluate_probabilities(

  y_test,

  prob_nn_reg_test
)


final_auc <- data.frame(

  Model = c(

    "LASSO",

    "Random Forest",

    "Original Neural Network",

    "Regularized Neural Network"
  ),


  ROC_AUC = c(

    lasso_final$roc_auc,

    rf_final$roc_auc,

    nn_final$roc_auc,

    nn_reg_final$roc_auc
  ),


  PR_AUC = c(

    lasso_final$pr_auc,

    rf_final$pr_auc,

    nn_final$pr_auc,

    nn_reg_final$pr_auc
  )
)


print(
  final_auc
)


# ============================================================
# 20. FINAL THRESHOLD-BASED TEST PERFORMANCE
# ============================================================

lasso_test_metrics <- evaluate_at_threshold(

  y_test,

  prob_lasso_test,

  lasso_thr
)


rf_test_metrics <- evaluate_at_threshold(

  y_test,

  prob_rf_test,

  rf_thr
)


nn_test_metrics <- evaluate_at_threshold(

  y_test,

  prob_nn_test,

  nn_thr
)


nn_reg_test_metrics <- evaluate_at_threshold(

  y_test,

  prob_nn_reg_test,

  nn_reg_thr
)


final_threshold_results <- rbind(

  data.frame(
    Model = "LASSO",
    lasso_test_metrics
  ),

  data.frame(
    Model = "Random Forest",
    rf_test_metrics
  ),

  data.frame(
    Model = "Original NN",
    nn_test_metrics
  ),

  data.frame(
    Model = "Regularized NN",
    nn_reg_test_metrics
  )
)


print(
  final_threshold_results
)


# ============================================================
# 21. NEURAL NETWORK PERMUTATION IMPORTANCE
# ============================================================
#
# Uses:
#   Regularized neural network
#   Validation dataset
#
# Importance is measured as:
#
# baseline PR-AUC - PR-AUC after permuting one feature
#
# WARNING:
# This is computationally expensive because each gene is
# permuted and evaluated independently.
#
# ============================================================

nn_permutation_importance <- function(

    model,

    X_matrix,

    y_true,

    gene_names
) {

  model$eval()


  X_original_t <- torch_tensor(

    X_matrix,

    dtype = torch_float()
  )


  pred_original <- torch_sigmoid(

    model(
      X_original_t
    )
  )


  pred_original <- as.numeric(

    as_array(
      pred_original
    )
  )


  baseline_pr <- pr.curve(

    scores.class0 =
      pred_original[
        y_true == 1
      ],

    scores.class1 =
      pred_original[
        y_true == 0
      ]

  )$auc.integral


  importance <- numeric(
    ncol(X_matrix)
  )


  for (
    j in seq_len(
      ncol(X_matrix)
    )
  ) {

    X_perm <- X_matrix


    # Break the association between one gene and the outcome

    X_perm[, j] <-
      sample(
        X_perm[, j]
      )


    X_perm_t <- torch_tensor(

      X_perm,

      dtype = torch_float()
    )


    pred_perm <- torch_sigmoid(

      model(
        X_perm_t
      )
    )


    pred_perm <- as.numeric(

      as_array(
        pred_perm
      )
    )


    perm_pr <- pr.curve(

      scores.class0 =
        pred_perm[
          y_true == 1
        ],

      scores.class1 =
        pred_perm[
          y_true == 0
        ]

    )$auc.integral


    importance[j] <-
      baseline_pr -
      perm_pr


    if (
      j %% 100 == 0
    ) {

      cat(

        "Permutation importance:",

        j,

        "/",

        ncol(X_matrix),

        "genes\n"
      )
    }
  }


  results <- data.frame(

    gene = gene_names,

    importance = importance
  )


  results[
    order(
      results$importance,
      decreasing = TRUE
    ),
  ]
}


# ------------------------------------------------------------
# Run neural-network permutation importance
# ------------------------------------------------------------

X_val_matrix <- as.matrix(
  X_val
)


nn_importance <- nn_permutation_importance(

  model = model_nn_reg,

  X_matrix = X_val_matrix,

  y_true = y_val,

  gene_names =
    colnames(
      X_val_matrix
    )
)


head(
  nn_importance,
  30
)


# ============================================================
# 22. ANNOTATE NEURAL NETWORK IMPORTANT GENES
# ============================================================

nn_importance_annotated <- merge(

  nn_importance,

  pan[
    ,
    .(
      Gene,
      `Non-unique Gene name`,
      Annotation
    )
  ],

  by.x = "gene",

  by.y = "Gene",

  all.x = TRUE
)


nn_importance_annotated <-
  nn_importance_annotated[
    order(
      nn_importance_annotated$importance,
      decreasing = TRUE
    ),
  ]


head(
  nn_importance_annotated,
  30
)


write.csv(

  nn_importance_annotated,

  file.path(
    results_dir,
    "neural_network_gene_importance.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 23. COMPARE TOP GENES ACROSS MODELS
# ============================================================

top_n <- 30


# ------------------------------------------------------------
# Top LASSO genes
# ------------------------------------------------------------

top_lasso <- head(

  selected_annotated$gene,

  top_n
)


# ------------------------------------------------------------
# Top Random Forest genes
# ------------------------------------------------------------

top_rf <- head(

  rf_importance$gene,

  top_n
)


# ------------------------------------------------------------
# Top neural-network genes
# ------------------------------------------------------------

top_nn <- head(

  nn_importance$gene,

  top_n
)


# ------------------------------------------------------------
# Pairwise overlaps
# ------------------------------------------------------------

cat(
  "\nLASSO vs Random Forest:\n"
)


print(

  intersect(

    top_lasso,

    top_rf
  )
)


cat(
  "\nLASSO vs Neural Network:\n"
)


print(

  intersect(

    top_lasso,

    top_nn
  )
)


cat(
  "\nRandom Forest vs Neural Network:\n"
)


print(

  intersect(

    top_rf,

    top_nn
  )
)


# ------------------------------------------------------------
# Features shared by all three model classes
# ------------------------------------------------------------

cat(
  "\nGenes identified by ALL THREE:\n"
)


common_all <- Reduce(

  intersect,

  list(

    top_lasso,

    top_rf,

    top_nn
  )
)


print(
  common_all
)


# ============================================================
# 24. SAVE FINAL RESULTS
# ============================================================

write.csv(

  final_auc,

  file.path(
    results_dir,
    "model_auc_comparison.csv"
  ),

  row.names = FALSE
)


write.csv(

  final_threshold_results,

  file.path(
    results_dir,
    "model_threshold_comparison.csv"
  ),

  row.names = FALSE
)


write.csv(

  data.frame(
    gene = common_all
  ),

  file.path(
    results_dir,
    "common_top_genes.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 25. SAVE MODEL OBJECTS
# ============================================================
#
# Optional, but useful for reproducibility so models do not
# need to be retrained every time.
#
# ============================================================

saveRDS(

  cvfit,

  file.path(
    results_dir,
    "lasso_model.rds"
  )
)


saveRDS(

  rf,

  file.path(
    results_dir,
    "random_forest_model.rds"
  )
)


# ============================================================
# END OF ANALYSIS
# ============================================================
