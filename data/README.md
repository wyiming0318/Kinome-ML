The folowing files are included in the folder:
kinase_subtype.txt:  annotation of kinase subtypes
nn_per_subtype_pi_nopi_all.csv: mean and std of prediction performance (AUROC, AUPR, BACC) of neural network under different random seeds. For each subtype, overall performance, performance of data with privileged info (if available) and performance of data without privileged info (if available) are provided.
gbdt_per_subtype_pi_nopi_all.csv: similar to nn_per_subtype_pi_nopi_all.csv. But GBDT is deterministic under my hyperparameter once the training split is fixed. So std here are all 0s.
Ensemble_per_subtype_pi_nopi_all.csv: averaged the predictions from GBDT and neural network and re-calculated the prediction performance.


- `UsedESM2_GBDT_pred_byrep.csv`
  - Predictions from the GBDT model.
- `UsedESM2_NN_pred_byrep.csv`
  - Predictions from the neural network model.
- `UsedESM2_ensemble_pred_byrep.csv`
  - Predictions from the NN+GBDT ensemble model.
    
- `Label`: ground-truth binary label (0/1).
- `Kinase`: kinase name.
- `mutation`: mutation identifier.
- `pred_rep0` to `pred_rep4`: prediction scores for CV repetitions 0 to 4.
