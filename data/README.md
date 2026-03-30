The folowing files are included in the folder:
kinase_subtype.txt:  annotation of kinase subtypes
nn_per_subtype_pi_nopi_all.csv: mean and std of prediction performance (AUROC, AUPR, BACC) of neural network under different random seeds. For each subtype, overall performance, performance of data with privileged info (if available) and performance of data without privileged info (if available) are provided.
gbdt_per_subtype_pi_nopi_all.csv: similar to nn_per_subtype_pi_nopi_all.csv. But GBDT is deterministic under my hyperparameter once the training split is fixed. So std here are all 0s.
Ensemble_per_subtype_pi_nopi_all.csv: averaged the predictions from GBDT and neural network and re-calculated the prediction performance.
