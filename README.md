# Kinome-ML_Files_2023
README: Supporting data for 2023 manuscript Wang et. al. The support files consist of two folders in which a MD folder includes the files to setup Molecular Dynamics simulations and post-processing scripts for feature extraction, as well as a Docking folder with input structure files for molecular docking and resulting poses of docked ternary complexes.
## Folders and Subfolders
### MD: 
This folder has the input .mdp scripts requried to perform Molecular Dynamics simulations (MD_input_mdp), intial input structures of active and inactive conformations of selected ALK, BRAF, EGFR, HER2, and MEK1 mutants (MD_initial_structures), and scripts used for feature extraction from MD trajectories (MD_Feature_Calculations) <br> <br>
#### MD_input_mdp:<br>
1. minim.mdp/md.mdp/md2.mdp - .mdp files required to perform Molecular Dynamics of mutant systems <br>
2. energy_calc - .mdp file for the rerun of trajectories for calculatiosn of subdomain interactions <br>
#### MD_initial_structures:<br>
Homology modeled structure of kinase mutants after engineered mutations and filling of missing residues from crystal structures
#### MD_Feature_Calculations: <br>
1. Hydrogen_Bond_Occupancies - .py scripts to perform hydrogen bond occupancy calculations using the MDAnalysis 1.0.1 package. The scripts and package are run in the Python 2.7 environments <br>
2. RMSF_SASA_RGYR - .sh scripts for calculations of RMSF, SASA, and radius of gyration of selected subdoamins P-loop, αC-helix, C-loop, and A-loop on the mutants using GROMACS <br>
3. Subdomain_Interactions - .sh scripts to calculate interaction potentials between the selected subdoamins <br>
### Docking:
This folder contains the input structures of selected kinase mutants, ATP, and peptide substrates used in molecular docking of the ternary complexes to model phosphorylation transition states. Docking simulations were performed using the Schrodinger 2019.4 suite.
 1. LastFramePDBs - .pdb files of 215 kinase mutant structures extracted from the last frame of MD trajectories. These structures were uses as receptor in docking simulations
 2. Peptide_Substrates - Four .mae structures of peptide substrates of five kinases. Peptide PEYIN was used for both EGFR and HER2 mutant systems
 3. ATP_qtot-4.mae - Structure of ATP<sup>4-</sup>
 4. TopPoses - .pdb files of top ranked poses of Kinase-ATP-peptide complexes for each mutant system. Poses were ranked using Glide Score and only poses that met the biological validitiy criteria were kept for ranking

### Machine learning:
classicML.py contains cross validation code to train and evaluate classic ML methods (GBDT, nGBDT and SVM).
NNmodel.py contains the neural network model.
NNtrainer.py contains cross validation code to train and evaluate neural networks.
The code was tested under 

numpy 1.26.4

sklearn 1.4.2

torch 2.6.0+cu124

xgboost 2.1.2

ngboost 0.5.1

### Data:
Raw data and processed data for ML can be downloaded from https://drive.google.com/file/d/1NCTE-2R0GyTNT-6PdBrAAqBrJ4ybY1oh/view?usp=sharing. To run ML code, unzip the file and put the data to the local folder.

./data/ESM2_first_token_emb_train.npy: CLS token embeddings of wild type kinase domain sequences extracted by ESM-2 (training data for PCA);

./data/ESM2_mutation_token_emb_train.npy: Token embeddings of all mutation sites from kinase domains extracted by ESM-2 (training data for PCA);

./data/data.csv: Curated kinase mutation dataset. Each row is a kinase mutation sample with feature vector, activation status and ID provided;

./processed_data/dataset.pkl: a pickle file containing data extracted from ./data/data.csv. Our ML code directly loads this file for model training and evaluation. Specifically, it is a list containing (1) feature matrix of kinase mutations, (2) a binary vector indicating the activation status of kinase mutations, (3) a vector indicating the confidence/fidelity of the activation status labels, (4) a binary vector indicating whether a kinase mutation has privileged features, (5) numpy array of ./data/ESM2_first_token_emb_train.npy (used for dimensionality reduction purpose), (6) numpy array of ./data/ESM2_mutation_token_emb_train.npy (used for dimensionality reduction purpose), (7) a list containing the dimension size of different types of features used in this study (i.e., ESM_first_dim, ESM_mutated_dim, group_code_dim, basic_feat_dim, MD_feat_dim, dock_feat_dim).

The [data/independent-validation](data/independent-validation) subfolder contains scored-variant lists and complete annotated prediction tables for FGFR, HER2 and ERK2/MAPK1, including experimental labels, model outputs, and SIFT/PolyPhen and hotspot-based evaluation scenarios. Its [README](data/independent-validation/README.md) provides dataset counts, column definitions, and the distinction between baseline predictions and retrospective hotspot-assisted results.

### ESM_feature_generation:

- Retained `wildtype_kinase_data.csv` and `mutant_df.csv` provide original reference sequences, domain coordinates and 1,003 experimental mutation records; they are not the complete training feature matrix.
- `01_KD_ESM_Feature_Generation.ipynb` validates sequence/coordinate identities and generates KD-only ESM-2 layer-33 mutant BOS and mutation-residue vectors (1,280 dimensions each), with out-of-KD rows explicitly unscorable.
- `02_Training_Data_and_PCA_Workflow.ipynb` documents the unchanged original training-data layout and two unstandardized reference PCAs: **50 components each, 100 combined**. Inference reuses model-associated saved preprocessing; classical no-PI models use 167 inputs, while released NN checkpoints use their recorded input dimensions (100 in the current released bundle).
- See [the folder README](ESM_feature_generation/README.md) for dependencies, optional execution settings and historical-provenance limitations. The superseded exploratory notebook and four inconsistent embedding/PCA exports were removed; the source CSVs and license remain unchanged. Newly generated embeddings do not replace historical training features.
