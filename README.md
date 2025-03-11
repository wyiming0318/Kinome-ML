# CancerAI_Files_2023
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
