import os
import MDAnalysis
import MDAnalysis.analysis.rms
import MDAnalysis.analysis.hbonds
import pandas as pd
import warnings

#########################################################################
# This step uses python package MDAnalysis to find all hydrogen bonds   #
# using input trajectories and structures. Hydrogen bonds on the        #
# αC-Helix, the A-loop, and in between these two subdoamins are found.  #
# Outputs csv files for each mutant trajectory. All hydrogen bonds are  #
# calculated using DHA-angle-cutoff = 150° and DA-distance-cutoff=3.2Å. #
#########################################################################

warnings.filterwarnings("ignore", category=DeprecationWarning)

dirname = input('Enter the full path of the directory where all structures(.gro/.tpr) and trajectories are stored:\n')
structure_fName = input('Now enter the name of the structure file (.gro/.tpr):\n')
trajectory_fName = input('Now enter the name of the trajectory file (.xtc/.trr):\n')
cnt = 1
print(trajectory_fName)
for files in os.listdir(dirname):

    subdir = dirname + '/' + files
    print(subdir)
    # print(subdir)
    if files.endswith('sh'):
        continue
    for subFile in os.listdir(subdir):
        if subFile.endswith(trajectory_fName):
            var_xtc = subdir + '/' + subFile  # trajectory file
        elif subFile.endswith(structure_fName):
            var_tpr = subdir + '/' + subFile  # structure file

    print(cnt)
    cnt=cnt+1

    print(var_xtc)
    print(var_tpr)

    u = MDAnalysis.Universe(var_tpr, var_xtc)

    ######################################### WHOLE PROTEIN ############################################################
    h_all = MDAnalysis.analysis.hbonds.hbond_analysis.HydrogenBondAnalysis(u, 'protein', 'protein',
                                                                           update_selection1=False,
                                                                           update_selection2=False,
                                                                           distance=3.2,
                                                                           angle=150,
                                                                           detect_hydrogens='distance',
                                                                           distance_type='heavy')

    h_all.run(start=501)

    h_all.generate_table()
    hbonds_occupancies_list = h_all.count_by_type().tolist()
    df2 = pd.DataFrame.from_records(h_all.table)
    df2.to_csv(files + '_all.csv', index=False, header=True)
    """
    ######################################## ahelix - aloop ############################################################
    h_ahelix_aloop = MDAnalysis.analysis.hbonds.hbond_analysis.HydrogenBondAnalysis(u, 'resid 72:90', 'resid 173:196',
                                                                                    update_selection1=False,
                                                                                    update_selection2=False,
                                                                                    distance=3.2,
                                                                                    angle=150,
                                                                                    detect_hydrogens='distance',
                                                                                    distance_type='heavy')

    h_ahelix_aloop.run(start=501)

    h_ahelix_aloop.generate_table()
    hbonds_occupancies_list = h_ahelix_aloop.count_by_type().tolist()
    df2 = pd.DataFrame.from_records(h_ahelix_aloop.table)
    df2.to_csv(files + '_ahelix_aloop.csv', index=False, header=True)

    ####################################### ahelix - ahelix ############################################################
    h_ahelix = MDAnalysis.analysis.hbonds.hbond_analysis.HydrogenBondAnalysis(u, 'resid 72:90', 'resid 72:90',
                                                                              update_selection1=False,
                                                                              update_selection2=False,
                                                                              distance=3.2,
                                                                              angle=150,
                                                                              detect_hydrogens='distance',
                                                                              distance_type='heavy')

    h_ahelix.run(start=501)

    h_ahelix.generate_table()
    hbonds_occupancies_list = h_ahelix.count_by_type().tolist()
    df2 = pd.DataFrame.from_records(h_ahelix.table)
    df2.to_csv(files + '_ahelix.csv', index=False, header=True)

    ######################################### aloop - aloop ############################################################
    h_aloop = MDAnalysis.analysis.hbonds.hbond_analysis.HydrogenBondAnalysis(u, 'resid 173:196', 'resid 173:196',
                                                                             update_selection1=False,
                                                                             update_selection2=False,
                                                                             distance=3.2,
                                                                             angle=150,
                                                                             detect_hydrogens='distance',
                                                                             distance_type='heavy')

    h_aloop.run(start=501)

    h_aloop.generate_table()
    hbonds_occupancies_list = h_aloop.count_by_type().tolist()
    df2 = pd.DataFrame.from_records(h_aloop.table)
    df2.to_csv(files + '_aloop.csv', index=False, header=True)
    """




