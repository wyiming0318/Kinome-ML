import os
import re
import pandas as pd
import numpy as np

###################################################################
# Each kinase residue can participate in 3 hydrogen bonds at max, #
# 1 main chain donor, 1 main chain acceptor, and side chain. This #
# step removes all overcounted hydrogen bonds and output file for #
# each mutant with unique hydrogen bonds.                         #
###################################################################

dirname = input('Input the full path of the directory of output files from last step:\n')
os.chdir(dirname)

inactive_names = []

ext = '.csv'
cnt = 1
for file in os.listdir(dirname):
    print(file)
    mutID = re.split("[_.]", file)[0]
    if mutID in inactive_names:
        continue
    else:
        inactive_names.append(mutID)

print(inactive_names)
timesteps = [*range(50100, 100100, 100)]  # Process 50 frames from 50ns to 100ns with 1ns steps
#################################### WHOLE ############################################################################
for mut in inactive_names:
    df_all = pd.read_csv(mut + '_all.csv')
    all_unique = pd.DataFrame(columns=['time', 'donor_index', 'acceptor_index',
                                         'donor_resnm', 'donor_resid', 'donor_atom',
                                         'acceptor_resnm', 'acceptor_resid', 'acceptor_atom',
                                         'distance', 'angle'])

    for t in timesteps:
        seg = df_all[df_all["time"] == t]

        seg = seg.drop_duplicates(subset=['donor_resnm', 'donor_resid', 'acceptor_resnm', 'acceptor_resid'],
                                  keep='first',
                                  ignore_index=True)
        seg = seg.drop_duplicates(subset=['acceptor_resnm', 'acceptor_resid', 'acceptor_atom'],
                                  keep='first',
                                  ignore_index=True)
        seg = seg.drop_duplicates(subset=['donor_resnm', 'donor_resid', 'donor_atom'],
                                  keep='first',
                                  ignore_index=True)

        dup_donor = np.where(np.array(seg.duplicated(subset=['donor_resnm', 'donor_resid'], keep=False).tolist()) == True)[
            0].tolist()
        temp = seg.iloc[dup_donor]
        mask = temp['donor_atom'].str.len() > 2
        last_change = 0
        for i in range(len(mask) - 1):
            if mask.iloc[i + 1] and mask.iloc[i] == False and last_change == 0:
                mask.iloc[i + 1] = False
                last_change = 1
            elif mask.iloc[i + 1] and mask.iloc[i] and last_change == 0:
                mask.iloc[i] = False
                last_change = 1
            else:
                last_change = 0

        mask = mask.to_frame()
        dup_idx = mask.index[mask['donor_atom'] == True].tolist()
        seg = seg.drop(dup_idx)

        all_unique = all_unique.append(seg)

    all_unique.to_csv(mut + '_all_unique.csv', index=False, header=True)

#################################### A-Loop ############################################################################
for mut in inactive_names:
    df_aloop = pd.read_csv(mut + '_aloop.csv')
    aloop_unique = pd.DataFrame(columns=['time', 'donor_index', 'acceptor_index',
                                         'donor_resnm', 'donor_resid', 'donor_atom',
                                         'acceptor_resnm', 'acceptor_resid', 'acceptor_atom',
                                         'distance', 'angle'])

    for t in timesteps:
        seg = df_aloop[df_aloop["time"] == t]

        seg = seg.drop_duplicates(subset=['donor_resnm', 'donor_resid', 'acceptor_resnm', 'acceptor_resid'],
                                  keep='first',
                                  ignore_index=True)
        seg = seg.drop_duplicates(subset=['acceptor_resnm', 'acceptor_resid', 'acceptor_atom'],
                                  keep='first',
                                  ignore_index=True)
        seg = seg.drop_duplicates(subset=['donor_resnm', 'donor_resid', 'donor_atom'],
                                  keep='first',
                                  ignore_index=True)

        dup_donor = np.where(np.array(seg.duplicated(subset=['donor_resnm', 'donor_resid'], keep=False).tolist()) == True)[
            0].tolist()
        temp = seg.iloc[dup_donor]
        mask = temp['donor_atom'].str.len() > 2
        last_change = 0
        for i in range(len(mask) - 1):
            if mask.iloc[i + 1] and mask.iloc[i] == False and last_change == 0:
                mask.iloc[i + 1] = False
                last_change = 1
            elif mask.iloc[i + 1] and mask.iloc[i] and last_change == 0:
                mask.iloc[i] = False
                last_change = 1
            else:
                last_change = 0

        mask = mask.to_frame()
        dup_idx = mask.index[mask['donor_atom'] == True].tolist()
        seg = seg.drop(dup_idx)

        aloop_unique = aloop_unique.append(seg)

    aloop_unique.to_csv(mut + '_aloop_unique.csv', index=False, header=True)

#################################### A-Helix ###########################################################################
for mut in inactive_names:
    df_ahelix = pd.read_csv(mut + '_ahelix.csv')
    ahelix_unique = pd.DataFrame(columns=['time', 'donor_index', 'acceptor_index',
                                         'donor_resnm', 'donor_resid', 'donor_atom',
                                         'acceptor_resnm', 'acceptor_resid', 'acceptor_atom',
                                         'distance', 'angle'])

    for t in timesteps:
        seg = df_ahelix[df_ahelix["time"] == t]

        seg = seg.drop_duplicates(subset=['donor_resnm', 'donor_resid', 'acceptor_resnm', 'acceptor_resid'],
                                  keep='first',
                                  ignore_index=True)
        seg = seg.drop_duplicates(subset=['acceptor_resnm', 'acceptor_resid', 'acceptor_atom'],
                                  keep='first',
                                  ignore_index=True)
        seg = seg.drop_duplicates(subset=['donor_resnm', 'donor_resid', 'donor_atom'],
                                  keep='first',
                                  ignore_index=True)

        dup_donor = np.where(np.array(seg.duplicated(subset=['donor_resnm', 'donor_resid'], keep=False).tolist()) == True)[
            0].tolist()

        temp = seg.iloc[dup_donor]
        mask = temp['donor_atom'].str.len() > 2
        last_change = 0
        for i in range(len(mask) - 1):
            if mask.iloc[i + 1] and mask.iloc[i] == False and last_change == 0:
                mask.iloc[i + 1] = False
                last_change = 1
            elif mask.iloc[i + 1] and mask.iloc[i] and last_change == 0:
                mask.iloc[i] = False
                last_change = 1
            else:
                last_change = 0

        mask = mask.to_frame()
        dup_idx = mask.index[mask['donor_atom'] == True].tolist()
        seg = seg.drop(dup_idx)

        ahelix_unique = ahelix_unique.append(seg)

    ahelix_unique.to_csv(mut + '_ahelix_unique.csv', index=False, header=True)
#################################### A-Helix & A-loop###################################################################
for mut in inactive_names:
    df_ahelix_aloop = pd.read_csv(mut + '_ahelix_aloop.csv')
    ahelix_aloop_unique = pd.DataFrame(columns=['time', 'donor_index', 'acceptor_index',
                                         'donor_resnm', 'donor_resid', 'donor_atom',
                                         'acceptor_resnm', 'acceptor_resid', 'acceptor_atom',
                                         'distance', 'angle'])

    for t in timesteps:
        seg = df_ahelix_aloop[df_ahelix_aloop["time"] == t]

        seg = seg.drop_duplicates(subset=['donor_resnm', 'donor_resid', 'acceptor_resnm', 'acceptor_resid'],
                                  keep='first',
                                  ignore_index=True)
        seg = seg.drop_duplicates(subset=['acceptor_resnm', 'acceptor_resid', 'acceptor_atom'],
                                  keep='first',
                                  ignore_index=True)
        seg = seg.drop_duplicates(subset=['donor_resnm', 'donor_resid', 'donor_atom'],
                                  keep='first',
                                  ignore_index=True)

        dup_donor = np.where(np.array(seg.duplicated(subset=['donor_resnm', 'donor_resid'], keep=False).tolist()) == True)[
            0].tolist()
        temp = seg.iloc[dup_donor]
        mask = temp['donor_atom'].str.len() > 2
        last_change = 0
        for i in range(len(mask) - 1):
            if mask.iloc[i + 1] and mask.iloc[i] == False and last_change == 0:
                mask.iloc[i + 1] = False
                last_change = 1
            elif mask.iloc[i + 1] and mask.iloc[i] and last_change == 0:
                mask.iloc[i] = False
                last_change = 1
            else:
                last_change = 0

        mask = mask.to_frame()
        dup_idx = mask.index[mask['donor_atom'] == True].tolist()
        seg = seg.drop(dup_idx)

        ahelix_aloop_unique = ahelix_aloop_unique.append(seg)

    ahelix_aloop_unique.to_csv(mut + '_ahelix_aloop_unique.csv', index=False, header=True)

###################################################################
# This step counts frequencies of residues participating in found #
# hydrogen bonds. Occupancies of each hydrogen bonds are          #
# calculated with respect to the total number of frames.          #
###################################################################
inactive_names = []

ext = '.csv'
cnt = 1

for file in os.listdir(dirname):
    mutID = re.split("[_.]", file)[0]
    if mutID in inactive_names:
        continue
    else:
        inactive_names.append(mutID)

print(inactive_names)

ahelix_resids = list(range(72, 91, 1))
aloop_resids = list(range(173, 197, 1))
resids = ahelix_resids + aloop_resids
data = pd.DataFrame(np.zeros((len(resids), len(inactive_names))), index=resids, columns=inactive_names)
# iterate through mutant list

for mut in inactive_names:
    df_ahelix = pd.read_csv(mut + '_ahelix_unique.csv')
    df_aloop = pd.read_csv(mut + '_aloop_unique.csv')
    df_ahelix_aloop = pd.read_csv(mut + '_ahelix_aloop_unique.csv')

    lastln = df_aloop.tail(1)  # get last line of data csv
    lastrow = lastln.values.tolist()  # convert to a list
    firstln = df_aloop.head(1)
    firstrow = firstln.values.tolist()
    print(mut)
    totalframe = lastrow[0][0] / 100 - firstrow[0][0]/100  # calculate total # of frames by 1ns in data csv
    print(totalframe)

    for resid in ahelix_resids:
        count_ahelix = (df_ahelix['donor_resid'].tolist().count(resid) + df_ahelix['acceptor_resid'].tolist().count(resid) +
                        df_ahelix_aloop['donor_resid'].tolist().count(resid) +
                        df_ahelix_aloop['acceptor_resid'].tolist().count(resid))

        data.loc[resid, mut] = count_ahelix/totalframe

    for resid in aloop_resids:
        count_aloop = (df_aloop['donor_resid'].tolist().count(resid) + df_aloop['acceptor_resid'].tolist().count(resid)
                       + df_ahelix_aloop['donor_resid'].tolist().count(resid) + df_ahelix_aloop[
                           'acceptor_resid'].tolist().count(resid))
        data.loc[resid, mut] = count_aloop / totalframe


print(data)
data.to_csv('Occupancies.csv', index=False, header=True)

