import os
import pandas as pd

print("###############################################################################################################")
rerun_MUT_path = ['/media/zchen/TOSHIBAEXT/MEK sims/priority/subdomain-interaction-calc/new-energy-calc-2files']
data_dict = {}
mut_lst = []
subdir = []
# SR_interactions_gmx_PL-ACH-CL-AL
print("")
print("###############################################################################################################")
print("")
print("")

# Get path to energy.csv.
MUT_missing_data = []
for i in range(len(rerun_MUT_path)):
    for file in os.listdir(rerun_MUT_path[i]):
        print(file)
        csv_path = rerun_MUT_path[i] + '/' + file + '/' + 'coul.csv'
        LJ_csv_path = rerun_MUT_path[i] + '/' + file + '/' + 'LJ.csv'
        LJ14_csv_path = rerun_MUT_path[i] + '/' + file + '/' + 'LJ14.csv'
        print(LJ_csv_path)

        # if not os.path.exists(csv_path):
        # MUT_missing_data.append(file)
        # else:
        # Get headers.
        header_rows = list(range(23, 43, 1))
        headers = pd.read_csv(csv_path, skiprows=lambda x: x not in header_rows, header=None, names=['col'])
        start = 'SR:'
        end = "'"
        header_lst = headers['col'].values.tolist()
        for j in range(len(header_lst)):
            header_lst[j] = header_lst[j].split('"')[1]
        header_lst.insert(0, 'Time')
        # Get rows data rows.
        df = pd.read_csv(csv_path, skiprows=43, delim_whitespace=True, header=None)
        df.columns = header_lst

        LJ_header_rows = list(range(23, 33, 1))
        LJ_headers = pd.read_csv(LJ_csv_path, skiprows=lambda x: x not in LJ_header_rows, header=None, names=['col'])
        LJ_start = 'SR:'
        LJ_end = "'"
        LJ_header_lst = LJ_headers['col'].values.tolist()
        for k in range(len(LJ_header_lst)):
            LJ_header_lst[k] = LJ_header_lst[k].split('"')[1]
        LJ_header_lst.insert(0, 'Time')
        # Get rows data rows.
        LJ_df = pd.read_csv(LJ_csv_path, skiprows=33, delim_whitespace=True, header=None)
        LJ_df.columns = LJ_header_lst

        LJ14_header_rows = list(range(23, 33, 1))
        LJ14_headers = pd.read_csv(LJ14_csv_path, skiprows=lambda x: x not in LJ14_header_rows, header=None,
                                   names=['col'])
        LJ14_start = '14:'
        LJ14_end = "'"
        LJ14_header_lst = LJ14_headers['col'].values.tolist()
        for k in range(len(LJ14_header_lst)):
            LJ14_header_lst[k] = LJ14_header_lst[k].split('"')[1]
        LJ14_header_lst.insert(0, 'Time')
        # Get rows data rows.
        LJ14_df = pd.read_csv(LJ14_csv_path, skiprows=33, delim_whitespace=True, header=None)
        LJ14_df.columns = LJ14_header_lst

        masterDF1 = df.join(LJ_df.iloc[:, 1:])
        masterDF = masterDF1.join(LJ14_df.iloc[:, 1:])
        # print(masterDF)
        data_dict[file] = masterDF
        mut_lst.append(file)
# print(len(data_dict))
# print(masterDF.iloc[0:2,1:10])


for i in range(len(MUT_missing_data)):
    print('Energy calculation data for ' + MUT_missing_data[i] + ' is missing!')
print("")

# print(len(header_lst))
# for i in range(0, len(header_lst)-1, 2):
#    print(i)

# mut_lst=['I1170S']
frame_WT = len(data_dict['WT'])
# print(frame_WT)
mean_dict = {}
std_dict = {}
for mutant in mut_lst:
    frame_mut = len(data_dict[mutant])
    diff_WT_mut = (data_dict[mutant].iloc[:min(frame_WT, frame_WT), 1:] -
                   data_dict['WT'].iloc[:min(frame_WT, frame_WT), 1:])
    # print(diff_WT_mut)
    # Sum up Coul-14 and Coul-SR interactions
    diff = pd.DataFrame()
    colnames = []
    for i in range(0, len(header_lst) - 1, 2):
        sum_header = header_lst[i + 1].split(':')[1]

        # print(diff_WT_mut.iloc[:, i])
        # print(diff_WT_mut.iloc[:, i+1])
        # print(diff_WT_mut.iloc[:,  int(20+i-0.5*i)])

        sum_clmblj = diff_WT_mut.iloc[:, i] + \
                     diff_WT_mut.iloc[:, i + 1] + \
                     diff_WT_mut.iloc[:, int(20 + i - 0.5 * i)] + \
                     diff_WT_mut.iloc[:, int(30 + i - 0.5 * i)]
        sum_clmblj = sum_clmblj.tolist()
        diff[sum_header] = sum_clmblj
        colnames.append(sum_header)

    diff_mean = []
    diff_std = []
    for i in range(len(colnames)):
        interaction = colnames[i]
        diff_mean.append(diff[interaction].mean(skipna=True))
        diff_std.append(diff[interaction].std())

    mean_dict[mutant] = diff_mean
    std_dict[mutant] = diff_std

energy_diff_mean = pd.DataFrame.from_dict(mean_dict, orient='index')
energy_diff_mean.columns = colnames
energy_diff_std = pd.DataFrame.from_dict(std_dict, orient='index')
energy_diff_std.columns = colnames

print(energy_diff_mean)
print(energy_diff_std)
energy_diff_mean.to_csv('/media/zchen/TOSHIBAEXT/MEK sims/priority/subdomain-interaction-calc/new-energy-calc-2files/mean.csv')
energy_diff_std.to_csv('/media/zchen/TOSHIBAEXT/MEK sims/priority/subdomain-interaction-calc/new-energy-calc-2files/std.csv')
