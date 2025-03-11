import numpy as np
import pandas as pd
import csv
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import balanced_accuracy_score as accuracy_score
import copy
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import math, copy, time
from torch.autograd import Variable
import pickle
from NNmodel import *
from sklearn.decomposition import PCA
from transformers import get_linear_schedule_with_warmup,get_cosine_schedule_with_warmup
from torchvision.ops import focal_loss, sigmoid_focal_loss

#Load dataset
f = open('./processed_data/dataset.pkl', 'rb')
dataset = pickle.load(f)
f.close()

data_feat, label, fidelity_vec, PI_vec, ESM_first_basis, ESM_mutated_basis, feat_dim_list = dataset
ESM_first_dim, ESM_mutated_dim, group_code_dim, basic_feat_dim, MD_feat_dim, dock_feat_dim = feat_dim_list

#PCA for ESM2 feats
PCA_dim = 50
PCA_first = PCA(n_components=PCA_dim, random_state=18)
PCA_mutated = PCA(n_components=PCA_dim, random_state=18)
PCA_first.fit(ESM_first_basis)
PCA_mutated.fit(ESM_mutated_basis)


#Replace original ESM2 feats with dimensionality-reduced ones
ESM_first_PCA = PCA_first.transform(data_feat[:, :ESM_first_dim])
ESM_mutated_PCA = PCA_mutated.transform(data_feat[:, ESM_first_dim:ESM_first_dim + ESM_mutated_dim])

X = copy.deepcopy(np.hstack([ESM_first_PCA, ESM_mutated_PCA, data_feat[:, ESM_first_dim + ESM_mutated_dim:]]))

PI_dim = MD_feat_dim + dock_feat_dim

use_basic = False
use_basic_joint = False

if use_basic:
	main_dim = PCA_dim + PCA_dim
else:
	main_dim = PCA_dim + PCA_dim + group_code_dim + basic_feat_dim


#add basic feats to x_joint
if use_basic_joint:
	X = np.hstack([X, X[:, :PCA_dim + PCA_dim + group_code_dim + basic_feat_dim]])
	joint_dim = PI_dim + PCA_dim + PCA_dim + group_code_dim + basic_feat_dim
else:
	X = np.hstack([X, X[:, :PCA_dim + PCA_dim]])
	joint_dim = PI_dim + PCA_dim + PCA_dim

joint_start_idx = PCA_dim + PCA_dim + group_code_dim + basic_feat_dim



def trainer(train_data, test_data, h_list_main, h_list_joint1, h_list_joint2, dropout, n_tree, tree_depth, tree_feature_rate, gamma, num_epoch):
	
	X_train, y_train, PI_train = train_data

	cutoff = 0.5 

	X_train_withPI = []
	y_train_withPI = []
	X_train_withoutPI = []
	y_train_withoutPI = []

	for i in range(len(PI_train)):
		if PI_train[i] == 1:
			X_train_withPI.append(X_train[i])
			y_train_withPI.append(y_train[i])
		else:
			X_train_withoutPI.append(X_train[i])
			y_train_withoutPI.append(y_train[i])


	X_main_mem = copy.deepcopy(np.array(X_train_withPI))[:, :main_dim]
	X_joint_mem = copy.deepcopy(np.array(X_train_withPI))[:, joint_start_idx:]

	X_train_withPI_main = np.array(X_train_withPI)[:, :main_dim]
	y_train_withPI = np.array(y_train_withPI).reshape(-1, 1)
	X_train_withoutPI_main = np.array(X_train_withoutPI)[:, :main_dim]
	y_train_withoutPI = np.array(y_train_withoutPI).reshape(-1, 1)


	model = MemNN(main_dim, joint_dim, h_list_main, h_list_joint1, h_list_joint2, dropout, n_tree, tree_depth, tree_feature_rate)
	pytorch_total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
	print ('Model initialized, total num params:', pytorch_total_params)
	model.cuda()
	optimizer = torch.optim.AdamW(model.parameters(), lr=0.001, betas=(0.9, 0.98), eps=1e-8, weight_decay=5*(1e-2))
	scheduler = get_linear_schedule_with_warmup(
				optimizer,
				num_warmup_steps = 50,
				num_training_steps = num_epoch)



	for i in range(num_epoch):
		model.train()
		optimizer.zero_grad()


		X_main_mem_gpu = torch.FloatTensor(X_main_mem).cuda()
		X_joint_mem_gpu = torch.FloatTensor(X_joint_mem).cuda()

		X_train_withPI_main_gpu = torch.FloatTensor(X_train_withPI_main).cuda()
		y_train_withPI_gpu = torch.FloatTensor(y_train_withPI).cuda()
		X_train_withoutPI_main_gpu = torch.FloatTensor(X_train_withoutPI_main).cuda()
		y_train_withoutPI_gpu = torch.FloatTensor(y_train_withoutPI).cuda()


		y_pred_withPI = model.forward_main(X_train_withPI_main_gpu, [True, None, X_joint_mem_gpu])
		y_pred_withoutPI = model.forward_main(X_train_withoutPI_main_gpu, [False, X_main_mem_gpu, X_joint_mem_gpu])
		y_pred_aux = model.forward_main(X_train_withPI_main_gpu, [False, X_main_mem_gpu, X_joint_mem_gpu])


		all_pred = torch.cat([y_pred_withPI, y_pred_withoutPI, y_pred_aux], dim=0)
		all_true = torch.cat([y_train_withPI_gpu, y_train_withoutPI_gpu, y_train_withPI_gpu], dim=0)

		#all_pred = torch.cat([y_pred_withoutPI, y_pred_aux], dim=0)
		#all_true = torch.cat([y_train_withoutPI_gpu, y_train_withPI_gpu], dim=0)

		loss = sigmoid_focal_loss(all_pred, all_true, reduction='mean', gamma=gamma)

		loss.backward()
		nn.utils.clip_grad_norm_(model.parameters(), 1)
		optimizer.step()
		scheduler.step()


		if (i+1) % 100 == 0:
			PI_result_pack_train, NO_PI_result_pack_train, ALL_result_pack_train, _, _  = eval_func(model, train_data, [X_main_mem, X_joint_mem], cutoff)
			PI_result_pack_test, NO_PI_result_pack_test, ALL_result_pack_test, PI_raw_test, NO_PI_raw_test = eval_func(model, test_data, [X_main_mem, X_joint_mem], cutoff)
			print ('Epoch', i+1)
			print ('With PI results:')
			print ('Train aupr, auc, bacc:', PI_result_pack_train[0], PI_result_pack_train[1], PI_result_pack_train[2])
			print ('Test aupr, auc, bacc:', PI_result_pack_test[0], PI_result_pack_test[1], PI_result_pack_test[2])
			print ('Without PI results:')
			print ('Train aupr, auc, bacc:', NO_PI_result_pack_train[0], NO_PI_result_pack_train[1], NO_PI_result_pack_train[2])
			print ('Test aupr, auc, bacc:', NO_PI_result_pack_test[0], NO_PI_result_pack_test[1], NO_PI_result_pack_test[2])
			print ('All results:')
			print ('Train aupr, auc, bacc:', ALL_result_pack_train[0], ALL_result_pack_train[1], ALL_result_pack_train[2])
			print ('Test aupr, auc, bacc:', ALL_result_pack_test[0], ALL_result_pack_test[1], ALL_result_pack_test[2])
			print ('=====================================')


	return PI_raw_test, NO_PI_raw_test

def eval_func(model, data, memory, cutoff):

	model.eval()

	X, y, PI_vec = data

	X_withPI = []
	y_withPI = []
	X_withoutPI = []
	y_withoutPI = []

	for i in range(len(PI_vec)):
		if PI_vec[i] == 1:
			X_withPI.append(X[i])
			y_withPI.append(y[i])
		else:
			X_withoutPI.append(X[i])
			y_withoutPI.append(y[i])

	X_withPI_main = np.array(X_withPI)[:, :main_dim]
	X_withPI_joint = np.array(X_withPI)[:, joint_start_idx:]

	y_withPI = np.array(y_withPI).reshape(-1, 1)
	X_withoutPI_main = np.array(X_withoutPI)[:, :main_dim]
	y_withoutPI = np.array(y_withoutPI).reshape(-1, 1)

	X_main_mem_gpu = torch.FloatTensor(memory[0]).cuda()
	X_joint_mem_gpu = torch.FloatTensor(memory[1]).cuda()

	X_withPI_main_gpu = torch.FloatTensor(X_withPI_main).cuda()
	X_withPI_joint_gpu = torch.FloatTensor(X_withPI_joint).cuda()
	X_withoutPI_main_gpu = torch.FloatTensor(X_withoutPI_main).cuda()


	y_pred_withoutPI = torch.sigmoid(model.forward_main(X_withoutPI_main_gpu, [False, X_main_mem_gpu, X_joint_mem_gpu])).cpu().detach().numpy()
	y_pred_withPI = torch.sigmoid(model.forward_main(X_withPI_main_gpu, [True, None, X_withPI_joint_gpu])).cpu().detach().numpy()

	aupr_withPI, auc_withPI, acc_withPI = get_metrics_out(y_pred_withPI.reshape(-1), y_withPI.reshape(-1), cutoff)
	aupr_withoutPI, auc_withoutPI, acc_withoutPI = get_metrics_out(y_pred_withoutPI.reshape(-1), y_withoutPI.reshape(-1), cutoff)

	y_pred_all = y_pred_withPI.reshape(-1).tolist() + y_pred_withoutPI.reshape(-1).tolist()
	y_all = y_withPI.reshape(-1).tolist() + y_withoutPI.reshape(-1).tolist()
	aupr_all, auc_all, acc_all = get_metrics_out(np.array(y_pred_all).reshape(-1), np.array(y_all).reshape(-1), cutoff)

	return [aupr_withPI, auc_withPI, acc_withPI], [aupr_withoutPI, auc_withoutPI, acc_withoutPI], [aupr_all, auc_all, acc_all], [y_pred_withPI, y_withPI], [y_pred_withoutPI, y_withoutPI]


def get_metrics_out(pred, y, cutoff):
	cutoff_array = np.ones(len(y))*cutoff
	auc = roc_auc_score(y, pred)
	aupr = average_precision_score(y, pred)
	y_pred_binary = (pred >= cutoff_array).astype(float)
	acc = accuracy_score(y, y_pred_binary.reshape(-1))
	return aupr, auc, acc


def CV_NN(data, method_param, num_fold, rs):

	X, y, PI_ind, feat_dim_list = data
	ESM_first_dim, ESM_mutated_dim, group_code_dim, basic_feat_dim, MD_feat_dim, dock_feat_dim = feat_dim_list

	h_list_main, h_list_joint1, h_list_joint2, dropout, n_tree, tree_depth, tree_feature_rate, gamma, num_epoch = method_param

	skf = StratifiedKFold(n_splits=num_fold, shuffle=True, random_state=rs)

	pred_data_withPI = []
	pred_data_withoutPI = []

	ens_num = 5

	test_index_list = []

	for i, (train_index, test_index) in enumerate(skf.split(X, label)):

		test_index_list.append(test_index)

		X_train, X_test = X[train_index], X[test_index]
		y_train, y_test = y[train_index], y[test_index]
		PI_train, PI_test = PI_ind[train_index], PI_ind[test_index]

		cutoff = 0.5

		#train ens_num NNs and average the predictions
		for j in range(ens_num):

			[y_pred_withPI, y_withPI], [y_pred_withoutPI, y_withoutPI] = trainer([X_train, y_train, PI_train], [X_test, y_test, PI_test], h_list_main, h_list_joint1, h_list_joint2, dropout, n_tree, tree_depth, tree_feature_rate, gamma, num_epoch)

			if j == 0:
				y_pred_withPI_sum = y_pred_withPI
				y_pred_withoutPI_sum = y_pred_withoutPI
			else:
				y_pred_withPI_sum += y_pred_withPI
				y_pred_withoutPI_sum += y_pred_withoutPI

		y_pred_withPI = y_pred_withPI_sum / float(ens_num)
		y_pred_withoutPI = y_pred_withoutPI_sum / float(ens_num)


		cutoff_withPI = np.ones(len(y_pred_withPI))*cutoff
		cutoff_withoutPI = np.ones(len(y_pred_withoutPI))*cutoff

		pred_data_withPI.append([y_pred_withPI, y_withPI, cutoff_withPI])
		pred_data_withoutPI.append([y_pred_withoutPI, y_withoutPI, cutoff_withoutPI])


	for i in range(len(pred_data_withPI)):
		if i == 0:
			y_pred_withPI_combined = pred_data_withPI[i][0].reshape(-1, 1)
			y_withPI_combined = pred_data_withPI[i][1].reshape(-1, 1)
			cutoff_withPI_combined = pred_data_withPI[i][2].reshape(-1, 1)

			y_pred_withoutPI_combined = pred_data_withoutPI[i][0].reshape(-1, 1)
			y_withoutPI_combined = pred_data_withoutPI[i][1].reshape(-1, 1)
			cutoff_withoutPI_combined = pred_data_withoutPI[i][2].reshape(-1, 1)
		else:
			y_pred_withPI_combined = np.vstack([y_pred_withPI_combined, pred_data_withPI[i][0].reshape(-1, 1)])
			y_withPI_combined = np.vstack([y_withPI_combined, pred_data_withPI[i][1].reshape(-1, 1)])
			cutoff_withPI_combined = np.vstack([cutoff_withPI_combined, pred_data_withPI[i][2].reshape(-1, 1)])

			y_pred_withoutPI_combined = np.vstack([y_pred_withoutPI_combined, pred_data_withoutPI[i][0].reshape(-1, 1)])
			y_withoutPI_combined = np.vstack([y_withoutPI_combined, pred_data_withoutPI[i][1].reshape(-1, 1)])
			cutoff_withoutPI_combined = np.vstack([cutoff_withoutPI_combined, pred_data_withoutPI[i][2].reshape(-1, 1)])


	auc_withPI = roc_auc_score(y_withPI_combined.reshape(-1), y_pred_withPI_combined.reshape(-1))
	aupr_withPI = average_precision_score(y_withPI_combined.reshape(-1), y_pred_withPI_combined.reshape(-1))
	y_pred_binary_withPI = (y_pred_withPI_combined >= cutoff_withPI_combined).astype(float)
	acc_withPI = accuracy_score(y_withPI_combined.reshape(-1), y_pred_binary_withPI.reshape(-1))


	auc_withoutPI = roc_auc_score(y_withoutPI_combined.reshape(-1), y_pred_withoutPI_combined.reshape(-1))
	aupr_withoutPI = average_precision_score(y_withoutPI_combined.reshape(-1), y_pred_withoutPI_combined.reshape(-1))
	y_pred_binary_withoutPI = (y_pred_withoutPI_combined >= cutoff_withoutPI_combined).astype(float)
	acc_withoutPI = accuracy_score(y_withoutPI_combined.reshape(-1), y_pred_binary_withoutPI.reshape(-1))

	y_pred_all = y_pred_withPI_combined.reshape(-1).tolist() + y_pred_withoutPI_combined.reshape(-1).tolist()
	y_pred_binary_all = y_pred_binary_withPI.reshape(-1).tolist() + y_pred_binary_withoutPI.reshape(-1).tolist()
	y_all = y_withPI_combined.reshape(-1).tolist() + y_withoutPI_combined.reshape(-1).tolist()

	auc_all = roc_auc_score(y_all, y_pred_all)
	aupr_all = average_precision_score(y_all, y_pred_all)
	acc_all = accuracy_score(y_all, y_pred_binary_all)

	print ('CV REP:', aupr_all, auc_all, acc_all)

	return [aupr_withPI, auc_withPI, acc_withPI], \
		   [aupr_withoutPI, auc_withoutPI, acc_withoutPI], \
		   [aupr_all, auc_all, acc_all], \
		   [y_pred_withPI_combined.reshape(-1), y_pred_binary_withPI.reshape(-1), y_withPI_combined.reshape(-1)], \
		   [y_pred_withoutPI_combined.reshape(-1), y_pred_binary_withoutPI.reshape(-1), y_withoutPI_combined.reshape(-1)]


def print_result(best_results):
	print ('All aupr mean and std', np.round(best_results[1][1][0],3), '±', np.round(best_results[1][1][1],3))
	print ('All auc mean and std', np.round(best_results[1][0][0],3), '±', np.round(best_results[1][0][1],3))
	print ('All bacc mean and std', np.round(best_results[1][2][0],3), '±', np.round(best_results[1][2][1],3))
	print ('-------------------------------------------------------------')
	print ('PI aupr mean and std', np.round(best_results[2][1][0],3), '±', np.round(best_results[2][1][1],3))
	print ('PI auc mean and std', np.round(best_results[2][0][0],3), '±', np.round(best_results[2][0][1],3))
	print ('PI bacc mean and std', np.round(best_results[2][2][0],3), '±', np.round(best_results[2][2][1],3))
	print ('-------------------------------------------------------------')
	print ('No PI aupr mean and std', np.round(best_results[3][1][0],3), '±', np.round(best_results[3][1][1],3))
	print ('No PI auc mean and std', np.round(best_results[3][0][0],3), '±', np.round(best_results[3][0][1],3))
	print ('No PI bacc mean and std', np.round(best_results[3][2][0],3), '±', np.round(best_results[3][2][1],3))
	print ('================================================')


def Repeat_CV(data, num_rep=5):

	h_list_main = [256, 1]
	h_list_joint1 = [256]
	h_list_joint2 = [256, 256]
	num_epoch = 1000
	dropout = 0.6
	n_tree = 512
	tree_depth = 7 
	tree_feature_rate = 0.7
	gamma = 10

	method_pack = [h_list_main, h_list_joint1, h_list_joint2, dropout, n_tree, tree_depth, tree_feature_rate, gamma, num_epoch]

	rep_auc_list = []
	rep_aupr_list = []
	rep_acc_list = []

	rep_auc_list_pi = []
	rep_aupr_list_pi = []
	rep_acc_list_pi = []

	rep_auc_list_nopi = []
	rep_aupr_list_nopi = []
	rep_acc_list_nopi = []

	rep_raw_pred_list = []

	for a_rep in range(num_rep):
		result_PI, result_noPI, result_all, raw_pred_withPI, raw_pred_withoutPI = CV_NN(data, method_param=method_pack, num_fold=5, rs=a_rep)
							
		rep_auc_list.append(result_all[1])
		rep_aupr_list.append(result_all[0])
		rep_acc_list.append(result_all[2])

		rep_raw_pred_list.append([raw_pred_withPI, raw_pred_withoutPI])

		rep_auc_list_pi.append(result_PI[1])
		rep_aupr_list_pi.append(result_PI[0])
		rep_acc_list_pi.append(result_PI[2])

		rep_auc_list_nopi.append(result_noPI[1])
		rep_aupr_list_nopi.append(result_noPI[0])
		rep_acc_list_nopi.append(result_noPI[2])



		auc_mean, auc_std = np.mean(rep_auc_list), np.std(rep_auc_list)
		aupr_mean, aupr_std = np.mean(rep_aupr_list), np.std(rep_aupr_list)
		acc_mean, acc_std = np.mean(rep_acc_list), np.std(rep_acc_list)

		auc_pi_mean, auc_pi_std = np.mean(rep_auc_list_pi), np.std(rep_auc_list_pi)
		aupr_pi_mean, aupr_pi_std = np.mean(rep_aupr_list_pi), np.std(rep_aupr_list_pi)
		acc_pi_mean, acc_pi_std = np.mean(rep_acc_list_pi), np.std(rep_acc_list_pi)

		auc_nopi_mean, auc_nopi_std = np.mean(rep_auc_list_nopi), np.std(rep_auc_list_nopi)	
		aupr_nopi_mean, aupr_nopi_std = np.mean(rep_aupr_list_nopi), np.std(rep_aupr_list_nopi)
		acc_nopi_mean, acc_nopi_std = np.mean(rep_acc_list_nopi), np.std(rep_acc_list_nopi)
							

		results_raw = [rep_raw_pred_list, \
			[[auc_mean, auc_std], [aupr_mean, aupr_std], [acc_mean, acc_std]], \
			[[auc_pi_mean, auc_pi_std], [aupr_pi_mean, aupr_pi_std], [acc_pi_mean, acc_pi_std]], \
			[[auc_nopi_mean, auc_nopi_std], [aupr_nopi_mean, aupr_nopi_std], [acc_nopi_mean, acc_nopi_std]]]

	print_result(results_raw)

	f = open('NN_results.pkl', 'wb')
	pickle.dump([results_raw, method_pack], f)
	f.close()


Repeat_CV([X, label, PI_vec, feat_dim_list], num_rep=5)


