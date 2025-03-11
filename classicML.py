import numpy as np
import scipy
import csv
import copy
import pickle
import math
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from ngboost import NGBClassifier as ngb
from ngboost.scores import MLE, LogScore, CRPScore
from ngboost.distns import k_categorical, Bernoulli ##Classifier
from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier
from xgboost import XGBClassifier
from sklearn.impute import KNNImputer
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.metrics import balanced_accuracy_score as accuracy_score

#Load dataset
f = open('./processed_data/dataset.pkl', 'rb')
dataset = pickle.load(f)
f.close()

data_feat, label, fidelity_vec, PI_vec, ESM_first_basis, ESM_mutated_basis, feat_dim_list = dataset
ESM_first_dim, ESM_mutated_dim, group_code_dim, basic_feat_dim, MD_feat_dim, dock_feat_dim = feat_dim_list

#Set K of KNN imputer 
knn_num = math.ceil(np.sqrt(len(data_feat)))

#PCA for ESM2 feats
PCA_dim = 50
PCA_first = PCA(n_components=PCA_dim, random_state=18)
PCA_mutated = PCA(n_components=PCA_dim, random_state=18)
PCA_first.fit(ESM_first_basis)
PCA_mutated.fit(ESM_mutated_basis)

#Replace original ESM2 feats with dimensionality-reduced ones
ESM_first_PCA = PCA_first.transform(data_feat[:, :ESM_first_dim])
ESM_mutated_PCA = PCA_mutated.transform(data_feat[:, ESM_first_dim:ESM_first_dim + ESM_mutated_dim])
X = copy.copy(np.hstack([ESM_first_PCA, ESM_mutated_PCA, data_feat[:, ESM_first_dim + ESM_mutated_dim:]]))


def CV_classicML(data, method_param, use_PI, use_imputation, normalization, num_fold, rs):

	X, y, PI_ind, feat_dim_list = data #PI_ind binary vector indicating whether the sample has PI features
	ESM_first_dim, ESM_mutated_dim, group_code_dim, basic_feat_dim, MD_feat_dim, dock_feat_dim = feat_dim_list


	skf = StratifiedKFold(n_splits=num_fold, shuffle=True, random_state=rs)

	if use_PI == False:
		X = X[:, :PCA_dim + PCA_dim + group_code_dim + basic_feat_dim]


	pred_data = []

	for i, (train_index, test_index) in enumerate(skf.split(X, label)):

		X_train, X_test = X[train_index], X[test_index]
		y_train, y_test = y[train_index], y[test_index]
		PI_train, PI_test = PI_ind[train_index], PI_ind[test_index]

		cutoff = 0.5 #binarization cutoff
		cutoff_array = np.ones(len(y_test))*cutoff

		if use_imputation == True:
			imp = KNNImputer(missing_values=np.nan, n_neighbors=knn_num, weights='distance')
			imp.fit(X_train)
			X_train = imp.transform(X_train)
			X_test = imp.transform(X_test)


		if normalization == True:
			normalizer = StandardScaler()
			normalizer.fit(X_train)
			X_train = normalizer.transform(X_train)
			X_test = normalizer.transform(X_test)

		method_name, method_params_dict = method_param


		if method_name == 'NGB':
			max_depth = method_params_dict['max_depth']
			n_estimators = method_params_dict['n_estimators']
			learner = DecisionTreeRegressor(criterion='friedman_mse', max_depth=max_depth, random_state = 18)
			clf = ngb(Base=learner, n_estimators = n_estimators, Score = LogScore, Dist = Bernoulli, random_state = 18)
			try:
				clf.fit(X_train, y_train.astype(int))
				y_pred = clf.predict_proba(X_test)[:, 1]
			except:
				y_pred = np.ones(len(X_test))*0.5

		if method_name == 'GBDT':
			max_depth = method_params_dict['max_depth']
			n_estimators = method_params_dict['n_estimators']
			clf = XGBClassifier(n_estimators=n_estimators, max_depth=max_depth, n_jobs=20, random_state=18)
			clf.fit(X_train, y_train)
			y_pred = clf.predict_proba(X_test)[:, 1]

		if method_name == 'SVM':
			C = method_params_dict['C']
			clf = SVC(C=C, random_state=18, probability=True)
			clf.fit(X_train, y_train)
			y_pred = clf.predict_proba(X_test)[:, 1]

		pred_data.append([y_pred, y_test, PI_test, cutoff_array])

	#Combine predictions from all test folds to calculate micro avg performance metrics
	for i in range(len(pred_data)):
		if i == 0:
			y_pred_combined = pred_data[i][0].reshape(-1, 1)
			y_test_combined = pred_data[i][1].reshape(-1, 1)
			PI_test_combined = pred_data[i][2].reshape(-1, 1)
			cutoff_combined = pred_data[i][3].reshape(-1, 1)
		else:
			y_pred_combined = np.vstack([y_pred_combined, pred_data[i][0].reshape(-1, 1)])
			y_test_combined = np.vstack([y_test_combined, pred_data[i][1].reshape(-1, 1)])
			PI_test_combined = np.vstack([PI_test_combined, pred_data[i][2].reshape(-1, 1)])
			cutoff_combined = np.vstack([cutoff_combined, pred_data[i][3].reshape(-1, 1)])


	auc = roc_auc_score(y_test_combined.reshape(-1), y_pred_combined.reshape(-1))
	aupr = average_precision_score(y_test_combined.reshape(-1), y_pred_combined.reshape(-1))
	y_pred_binary = (y_pred_combined >= cutoff_combined).astype(float)
	acc = accuracy_score(y_test_combined.reshape(-1), y_pred_binary.reshape(-1))


	index_pi = np.where(PI_test_combined.reshape(-1) == 1)
	index_withoutpi = np.where(PI_test_combined.reshape(-1) == 0)

	auc_pi = roc_auc_score(y_test_combined[index_pi].reshape(-1), y_pred_combined[index_pi].reshape(-1))
	aupr_pi = average_precision_score(y_test_combined[index_pi].reshape(-1), y_pred_combined[index_pi].reshape(-1))
	y_pred_binary = (y_pred_combined[index_pi] >= cutoff_combined[index_pi]).astype(float)
	acc_pi = accuracy_score(y_test_combined[index_pi].reshape(-1), y_pred_binary.reshape(-1))

	auc_withoutpi = roc_auc_score(y_test_combined[index_withoutpi].reshape(-1), y_pred_combined[index_withoutpi].reshape(-1))
	aupr_withoutpi = average_precision_score(y_test_combined[index_withoutpi].reshape(-1), y_pred_combined[index_withoutpi].reshape(-1))
	y_pred_binary = (y_pred_combined[index_withoutpi] >= cutoff_combined[index_withoutpi]).astype(float)
	acc_withoutpi = accuracy_score(y_test_combined[index_withoutpi].reshape(-1), y_pred_binary.reshape(-1))

	return [[aupr, auc, acc], [aupr_pi, auc_pi, acc_pi], [aupr_withoutpi, auc_withoutpi, acc_withoutpi], [y_pred_combined.reshape(-1), y_test_combined.reshape(-1), y_pred_binary.reshape(-1), PI_test_combined.reshape(-1)]]


def Repeat_CV(data, num_rep, method_param, use_PI, use_imputation, normalization, num_fold):

	rep_auc_list = []
	rep_aupr_list = []
	rep_acc_list = []

	rep_auc_list_pi = []
	rep_aupr_list_pi = []
	rep_acc_list_pi = []

	rep_auc_list_nopi = []
	rep_aupr_list_nopi = []
	rep_acc_list_nopi = []

	rep_pred_list = []

	for a_rep in range(num_rep):
		result_pack = CV_classicML(data=data, method_param=method_param, use_PI=use_PI, use_imputation=use_imputation, normalization=normalization, num_fold=num_fold, rs=a_rep)
		[aupr, auc, acc] = result_pack[0]
		[aupr_pi, auc_pi, acc_pi] = result_pack[1]
		[aupr_withoutpi, auc_withoutpi, acc_withoutpi] = result_pack[2]
		[y_pred_all, y_test_all, y_pred_binary_all, PI_all] = result_pack[3]

		rep_pred_list.append([y_pred_all, y_test_all, y_pred_binary_all, PI_all])

		rep_auc_list.append(auc)
		rep_aupr_list.append(aupr)
		rep_acc_list.append(acc)

		rep_auc_list_pi.append(auc_pi)
		rep_aupr_list_pi.append(aupr_pi)
		rep_acc_list_pi.append(acc_pi)

		rep_auc_list_nopi.append(auc_withoutpi)
		rep_aupr_list_nopi.append(aupr_withoutpi)
		rep_acc_list_nopi.append(acc_withoutpi)

	auc_mean, auc_std = np.mean(rep_auc_list), np.std(rep_auc_list)
	aupr_mean, aupr_std = np.mean(rep_aupr_list), np.std(rep_aupr_list)
	acc_mean, acc_std = np.mean(rep_acc_list), np.std(rep_acc_list)

	auc_pi_mean, auc_pi_std = np.mean(rep_auc_list_pi), np.std(rep_auc_list_pi)
	aupr_pi_mean, aupr_pi_std = np.mean(rep_aupr_list_pi), np.std(rep_aupr_list_pi)
	acc_pi_mean, acc_pi_std = np.mean(rep_acc_list_pi), np.std(rep_acc_list_pi)

	auc_nopi_mean, auc_nopi_std = np.mean(rep_auc_list_nopi), np.std(rep_auc_list_nopi)
	aupr_nopi_mean, aupr_nopi_std = np.mean(rep_aupr_list_nopi), np.std(rep_aupr_list_nopi)
	acc_nopi_mean, acc_nopi_std = np.mean(rep_acc_list_nopi), np.std(rep_acc_list_nopi)

	return [rep_pred_list, \
	[[auc_mean, auc_std], [aupr_mean, aupr_std], [acc_mean, acc_std]], \
	[[auc_pi_mean, auc_pi_std], [aupr_pi_mean, aupr_pi_std], [acc_pi_mean, acc_pi_std]], \
	[[auc_nopi_mean, auc_nopi_std], [aupr_nopi_mean, aupr_nopi_std], [acc_nopi_mean, acc_nopi_std]]]


def hyper_tune(data, num_rep):


	max_bacc = 0.0
	for use_PI in [True, False]:
		for use_imputation in [True, False]:
			for i in [32, 128, 512]:
				for j in [4, 8, 16, 32, None]:
					param_dict = {}
					param_dict['n_estimators'] = i
					param_dict['max_depth'] = j
					method_pack = ['NGB', param_dict]

					results = Repeat_CV(data, num_rep, method_param=method_pack, use_PI=use_PI, use_imputation=use_imputation, normalization=True, num_fold=5)

					bacc_mean = results[1][0][0]

					if bacc_mean > max_bacc:
						max_bacc = bacc_mean
						best_results = results
						best_params = [use_PI, use_imputation, method_pack]
						f = open('NGB_best_results.pkl', 'wb')
						pickle.dump([best_results, best_params], f)
						f.close()

	max_bacc = 0.0
	for use_PI in [True, False]:
		for use_imputation in [True, False]:
			for i in [32, 128, 512]:
				for j in [4, 8, 16, 32, None]:
					param_dict = {}
					param_dict['n_estimators'] = i
					param_dict['max_depth'] = j
					method_pack = ['GBDT', param_dict]

					results = Repeat_CV(data, num_rep, method_param=method_pack, use_PI=use_PI, use_imputation=use_imputation, normalization=True, num_fold=5)

					bacc_mean = results[1][0][0]

					if bacc_mean > max_bacc:
						max_bacc = bacc_mean
						best_results = results
						best_params = [use_PI, use_imputation, method_pack]
						f = open('GBDT_best_results.pkl', 'wb')
						pickle.dump([best_results, best_params], f)
						f.close()

						print (best_results)

	max_bacc = 0.0
	for C in [10, 1, 0.1, 0.01]:
		for use_PI in [True, False]:
			for use_imputation in [True, False]:
				param_dict = {}
				param_dict['C'] = C
				method_pack = ['SVM', param_dict]

				results = Repeat_CV(data, num_rep, method_param=method_pack, use_PI=use_PI, imputation_param=use_imputation, normalization=True, num_fold=5)

				bacc_mean = results[1][0][0]

				if bacc_mean > max_bacc:
					max_bacc = bacc_mean
					best_results = results
					best_params = [use_PI, use_imputation, method_pack]
					f = open('SVM_best_results.pkl', 'wb')
					pickle.dump([best_results, best_params], f)
					f.close()
					
					print (best_results)


hyper_tune([X, label, PI_vec, feat_dim_list], num_rep=5)

