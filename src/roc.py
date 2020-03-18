import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle

from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
from sklearn.metrics import roc_auc_score
from imblearn.over_sampling import SMOTE

from sklearn.model_selection import GridSearchCV

import sys
import csv

pattern = sys.argv[1]
data_r = pd.read_csv(pattern)
columns = list(data_r)

for i in columns[:-1]:
	print(i)
	data = data_r[[i, 'label']]

	MU = data['label'] == "group1"
	WT = data['label'] == "group2"

	#ADDITIONAL CODE
	X_MU = data[MU].drop('label', axis = 1)
	X_WT = data[WT].drop('label', axis = 1)
	group1SampleSize = len(X_MU)
	group2SampleSize = len(X_WT)

	y_MU = data[MU] .label
	y_WT = data[WT] .label

	y_MU = label_binarize(y_MU, classes=["group1", "group2"])
	y_WT = label_binarize(y_WT, classes=["group1", "group2"])

	#SMOTE PARAMETERIZATION
	X_MU_train, X_MU_test, y_MU_train, y_MU_test = train_test_split(X_MU, y_MU,test_size=0.3)
	X_WT_train, X_WT_test, y_WT_train, y_WT_test = train_test_split(X_WT, y_WT,test_size=0.3)

	#TRAINING
	y_train = np.concatenate([y_MU_train, y_WT_train])
	X_train = np.concatenate([X_MU_train, X_WT_train])
	y_test = np.concatenate([y_MU_test, y_WT_test])
	X_test = np.concatenate([X_MU_test, X_WT_test])

	if group1SampleSize == group2SampleSize:
		X_train_res = X_train
		y_train_res = y_train
	else:
		sm = SMOTE(ratio = 1.0, k_neighbors = 2)
		X_train_res, y_train_res = sm.fit_sample(X_train, y_train.ravel())

	#DEFINING PARAMETER RANGE FOR CLF OPTIMIZATION
	param_grid = {'C': [0.01, 0.1, 0.05, 1, 1.05, 1.1, 10, 100, 1000], 'gamma': [10, 1, 0.1, 0.01, 0.001, 0.0001], 'kernel': ['rbf']}
	grid = GridSearchCV(svm.SVC(), param_grid, refit = True, verbose = False, iid = True, cv = 5)
	grid = grid.fit(X_train_res, y_train_res.ravel())

	#CLF RESULTS
	clf = svm.SVC(kernel = 'rbf', C = grid.best_params_['C'], gamma = grid.best_params_['gamma']).fit(X_train_res, y_train_res.ravel())
	#y_score = clf.predict(X_test)
	y_score = clf.fit(X_train_res, y_train_res).decision_function(X_test)

	#SEGMENTING X AS DATA AND Y AS TARGET
	X = data.drop('label', axis = 1)
	y = data['label']

	# Binarize the output
	y = label_binarize(y, classes=["group1", "group2"])
	n_classes = y.shape[1]

	# Compute ROC curve and ROC area for each class
	fpr = dict()
	tpr = dict()
	roc_auc = dict()

	print("Examining y score")
	print(y_score)

	#for i in range(n_classes):
	fpr[0], tpr[0], _ = roc_curve(y_test, y_score)
	roc_auc[0] = auc(fpr[0], tpr[0])

	# Compute micro-average ROC curve and ROC area
	fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
	roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

	plt.figure()
	lw = 2

	print(roc_auc)

	plt.plot(fpr[0], tpr[0], color='darkorange', lw=lw, label='ROC curve (area = %0.2f)' % roc_auc[0])
	plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.02])
	plt.xlabel('False Positive Rate - 1 minus Specificity')
	plt.ylabel('True Positive Rate - Sensitivity')
	plt.title('Receiver Operating Characteristic for ' + i)
	plt.legend(loc="lower right")
	plt.savefig(i + 'roc_curve.png')
	plt.close()
