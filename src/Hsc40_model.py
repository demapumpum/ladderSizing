# Training machine learning model for Hsc40 data
# Classifiers used:
# 					Random Forest
# 					KNeighbors

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
import pickle
from imblearn.over_sampling import SMOTE
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.neighbors import KNeighborsClassifier

def randomForestClassifier(X_train, y_train, n_fold):
	parameter_RandomForest = {
		'n_estimators': (10, 30, 50, 100, 200, 300, 400, 500, 700, 800, 1000),
		'max_features': ('auto', 'sqrt', 'log2', None),
		'criterion': ('gini', 'entropy')
	}
	gs_RandomForest = GridSearchCV(RandomForestClassifier(), parameter_RandomForest, verbose=1, cv=n_fold, n_jobs=-1)
	model = gs_RandomForest.fit(X_train, y_train)
	return model

def kNeighborsClassifier(X_train, y_train, n_fold):
	parameter_KNeighbors = {
		'n_neighbors': (3, 5, 7, 9, 11, 13, 15),
		'weights': ('uniform', 'distance'),
		'algorithm': ('ball_tree', 'kd_tree', 'brute', 'auto')
	}
	gs_KNeighbors = GridSearchCV(KNeighborsClassifier(), parameter_KNeighbors, verbose=1, cv=n_fold, n_jobs=-1)
	model = gs_KNeighbors.fit(X_train, y_train)
	return model

file = pd.read_csv('training.csv')			# Access data for feeding the model

filename = []
label = []
area = []
length = []
height = []

for i in range(len(file)):										# Read csv file and store features in array
	filename.append(file.iat[i, 0])
	label.append(file.iat[i, 1])
	area.append(file.iat[i, 2])
	length.append(file.iat[i, 3])
	height.append(file.iat[i, 4])

matrix = [area, length, height]									# put features in a matrix
X_train = np.column_stack(matrix)
y_train = label													# label (1 or 0)

# X_train, X_test, y_train, y_test = train_test_split(with_filename, y, test_size=0.25)	# Split data 75% training 25% testing

# test_matrix = [X_test, y_test]								# Put test data
# test_data = np.column_stack(test_matrix)						# in a matrix and
# pd.DataFrame(test_data).to_csv("test_data.csv", index=False)	# store in a csv file (to know which data are unseen by the model)

# X_test = X_test[:, 1:]										# Drop 'filename' in matrix
# X_train = X_train[:, 1:]										# used only to know from which files testing data comes from

# pd.DataFrame(X_train).to_csv("X_train.csv", index=False)		# Store training and testing data in a csv file
# pd.DataFrame(X_test).to_csv("X_test.csv", index=False)		# These files will be used in creating ROC and confusion matrix
# pd.DataFrame(y_train).to_csv("y_train.csv", index=False)
# pd.DataFrame(y_test).to_csv("y_test.csv", index=False)

model1 = randomForestClassifier(X_train, y_train, 25)

model2 = kNeighborsClassifier(X_train, y_train, 25)

																# Print confusion matrix, classification report, and accuracy score for both models
print('\nFor Random Forest:\n')									# Also prints the model with highest score, best parameter for the models from GridSearchCV
print(model1.best_score_)
print(model1.best_estimator_)
print(model1.best_params_)
print('\n\n')

print('\nFor KNeighbors:\n')
print(model2.best_score_)
print(model2.best_estimator_)
print(model2.best_params_)
print('\n\n')

filename1 = 'Hsc40_model_RandomForest_improved.sav'
filename2 = 'Hsc40_model_KNeighbors_improved.sav'
pickle.dump(model1, open(filename1, 'wb'))						# Store Random Forest model in 'Hsc40_model_RandomForest' through pickle
pickle.dump(model2, open(filename2, 'wb'))						# Store KNeighbors model in 'Hsc40_model_KNeighbors'