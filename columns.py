import math
import linecache
import numpy as np
from operator import itemgetter
from sklearn import linear_model
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor

file_sample_1 = "data/intersected_final_chr1_cutoff_20_sample.bed"
file_test_1 = "data/intersected_final_chr1_cutoff_20_test.bed"
file_train_1 = "data/intersected_final_chr1_cutoff_20_train_revised.bed"
file_bed_1 = "data/intersected_final_chr1_cutoff_20.bed"

file_sample_10 = "data/intersected_final_chr10_cutoff_20_sample.bed"
file_test_10 = "data/intersected_final_chr10_cutoff_20_test.bed"
file_train_10 = "data/intersected_final_chr10_cutoff_20_train_revised.bed"
file_bed_10 = "data/intersected_final_chr10_cutoff_20.bed"

file_sample_21 = "data/intersected_final_chr21_cutoff_20_sample.bed"
file_test_21 = "data/intersected_final_chr21_cutoff_20_test.bed"
file_train_21 = "data/intersected_final_chr21_cutoff_20_train_revised.bed"
file_bed_21 = "data/intersected_final_chr21_cutoff_20.bed"

# Assign file for current running
file_sample = file_sample_21;
file_test = file_test_21;
file_train = file_train_21;
file_bed = file_bed_21;

# def knn(n):
#     print "Starting KNN calculation with n=%d" % n

#     f_train = open(file_train)
#     X_train = []
#     Y = []
#     with open(file_sample) as f:
#         for i, line in enumerate(f, 1):
#             sample = line.strip().split()
#             if sample[4] != 'nan':
#                 Y.append(sample[4])
#                 row = linecache.getline(file_train, i).strip().split()
#                 train_data = row[4: -1]
#                 neighbors = []
#                 total = 0
#                 count = 0
#                 for j in train_data:
#                     if j != 'nan':
#                         total = total + float(j)
#                         count = count + 1
#                 avg = total / float(count)
#                 train_data2 = []
#                 for j in train_data:
#                     if j == 'nan':
#                         train_data2.append(avg)
#                     else:
#                         train_data2.append(float(j))
#                 for d in train_data2:
#                     diff = ((float(sample[4]) - float(d)) ** 2)
#                     if len(neighbors < n):
#                         neighbors.append((d, diff))
#                     else:
#                         sorted(neighbors, key=itemgetter(1), reverse=True)
#                         if diff < neighbors[0][1]:
#                             neighbors = neighbors[1:]
#                             neighbors.append((d, diff))
#                 n_closest = []
#                 for pair in neighbors:
#                     n_closest.append(pair[0])
#                 X_train.append(n_closest)
#             elif

# Get non-NAN sites from sample file
def read_sample_nonNan():
    print "Getting started"
    f_train = open(file_train)
    X_train = []
    X_test = []
    Y = []
    with open(file_sample) as f:
        for i, line in enumerate(f, 1):
            sample = line.strip().split()
            if (sample[4] != 'nan'):
                Y.append(sample[4])
                row = linecache.getline(file_train, i).strip().split()
                train_data = row[4: -1]
                total = 0
                count = 0
                for j in train_data:
                    if j != 'nan':
                        total = total + float(j)
                        count = count + 1
                avg = total / float(count)
                train_data2 = []
                for j in train_data:
                    if j == 'nan':
                        train_data2.append(avg)
                    else:
                        train_data2.append(float(j))
                X_train.append(train_data2)
            elif (sample[5] == '0'):
                row = linecache.getline(file_train, i).strip().split()
                train_data = row[4: -1]
                total = 0
                count = 0
                for j in train_data:
                    if j != 'nan':
                        total = total + float(j)
                        count = count + 1
                avg = total / float(count)
                train_data2 = []
                for j in train_data:
                    if j == 'nan':
                        train_data2.append(avg)
                    else:
                        train_data2.append(float(j))
                X_test.append(train_data2)
                
    clf = linear_model.LinearRegression()
    clf.fit(np.array(X_train, 'float_'), np.array(Y, 'float_'))
    print clf.coef_

    Y_test = []

    X_test_final = []
    Y_test_final = []


    with open(file_test) as f_test:
        for k, line_test in enumerate(f_test, 1):
            test_sample = line_test.strip().split()
            if (test_sample[5] == '0'):
                Y_test.append(test_sample[4])

    for i, j in enumerate(Y_test):
        if j != 'nan':
            Y_test_final.append(j)
            X_test_final.append(X_test[i])

    print len(X_test_final)
    print len(Y_test_final)

    # The mean square error
    rss = np.mean((clf.predict(np.array(X_test_final, 'float_')) - np.array(Y_test_final, 'float_')) ** 2)
    print("Residual sum of squares: %.5f" % rss)

    # The root mean square error
    print("Residual root sum of squares: %.5f" % math.sqrt(rss))

    # Ridge
    clf_ridge = linear_model.Ridge (alpha = .5)
    clf_ridge.fit(np.array(X_train, 'float_'), np.array(Y, 'float_'))
    print clf_ridge.coef_

    # The mean square error
    rss_ridge = np.mean((clf_ridge.predict(np.array(X_test_final, 'float_')) - np.array(Y_test_final, 'float_')) ** 2)
    print("Residual sum of squares (Ridge): %.5f" % rss_ridge)

    # The root mean square error
    print("Residual root sum of squares (Ridge): %.5f" % math.sqrt(rss_ridge))

    # Lasso
    clf_lasso = linear_model.Lasso (alpha = .1)
    clf_lasso.fit(np.array(X_train, 'float_'), np.array(Y, 'float_'))
    
    # The mean square error
    rss_lasso = np.mean((clf_lasso.predict(np.array(X_test_final, 'float_')) - np.array(Y_test_final, 'float_')) ** 2)
    print("Residual sum of squares (Lasso): %.5f" % rss_lasso)

    # The root mean square error
    print("Residual root sum of squares (Lasso): %.5f" % math.sqrt(rss_lasso))

    # Bayesian Ridge Regression
    clf_bayes = linear_model.BayesianRidge()
    clf_bayes.fit(np.array(X_train, 'float_'), np.array(Y, 'float_'))
    
    # The mean square error
    rss_bayes = np.mean((clf_bayes.predict(np.array(X_test_final, 'float_')) - np.array(Y_test_final, 'float_')) ** 2)
    print("Residual sum of squares (Bayes): %.5f" % rss_bayes)

    # The root mean square error
    print("Residual root sum of squares (Bayes): %.5f" % math.sqrt(rss_bayes))   

    # Decision Tree Regressor
    clf_decision = DecisionTreeRegressor(max_depth=2)
    clf_decision.fit(np.array(X_train, 'float_'), np.array(Y, 'float_'))
    
    # The mean square error
    rss_decision = np.mean((clf_decision.predict(np.array(X_test_final, 'float_')) - np.array(Y_test_final, 'float_')) ** 2)
    print("Residual sum of squares (Decision): %.5f" % rss_decision)

    # The root mean square error
    print("Residual root sum of squares (Decision): %.5f" % math.sqrt(rss_decision))

    # Random Forest Regressor
    num_trees = 100
    clf_forest = RandomForestRegressor(n_estimators = num_trees)
    clf_forest.fit(np.array(X_train, 'float_'), np.array(Y, 'float_'))
    
    # The mean square error
    rss_forest = np.mean((clf_forest.predict(np.array(X_test_final, 'float_')) - np.array(Y_test_final, 'float_')) ** 2)
    print("Residual sum of squares (Random Forest): %.5f" % rss_forest)

    # The root mean square error
    print("Residual root sum of squares (Random Forest): %.5f" % math.sqrt(rss_forest))

def main():
    read_sample_nonNan()
    # knn(5)

main()
