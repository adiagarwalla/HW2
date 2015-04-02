import math
import linecache
import numpy as np
from operator import itemgetter
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor

file_sample_1 = "data/intersected_final_chr1_cutoff_20_sample.bed"
file_test_1 = "data/intersected_final_chr1_cutoff_20_test.bed"
file_train_1 = "data/intersected_final_chr1_cutoff_20_train_revised.bed"

file_sample_10 = "data/intersected_final_chr10_cutoff_20_sample.bed"
file_test_10 = "data/intersected_final_chr10_cutoff_20_test.bed"
file_train_10 = "data/intersected_final_chr10_cutoff_20_train_revised.bed"

file_sample_21 = "data/intersected_final_chr21_cutoff_20_sample.bed"
file_test_21 = "data/intersected_final_chr21_cutoff_20_test.bed"
file_train_21 = "data/intersected_final_chr21_cutoff_20_train_revised.bed"

# Assign file for current running
file_sample = file_sample_1;
file_test = file_test_1;
file_train = file_train_1;

def method_error(X_train, Y_train, X_test, Y_test):
    # Linear Regression
    print "Linear Regression:"
    clf = linear_model.LinearRegression()
    clf.fit(X_train, Y_train)
    Y_pred = clf.predict(X_test)

    # R^2 calculation
    r2 = r2_score(Y_test, Y_pred)
    print("R^2 regression score: %.5f" % r2)

    # The root mean square error
    rmse = mean_squared_error(Y_test, Y_pred)**0.5
    print("Root mean sum of squares: %.5f" % rmse)
    print

    # Ridge
    print "Ridge regression:"
    clf_ridge = linear_model.Ridge(alpha = .5)
    clf_ridge.fit(X_train, Y_train)
    Y_pred = clf_ridge.predict(X_test)

    # R^2 calculation
    r2 = r2_score(Y_test, Y_pred)
    print("R^2 regression score: %.5f" % r2)

    # The root mean square error
    rmse = mean_squared_error(Y_test, Y_pred)**0.5
    print("Root mean sum of squares: %.5f" % rmse)
    print

    # Lasso
    print "Lasso:"
    clf_lasso = linear_model.Lasso (alpha = .001)
    clf_lasso.fit(X_train, Y_train)
    Y_pred = clf_lasso.predict(X_test)

    # R^2 calculation
    r2 = r2_score(Y_test, Y_pred)
    print("R^2 regression score: %.5f" % r2)

    # The root mean square error
    rmse = mean_squared_error(Y_test, Y_pred)**0.5
    print("Root mean sum of squares: %.5f" % rmse)
    print

    # Bayesian Ridge Regression
    print "Bayesian Ridge Regression:"
    clf_bayes = linear_model.BayesianRidge()
    clf_bayes.fit(X_train, Y_train)
    Y_pred = clf_bayes.predict(X_test)

    # R^2 calculation
    r2 = r2_score(Y_test, Y_pred)
    print("R^2 regression score: %.5f" % r2)

    # The root mean square error
    rmse = mean_squared_error(Y_test, Y_pred)**0.5
    print("Root mean sum of squares: %.5f" % rmse)
    print

    # Decision Tree Regressor
    print "Decision Tree Regressor:"
    clf_decision = DecisionTreeRegressor(max_depth=2)
    clf_decision.fit(X_train, Y_train)
    Y_pred = clf_decision.predict(X_test)

    # R^2 calculation
    r2 = r2_score(Y_test, Y_pred)
    print("R^2 regression score: %.5f" % r2)

    # The root mean square error
    rmse = mean_squared_error(Y_test, Y_pred)**0.5
    print("Root mean sum of squares: %.5f" % rmse)
    print

    # Random Forest Regressor
    print "Random Forest Regressor:"
    num_trees = 100
    clf_forest = RandomForestRegressor(n_estimators = num_trees)
    clf_forest.fit(X_train, Y_train)
    Y_pred = clf_forest.predict(X_test)

    # R^2 calculation
    r2 = r2_score(Y_test, Y_pred)
    print("R^2 regression score: %.5f" % r2)

    # The root mean square error
    rmse = mean_squared_error(Y_test, Y_pred)**0.5
    print("Root mean sum of squares: %.5f" % rmse)

def knn(k):
    print "Running analysis with KNN where K=%d" % k

    X_train = []
    Y_train = []
    X_test  = []

    f_train = open(file_train)
    with open(file_sample) as f:

        # Construct X_train, X_RSS, and Y_train
        for i, line in enumerate(f, 1):
            sample = line.strip().split()

            # Construct X_train and Y_train
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
            if sample[4] != 'nan':
                Y_train.append(float(sample[4]))
                X_train.append(train_data2)
            # Construct X_test
            elif sample[5] == '0':
                X_test.append(train_data2)
                
        # Construct X_RSS
        X_RSS = []
        for j, y in enumerate(Y_train):
            rss_row = []
            for l in X_train[j]:
                rss = ((y - l) ** 2)
                rss_row.append(rss)
            X_RSS.append(rss_row)

        # Get X_RSS column sums
        X_RSS_sums = []
        for i, rss_row in enumerate(X_RSS):
            if i == 0:
                for j, rss in enumerate(rss_row):
                    X_RSS_sums.append([rss, j])
            else:
                for j, rss in enumerate(rss_row):
                    X_RSS_sums[j][0] += rss
        X_RSS_sums.sort(key=lambda x: x[0])

        # Find the k closest tissues
        k_closest_columns = []
        for i, pair in enumerate(X_RSS_sums):
            if i < k:
                k_closest_columns.append(pair[1])

        # Construct new training and testing sets with only those tissues
        X_train_pared = []
        X_test_pared  = []

        for row in X_train:
            new_row = []
            for i in k_closest_columns: 
                new_row.append(row[i])
            X_train_pared.append(new_row)

        for row in X_test:
            new_row = []
            for i in k_closest_columns:
                new_row.append(row[i])
            X_test_pared.append(new_row)

    Y_test  = []
    X_test_final = []
    Y_test_final = []

    with open(file_test) as f_test:
        for k, line_test in enumerate(f_test, 1):
            test_sample = line_test.strip().split()
            if (test_sample[5] == '0'):
                Y_test.append(test_sample[4])

    for i, j in enumerate(Y_test):
        if j != 'nan':
            Y_test_final.append(float(j))
            X_test_final.append(X_test_pared[i])

    method_error(X_train_pared, Y_train, X_test_final, Y_test_final)

# Get non-NAN sites from sample file
def allSamples():
    print "Running analysis for all samples"
    f_train = open(file_train)
    X_train = []
    X_test  = []
    Y_train = []
    with open(file_sample) as f:
        for i, line in enumerate(f, 1):
            sample = line.strip().split()
            if (sample[4] != 'nan'):
                Y_train.append(float(sample[4]))
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
            Y_test_final.append(float(j))
            X_test_final.append(X_test[i])

    method_error(X_train, Y_train, X_test_final, Y_test_final)
    print
    
def main():
    allSamples()
    for i in range(1, 16):
        knn(i)

main()
