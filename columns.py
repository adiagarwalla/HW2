import math
import linecache
import numpy as np
from sklearn import linear_model

file_sample = "data/intersected_final_chr1_cutoff_20_sample.bed"
file_test = "data/intersected_final_chr1_cutoff_20_test.bed"
file_train = "data/intersected_final_chr1_cutoff_20_train_revised.bed"
file_bed = "data/intersected_final_chr1_cutoff_20.bed"

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

def main():
    read_sample_nonNan()

main()
