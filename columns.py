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
    X = []
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
                X.append(train_data2)
    print X[0:10]
                
    # XMatrix = np.array(X)
    # YMatrix = np.array(Y)
    
    clf = linear_model.LinearRegression()
    clf.fit(np.array(X, 'float_'), np.array(Y, 'float_'))
    print clf.coef_

    Y_test = []

    with open(file_test) as f_test:
    	for k, line_test in enumerate(f_test, 1):
    		test_sample = line_test.strip().split()
    		if (test_sample[5] == '0'):
    			Y_test.append(test_sample[5])

    # The mean square error
    print("Residual sum of squares: %.2f"
      % np.mean((clf.predict(np.array(X, 'float_')) - np.array(Y_test, 'float_')) ** 2))

    # The root mean square error
    print("Residual root sum of squares: %.2f"
      % math.sqrt(np.mean((clf.predict(np.array(X, 'float_')) - np.array(Y_test, 'float_')) ** 2)))

def main():
    read_sample_nonNan()

main()
