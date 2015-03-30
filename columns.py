import numpy as np
from sklearn import linear_model

file_sample = "../intersected_final_chr1_cutoff_20_sample.bed"
file_test = "../intersected_final_chr1_cutoff_20_test.bed"
file_train = "../intersected_final_chr1_cutoff_20_train_revised.bed"
file_bed = "../intersected_final_chr1_cutoff_20.bed"

# Get non-NAN sites from sample file
def read_sample_nonNan():
	print "Getting started"
	f_train = open(file_train)
	test = []
	with open(file_sample) as f:
		for i, line in enumerate(f, 1):
			output_line = line.strip().split()
			if (output_line[5] == '1'):
				clf = linear_model.LinearRegression()


	print test

def main():
	read_sample_nonNan()

main()