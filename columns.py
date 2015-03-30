file_sample = "../intersected_final_chr1_cutoff_20_sample.bed"
file_test = "../intersected_final_chr1_cutoff_20_test.bed"
file_train = "../intersected_final_chr1_cutoff_20_train_revised.bed"
file_bed = "../intersected_final_chr1_cutoff_20.bed"

# Get non-NAN sites from sample file
def read_sample_nonNan():
	print "Getting started"
	test = []
	with open(file_sample) as f:
		test = f.readLines()

	print test

def main():
	read_sample_nonNan()

main()


    
