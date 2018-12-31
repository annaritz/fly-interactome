from __future__ import print_function
import matplotlib.pyplot as plt
import sys


def main(args):
	infile = args[1]
	col = int(args[2])
	num_bins = int(args[3])
	title = args[4]
	xlabel = args[5]
	ylabel = args[6]
	outfile = args[7]

	lines = open(infile).readlines()
	data = [float(line.strip().split()[col]) for line in lines if line[0] != '#']

	fig, ax = plt.subplots()
	ax.hist(data, num_bins)
	ax.set_title(title)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)

	plt.tight_layout()
	plt.savefig(outfile)
	print('Wrote to ',outfile)

if __name__ == '__main__':
	if len(sys.argv) != 8:
		sys.exit('USAGE: python plot-hist.py <INFILE> <COL> <NUM_BINS> <TITLE> <XLABEL> <YLABEL> <OUTFILE>')
	main(sys.argv)