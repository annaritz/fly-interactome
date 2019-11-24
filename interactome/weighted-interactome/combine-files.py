ORIG_FILE = '../interactome-flybase-collapsed-evidence.txt'
WEIGHTED_FILE = 'final-weighted-interactome.txt'
OUTFILE = 'interactome-flybase-collapsed-weighted.txt'

## simple script to insert a 3rd column into the original file that includes the weight.
## TODO make these not hard-coded.
## TODO incorporate this into the previous scripts that weights the interactome.

def main():
	## combine infile and weighted file:
	out = open(OUTFILE,'w')
	out.write('#symbol1\tsymbol2\tweight\tPubMedIDs\tFlyBase1\tFlyBase2\tDBs\tEvidence\n')

	## get weights
	weights = {}
	with open(WEIGHTED_FILE) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			weights[(row[0],row[1])] = float(row[2])
			weights[(row[1],row[0])] = float(row[2])
	
	## combine files
	self_edge = 0
	other_skipped = 0
	with open(ORIG_FILE) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			edge = (row[0],row[1])
			if edge not in weights:
				if row[0] == row[1]:
					self_edge += 1
				else:
					other_skipped +=1
					print('WARNING: (%s,%s) not in weights' % (row[0],row[1]))
				continue
			new_row = row[0:2] + [str(weights[(row[0],row[1])])] + row[2:]
			out.write('\t'.join(new_row)+'\n')
	out.close()
	print('wrote to %s' % (OUTFILE))
	print('skipped %d self edges and %d other edges' % (self_edge,other_skipped))

	return

if __name__ == '__main__':
	main()