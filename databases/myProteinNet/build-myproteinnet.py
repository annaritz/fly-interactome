import sys
'''
Parses the FlyBase version of the myProteinNet.  
Replaces PubMed IDs with "NA"
'''

def main():
	infile = 'InteractionsDetectionMethods.tsv'
	evidencefile = 'mi.owl'

	## parse evidence file into a key of (MI,name) pairs.
	evmap = {}
	mi=None
	with open(evidencefile) as fin:
		for line in fin:
			if 'id: MI:' in line:
				mi = line.strip().split(' ')[1]
			elif 'name: ' in line and mi:
				name = '-'.join(line.strip().split(' ')[1:])
				evmap[mi]=name
				mi = None

	print '%d MI terms mapped.' % (len(evmap))

	## generate dictionary of (key,value) pairs, where key is a tuple
	## that represents an edge and the value is a set of evidence names.
	rows = {}
	with open(infile) as fin:
		for line in fin:
			row = line.strip().split()
			e = tuple([row[0],row[1]])
			ev = evmap.get(row[2],row[2])
			if e not in rows:
				rows[e] = set()
			rows[e].add(ev)
	print '%d edges read' % (len(rows))

	outfile = 'myProteinNet-flybase.txt'
	out = open(outfile,'w')
	out.write('#FlyBase1\tFlyBase2\tPubMedIDs\tEvidence\n')
	for e in rows:
		out.write('%s\t%s\tNA\t%s\n' % (e[0],e[1],';'.join(rows[e])))
	out.close()
	print 'wrote to %s' % (outfile)



if __name__ == '__main__':
	main()