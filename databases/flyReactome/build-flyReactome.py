import sys
'''
Rearranges columns in the FlyReactome tab-delimited file.  
Keeps them in UniProt IDs.
'''

def main():
	infile = 'drosophila_melanogaster.interactions.txt'

	## lines is a dictionary of (key,value) pairs, where the keys are tuples
	## representing edges and the values are lists of (a) pubmedids and (b) 
	## interaction type.
	lines = {}
	with open(infile) as fin:
		for line in fin:
			if '#' in line: # skip header
				continue
			row = line.strip().split('\t') 
			if len(row) < 9: # require PubMed IDs
				continue
			print row
			node1 = row[0].split(':')[1]
			node2 = row[3].split(':')[1]
			e = tuple([node1,node2])
			if e not in lines:
				lines[e] = [set(),set()] # pubmedids, type
			lines[e][0].update(row[8].split(','))
			lines[e][1].update(row[6].split(','))
	print '%d edges read' % (len(lines))


	outfile = 'flyReactome-uniprot.txt'
	out = open(outfile,'w')
	out.write('#UniProt1\tUniProt2\tPubMedIDs\tType\n')
	for e in lines:
		out.write('%s\t%s\t%s\t%s\n' % (e[0],e[1],';'.join(lines[e][0]),';'.join(lines[e][1])))
	out.close()
	print 'wrote to %s' % (outfile)

if __name__ == '__main__':
	main()
