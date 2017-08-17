import sys
'''
Parses the mentha file by re-ordering columns.
'''

def main():
	filename = '7227'
	rows = []
	with open(filename) as fin:
		for line in fin:
			if 'Protein A' in line: # skip header
				continue
			if line.strip() == '': # skip empty lines
				continue
			rows.append(line.strip().split(';')) # Protein A;Gene A;Protein B;Gene B;Score;PMID
			rows[-1][5] = [pmid for pmid in rows[-1][5].split(' ') if 'unassigned' not in pmid]

	outname = 'mentha-uniprot.txt' 
	out = open(outname,'w')
	out.write('#UniProtID1\tUniProtID2\tPubMedIds\tName1\tName2\tScore\n')
	for r in rows:
		out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % 
			(r[0],r[2],';'.join(r[5]),r[1],r[3],r[4]))
	out.close()
	print 'Wrote %d lines to %s' % (len(rows),outname)

if __name__ == '__main__':
	main()