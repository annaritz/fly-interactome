import sys
'''
Rearranges the columns to generate FlyBase IDs.
'''

def main():
	infile = 'flymine-query.txt'
	lines = {}
	with open(infile) as fin:
		for line in fin:
			if '#' in line: # skip header
				continue
			row = line.strip().split('\t') 
			e = tuple([row[0],row[1]])
			if e not in lines:
				lines[e] = [set(),set()] # pubmedids, name
			lines[e][0].add(row[2])
			lines[e][1].add(row[3])
	print '%d edges read' % (len(lines))


	outfile = 'flymine-flybase.txt'
	out = open(outfile,'w')
	out.write('#FlyBase1\tFlyBase2\tPubMedIDs\tInteractionName\n')
	for e in lines:
		# replace any spaces with dashes
		lines[e][0] = [s.replace(' ','-') for s in lines[e][0]]
		lines[e][1] = [s.replace(' ','-') for s in lines[e][1]]
		out.write('%s\t%s\t%s\t%s\n' % (e[0],e[1],';'.join(lines[e][0]),';'.join(lines[e][1])))
	out.close()
	print 'wrote to %s' % (outfile)

if __name__ == '__main__':
	main()