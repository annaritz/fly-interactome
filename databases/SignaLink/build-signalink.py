import sys
'''
Parse the SignaLink database, downloaded as described in the README
and saved as signalink-download.csv.
'''

def main(): 
	infile = 'signalink-download.csv'

	## unilines and flylines are dictionaries of (key,value) pairs, where the keys
	## are tuples representing edges and the values are lists of three sets:
	## (1) pubmedids
	## (2) interaction type
	## (3) directedness
	unilines = {}
	flylines = {}
	with open(infile) as fin:
		for line in fin:
			if 'source_name' in line: # skip header
				continue
			row = line.strip().split(';')
			uniedge = tuple([row[1],row[7]])
			flyedge = tuple([row[2],row[8]])
			if uniedge not in unilines:
				unilines[uniedge] = [set(),set(),set()] # pubmedids, interactiontype, directedness
			if flyedge not in flylines:
				flylines[flyedge] = [set(),set(),set()] # pubmedids, interactiontype, directedness

			pubmedupdate = set(row[16].split('|'))
			interactiontypeupdate = row[12].replace(' ','_')
			directednessupdate = row[13].replace(' ','_')

			unilines[uniedge][0].update(pubmedupdate)
			unilines[uniedge][1].add(interactiontypeupdate)
			unilines[uniedge][2].add(directednessupdate)

			flylines[flyedge][0].update(pubmedupdate)
			flylines[flyedge][1].add(interactiontypeupdate)
			flylines[flyedge][2].add(directednessupdate)

	print '%d uniprot edges and %d flybase edges' % (len(unilines),len(flylines))

	uniprotfile = 'signalink-uniprot.txt'
	flybasefile = 'signalink-flybase.txt'

	uniout = open(uniprotfile,'w')
	flyout = open(flybasefile,'w')
	uniout.write('#UniProt1\tUniProt2\tPubMedIDs\tInteractionType\tDirectedness\n')
	flyout.write('#FlyBase1\tFlyBase2\tPubMedIDs\tInteractionType\tDirectedness\n')

	for e in unilines:
		uniout.write('%s\t%s\t%s\t%s\t%s\n' % (e[0],e[1],';'.join(unilines[e][0]),';'.join(unilines[e][1]),';'.join(unilines[e][2])))
	for e in flylines:
		flyout.write('%s\t%s\t%s\t%s\t%s\n' % (e[0],e[1],';'.join(flylines[e][0]),';'.join(flylines[e][1]),';'.join(flylines[e][2])))
		
	uniout.close()
	flyout.close()
	print 'Wrote to %s and %s' % (uniprotfile,flybasefile)

if __name__ == '__main__':
	main()
