import sys
import mygene
import argparse
'''Utility script to generate node file with mappings from edgefile.'''

def parse_args():
	parser = argparse.ArgumentParser(description='Map Edgefile to Nodes')
	parser.add_argument('infile', metavar='INFILE', type=str, 
                    help='Input edge file with Node1,Node2 as first two columns.')
	parser.add_argument('outfile', metavar='OUTFILE', type=str, 
                    help='Output file to write to.')
	parser.add_argument('--mapfrom', dest='mapfrom', type=str, default='uniprot',
                    help='Map from namespace (default=uniprot)')
	parser.add_argument('--mapto', dest='mapto', type=str, default='FLYBASE',
                    help='Map to namespace (default=FLYBASE)')
	parser.add_argument('--species',dest='species',type=str,default='7227',
					help='Species taxon ID (default=7227).')
	args = parser.parse_args()
	print args
	return args

def main():
	## parse arguments
	args = parse_args()
	
	## read original file
	ids = set()
	with open(args.infile) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			ids.update(line.strip().split('\t')[:2])
	
	## get gene IDs
	mg = mygene.MyGeneInfo()
	result = mg.querymany(ids,scopes=args.mapfrom,fields=[args.mapto,'symbol','alias'],species=args.species)

	## write mapped file
	out = open(args.outfile,'w')
	out.write('%s\t%s\t%s\t%s\n' % (args.mapfrom,args.mapto,'Symbol','Aliases'))

	missing = 0
	tot=0
	for r in result:
		tot+=1
		if 'notfound' in r and r['notfound']==True or args.mapto not in r or 'symbol' not in r or 'alias' not in r:
			# skip if query wasn't mapped.
			missing+=1
			continue
		mapfrom = r['query']
		mapto = r[args.mapto]
		symbol = r['symbol']
		alias = r['alias']
		#print(mapfrom,mapto,alias)
		if isinstance(mapto,dict):
			if args.mapto == 'uniprot': 
				if 'Swiss-Prot' in mapto: # take reviewed (swiss-prot) ids
					mapto = mapto['Swiss-Prot']
				elif 'TrEMBL' in mapto: # otherwise, take unreviewed (TrEMBL) ids
					mapto = mapto['TrEMBL']
				else:
					sys.exit('Error: namespace %s has new options: %s\n', \
						args.mapto,';'.join(mapto.keys()))
			else:
				sys.exit('Error: namespace %s has multiple options: %s\n', \
					args.mapto,';'.join(mapto.keys()))
		if not isinstance(mapto,list):
			mapto = [mapto]
		if not isinstance(symbol,list):
			symbol = [symbol]
		if not isinstance(alias,list):
			alias = [alias]
		out.write('%s\t%s\t%s\t%s\n' % (mapfrom,';'.join(mapto),';'.join(symbol),';'.join(alias)))
		
	out.close()

	print '%d nodes in %s'  % (len(ids),args.infile)
	print '%d lines written to %s' % (tot,args.outfile)
	print '%d lines were missing an id' % (missing)

if __name__ == '__main__':
	main()