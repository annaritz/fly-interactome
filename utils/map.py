import sys
import mygene
import argparse
'''
Utility script that maps different namespaces, using the MyGene.info REST web service.
Used to convert FlyBase IDs to UniProtKB IDs and vice versa.
'''

def parse_args():
	parser = argparse.ArgumentParser(description='Map Edgefile to Different Namespace')
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
	parser.add_argument('--pubmedcolid',dest='colid',type=int,default=2,
					help='0-indexed column id of PubMed evidence (default=2).')
	parser.add_argument('--singlecol',dest='singlecol',action='store_true',default=False,
					help='Input file is single column of identifiers; just map these (default=False).')
	parser.add_argument('--retain_orig',dest='retain',action='store_true',default=False,
					help='Include original identifier when using --singlecol mode (default=False).')
	parser.add_argument('--merge_ids',dest='merge',action='store_true',default=False,
					help='If multiple IDs exist in the mapped namespace, merge them into one identifier (default=False).')
	args = parser.parse_args()
	print args
	return args

def mapSingleColumn(args):
	## read original file
	ids = set([line.strip() for line in open(args.infile).readlines()])
	print 'Querying %d total IDs:' % (len(ids))

	## get gene IDs
	mg = mygene.MyGeneInfo()
	result = mg.querymany(ids,scopes=args.mapfrom,fields=args.mapto,species=args.species)

	## build mapped dictionary.
	mapped = build_mapped_dict(result,args)
		
	## write mapped file
	out = open(args.outfile,'w')
	missing = 0
	tot=0
	for thisid in ids:
		if thisid not in mapped:
			missing += 1
			continue
		nodes = mapped[thisid]
		for node in nodes:
			if args.retain_orig:
				out.write('%s\t%s\n' % (thisid,node))
			else:
				out.write('%s\n' % (node))
			tot+=1

	print '%d ids in original file %s'  % (len(ids),args.infile)
	print '%d ids written to %s' % (tot,args.outfile)
	print '%d ids were missing' % (missing)
	return

def build_mapped_dict(result,args):
	'''
	Build mapped dictionary.
	'''

	## mapped is a dictionary of (key,value) pairs, where the key is the original 
	## namespace and the value is a set of IDs in the new mapped namespace.
	mapped = {}  
	for r in result:
		mapfrom = r['query']
		if 'notfound' in r and r['notfound']==True or args.mapto not in r:
			# skip if query wasn't mapped.
			continue
		else:
			mapto = r[args.mapto]
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
		if mapfrom not in mapped:
			mapped[mapfrom] = set()

		if isinstance(mapto,list):
			mapped[mapfrom].update(set([m for m in mapto]))
		else:
			mapped[mapfrom].add(mapto)
	return mapped

def main():
	## parse arguments
	args = parse_args()

	if args.singlecol:
		mapSingleColumn(args)
		return

	## read original file
	lines = []
	with open(args.infile) as fin:
		for line in fin:
			if line[0] == '#':
				header = line[1:].strip().split('\t')
			else:
				lines.append(line.strip().split('\t'))
	ids = set([r[0] for r in lines]).union(set([r[1] for r in lines]))
	print 'Querying %d total IDs:' % (len(ids))

	## get gene IDs
	mg = mygene.MyGeneInfo()
	result = mg.querymany(ids,scopes=args.mapfrom,fields=args.mapto,species=args.species)

	## build mapped dictionary.
	mapped = build_mapped_dict(result,args)

	if args.merge: # redefine mapped dictionary to be single-element items
		for key in mapped:
			# remove the ':CGXXXX' from the names (use common name only)
			mapped[key] = set([val.split(':')[0] for val in mapped[key]])
			# join the names into one string.
			mapped[key] = set([';'.join(sorted(mapped[key]))])
	
	## write mapped file
	out = open(args.outfile,'w')
	out.write('#%s1\t%s2\t%s\t%s\n' %  \
			(args.mapto,args.mapto,header[args.colid], \
			'\t'.join([header[i] for i in range(len(header)) if i != args.colid])))

	to_write = {}
	missing = 0
	tot=0
	for line in lines:
		if line[0] not in mapped or line[1] not in mapped:
			missing += 1
			continue
		nodes1 = mapped[line[0]]
		nodes2 = mapped[line[1]]
		for node1 in nodes1:
			for node2 in nodes2:
				# sort nodes since the graph is undirected
				u,v = sorted([node1,node2])

				if u not in to_write:
					to_write[u] = {}
				if v not in to_write[u]:
					if u == node1: # if sorting didn't make a difference, add as-is.
						to_write[u][v] = [e.replace(' ','-') for e in line]
					else: # sorting flipped nodes; flip first two entries.
						to_write[u][v] = [line[1].replace(' ','-'),line[0].replace(' ','-')] + \
							[e.replace(' ','-') for e in line[2:]]

				else: # append to the existing edge.
					if u == node1: # if sorting didn't make a difference, add as-is.
						for i in range(len(line)):
							if line[i] not in to_write[u][v][i]:
								to_write[u][v][i] += ';'+line[i].replace(' ','-')
					else: # sorting flipped nodes; flip first two entries.
						if line[1] not in to_write[u][v][0]:
							to_write[u][v][0] += ';'+line[1].replace(' ','-')
						if line[0] not in to_write[u][v][1]:
							to_write[u][v][1] += ';'+line[0].replace(' ','-')
						for i in range(2,len(line)):
							if line[i] not in to_write[u][v][i]:
								to_write[u][v][i] += ';'+line[i].replace(' ','-')

	# go through to_write dictionary and write the lines.
	for node1 in sorted(to_write.keys()):
		for node2 in sorted(to_write[node1].keys()):
			# remove redundant entries (these can occur because we combined multiple edges that represented the same interaction)
			for i in range(len(to_write[node1][node2])):
				to_write[node1][node2][i] = ';'.join(sorted(list(set(to_write[node1][node2][i].split(';')))))
				if to_write[node1][node2][i] == '':
					#print('EMPTY LINE')
					to_write[node1][node2][i] = 'NaN'

			## double check that the list is 4 elements long - die otherwise.
			if len(to_write[node1][node2]) < 4:
				print(node1,node2,to_write[node1][node2])
				sys.exit()

			# finally, write the line. 
			out.write('%s\t%s\t%s\t%s\n' %  \
				(node1,node2,to_write[node1][node2][args.colid], \
				'\t'.join([to_write[node1][node2][i] for i in range(len(to_write[node1][node2])) if i != args.colid])))
			
			tot+=1

	print '%d lines in original file %s'  % (len(lines),args.infile)
	print '%d lines written to %s' % (tot,args.outfile)
	print '%d lines were missing an id' % (missing)

if __name__ == '__main__':
	main()