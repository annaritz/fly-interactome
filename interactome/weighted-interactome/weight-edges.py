from __future__ import print_function
from optparse import OptionParser   
import sys						  
from math import exp
## NetworkX imports
import networkx as nx			   
from graph import Graph



def main(args):
	# parse arguments
	opts = parse_arguments(args)

	MIN_WEIGHT = 0.01 # minimum weight w1, w2, w3

	print('\nReading the network from %s' %(opts.network))
	net = Graph()
	net.read(opts.network,opts.collapsed)
	print(nx.info(net))

	print('\nReading the evidence-weighted edges from %s' % (opts.evidence))
	evidence_weights = read_evidence_weights(opts.evidence)
	
	meta_out = open(opts.outprefix+'_filenames.txt','w')
	meta_out.write('#A1\tA2\tW1\tW2\tW3\tFilename\n')
	for a1 in opts.a1:
		for a2 in opts.a2:
			for w1 in opts.w1:
				for w2 in opts.w2:
					
					## make sure this param combo is OK
					if a1 < 0 or a2 < 0:
						print('\nWARNING: skipping non-positive steepness parameter combo A1=%.f and A2=%f' % (a1,a2))
						continue
					w3 = 1-w1-w2
					
					if w1 < MIN_WEIGHT or w1 > 1 or w2 < MIN_WEIGHT or w2 > 1 or w3 < MIN_WEIGHT or w3 > 1:
						# this combination of W1/W2/W3 not OK; ignore.
						continue
					
					outfile = opts.outprefix+'_w1_%.3f_w2_%.3f.txt' % (w1,w2)
					meta_out.write('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%s\n' % (a1,a2,w1,w2,w3,outfile))
					out = open(outfile,'w')
					out.write('#tail\thead\tedge_weight\n')
					for (u,v,attrs) in net.edges(data=True):
						num_pmids = len([a for a in attrs['pmids'] if a != 'unpublished'])
						num_dbs = len(attrs['dbs'])
						if (u,v) in evidence_weights:
							ev_weight = evidence_weights[(u,v)]
						else:
							ev_weight = evidence_weights[(v,u)]
						final_weight = w1*non_linear_saturating_function(num_pmids,a1) + \
									w2 * non_linear_saturating_function(num_dbs,a2) + \
									w3 * ev_weight
						out.write('%s\t%s\t%e\n' % (u,v,final_weight))
					out.flush()
					out.close()
					#print('Wrote to %s' % (outfile))
	meta_out.close()
	print('Wrote meta information to "%s"' % (opts.outprefix+'_filenames.txt'))
	return

def non_linear_saturating_function(n,a):
	## from HIPPIE
	## https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0031826
	return float(2)/(1+exp(-a*n))-1

def read_evidence_weights(evidence):
	evidence_weights = {}
	with open(evidence) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			evidence_weights[(row[0],row[1])] = float(row[2])
	print('Read evidence weights for %d edges' % (len(evidence_weights)))
	return evidence_weights

def parse_arguments(args):
	usage = "python weight-edges.py [OPTIONS]"
	parser = OptionParser(usage=usage)
	parser.add_option('-n','--network',type='string',metavar='STR',\
						help='protein-protein interactome (tab-delimited). Required.')
	parser.add_option('-c','--collapsed',action='store_true',default=False,
						help='Interactome is collapsed by common name (has 7 columns).')
	parser.add_option('-o','--outprefix',type='string',metavar='STR',\
						help='output file prefix. Required.')
	parser.add_option('-e','--evidence',type='string',metavar='STR',\
						help='Evidence weights (from weight-edges-by-evidence.py)')
	parser.add_option('','--a1',type='float',metavar='FLOAT',action='append', default=[], \
						help='Parameter that controls the steepness of s1 (study-based term). Can pass multiple values. Default = 1.')
	parser.add_option('','--a2',type='float',metavar='FLOAT',action='append', default=[], \
						help='Parameter that controls the steepness of s2 (db-based term). Can pass multiple values. Default = 1')
	parser.add_option('','--w1',type='float',metavar='FLOAT',action='append', default=[], \
						help='Weight of s1 (study-based term).  Float between 0 and 1; w3 = 1-w1-w2. Can pass multiple values. Default = 0.33')
	parser.add_option('','--w2',type='float',metavar='FLOAT',action='append', default=[], \
						help='Weight of s2 (db-based term).  Float between 0 and 1; w3 = 1-w1-w2. Can pass multiple values. Default = 0.33')

	(opts, args) = parser.parse_args()
	print('ARGUMENTS ARE:',opts)

	if opts.network == None:
		parser.print_help()
		sys.exit('\nERROR: network file required')

	if opts.outprefix == None:
		parser.print_help()
		sys.exit('\nERROR: output filename prefix required')

	if opts.evidence == None:
		parser.print_help()
		sys.exit('\nERROR: evidence file required (generated from weight-edges-by-evidence.py).')


	# add default params here if not specified:
	if opts.a1 == []:
		opts.a1 = [1.0]
	if opts.a2 == []:
		opts.a2 = [1.0]
	if opts.w1 == []:
		opts.w1 = [0.33]
	if opts.w2 == []:
		opts.w2 = [0.33]
	
	return opts

if __name__ == '__main__':
	main(sys.argv)
