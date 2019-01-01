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

	print('\nReading the network from %s' %(opts.network))
	net = Graph()
	net.read(opts.network,opts.collapsed)
	print(nx.info(net))

	print('\nReading the evidence-weighted edges from %s' % (opts.evidence))
	evidence_weights = read_evidence_weights(opts.evidence)
	
	out = open(opts.outfile,'w')
	out.write('#tail\thead\tedge_weight\n')
	for (u,v,attrs) in net.edges(data=True):
		num_pmids = len([a for a in attrs['pmids'] if a != 'unpublished'])
		num_dbs = len(attrs['dbs'])

		final_weight = opts.w1*non_linear_saturating_function(num_pmids,opts.a1) + \
					opts.w2 * non_linear_saturating_function(num_dbs,opts.a2) + \
					opts.w3 * evidence_weights[(u,v)]
		out.write('%s\t%s\t%e\n' % (u,v,final_weight))
	out.flush()
	out.close()
	print('Wrote to %s' % (opts.outfile))
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
	parser.add_option('-o','--outfile',type='string',metavar='STR',\
						help='output file prefix. Required.')
	parser.add_option('-e','--evidence',type='string',metavar='STR',\
						help='Evidence weights (from weight-edges-by-evidence.py)')
	parser.add_option('','--a1',type='float',metavar='FLOAT',default=1.0, \
						help='Parameter that controls the steepness of s1 (study-based term). Default = 1.')
	parser.add_option('','--a2',type='float',metavar='FLOAT',default=1.0, \
						help='Parameter that controls the steepness of s2 (db-based term). Default = 1')
   	parser.add_option('','--w1',type='float',metavar='FLOAT',default=0.33, \
						help='Weight of s1 (study-based term).  Float between 0 and 1; w3 = 1-w1-w2. Default = 0.33')
   	parser.add_option('','--w2',type='float',metavar='FLOAT',default=0.33, \
						help='Weight of s2 (db-based term).  Float between 0 and 1; w3 = 1-w1-w2. Default = 0.33')

   	(opts, args) = parser.parse_args()
	print('ARGUMENTS ARE:',opts)

	if opts.network == None:
		parser.print_help()
		sys.exit('\nERROR: network file required')

	if opts.outfile == None:
		parser.print_help()
		sys.exit('\nERROR: output filename required')

	if opts.evidence == None:
		parser.print_help()
		sys.exit('\nERROR: evidence file required (generated from weight-edges-by-evidence.py).')

	if opts.a1 < 0 or opts.a2 < 0:
		parser.print_help()
		sys.exit('\nERROR: steepness parameters must be positive floats.')

	opts.w3 = 1.0-opts.w1-opts.w2
	if opts.w1 < 0 or opts.w1 > 1 or opts.w2 < 0 or opts.w2 > 1 or opts.w3 < 0 or opts.w3 > 1:
		parser.print_help()
		sys.exit('\nERROR: weights are not between 0 and 1 (remember that w3=1-w1-w2).')

	return opts

if __name__ == '__main__':
	main(sys.argv)
