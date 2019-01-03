from __future__ import print_function
import sys
import os
## NetworkX imports
import numpy as np
import networkx as nx			   
from graph import Graph
import matplotlib
import matplotlib.pyplot as plt

A_PARAMS = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 3]
W_PARAMS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

MIN_INTERACTIONS = 50

# skip interolog mapping
# FB2015_05 is the same as flybase_ppi
# hub-spoke-model is the same as anti-tag-coip & perrimon_coapcomplex
FILTERS = ['interolog','FB2015_05','hub-spoke-model','perrimon_coapcomplex']

def main(network,probfile,ground_truth_db):
	

	print('\nReading the network from %s' %(network))
	net = Graph()
	net.read(network,strip_db=False,collapsed=True,)
	print(nx.info(net))
	
	# get all evidence types for that db:
	ev_types = set()
	for (u,v,attrs) in net.edges(data=True):
		for ev_type in attrs['types']:
			if ground_truth_db in ev_type and not any([f in ev_type for f in FILTERS]):
				ev_types.add(ev_type)
	print('%d %s evidence types' % (len(ev_types),ground_truth_db))

	# filter those that don't have enough edges after identifying those with mutliple edges.
	ev_types_filtered = {}
	for ev_type in ev_types:
		# get edges with that attribute:
		edges_with_ev = [(u,v) for (u,v,attrs) in net.edges(data=True) if ev_type in attrs['types']]
		mult_ev_edges = [e for e in edges_with_ev if len(net[e[0]][e[1]]['types']) > 1]
		if len(mult_ev_edges) >= MIN_INTERACTIONS:
			ev_types_filtered[ev_type] = mult_ev_edges
			print('  %d edges with evidence type "%s"; %d have mult evidence types' % (len(edges_with_ev),ev_type,len(mult_ev_edges)))
	print('%d %s evidence types that passed filtering' % (len(ev_types_filtered),ground_truth_db))

	test_dir = 'param-sweep-'+ground_truth_db
	if not os.path.isdir(test_dir):
		print('making directory')
		print('mkdir %s' % (test_dir))
		os.system('mkdir %s' % (test_dir))
	else:
		print('NOT clearing directory.')
		#print('clearing directory')
		#print('rm -fr %s/*' % (test_dir))
		#os.system('rm -fr %s/*' % (test_dir))

	scores = {} # scores[ev_type][a1][a2][w1][w2] = Jaccard overlap of top quartile.
	for ev_type in ev_types_filtered:
		edges_to_keep = set(ev_types_filtered[ev_type])
		scores[ev_type] = {}
		

		interactome_prefix =  test_dir+'/'+ev_type.replace(':','-')
		print('mkdir %s' % (interactome_prefix))
		os.system('mkdir %s' % (interactome_prefix))
		
		# Generate interactome with JUST the edges in edges_to_keep 
		# with all the evidence types EXCEPT ev_type.
		# TODO: remove ground_truth_db if there are no other evtypes from that db
		# TODO (long term) -- remove pubmed IDs if we can link them back to evtype?
		interactome_file = interactome_prefix+'/interactome.txt'
		out = open(interactome_file,'w')
		out.write('#symbol1\tsymbol2\tPubMedIDs\tFlyBase1\tFlyBase2\tDBs\tEvidence\n')
		for (u,v) in edges_to_keep:
			attrs = net[u][v]
			id1 = net.node[u]['id']
			id2 = net.node[v]['id']

			out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
				(u,v,';'.join(attrs['pmids']), id1,id2,';'.join(attrs['dbs']), \
					';'.join([ev for ev in attrs['types'] if ev != ev_type])
				))
		out.close()
		print('Wrote interactome file to "%s"' % (interactome_file))

		## Generate evidence-weighted probabilities. Use --prob value.
		weighted_interactome_prefix = interactome_prefix+'/interactome-weighted'
		weighted_interactome_file = weighted_interactome_prefix+'.txt'
		cmd = 'python weight-edges-by-evidence.py -n %s --probs %s -o %s' % \
			(interactome_file,probfile,weighted_interactome_prefix)
		print(cmd)
		os.system(cmd)
		
		# Run parameter sweeps
		w_str = ''
		for w in W_PARAMS:
			w_str += ' --w1 %.3f --w2 %.3f' % (w,w)		
		i=1
		for a1 in A_PARAMS:
			scores[ev_type][a1] = {}
			for a2 in A_PARAMS:
				scores[ev_type][a1][a2] = {}

				out_prefix = interactome_prefix+'/a1_%.3f_a2_%.3f' % (a1,a2)
				score_file = out_prefix+'_scores.txt'
				if os.path.isfile(score_file):
					print('SKIPPING MAKING THIS FILE AGAIN: "%s"' % (score_file))
					continue

				print('\n%d of %d: A1 = %f, A2 = %f' % (i,len(A_PARAMS)**2,a1,a2))
				i+=1
				cmd = 'python weight-edges.py -n %s -c -e %s -o %s --a1 %.3f --a2 %.3f %s\n' % \
							(interactome_file,weighted_interactome_file,out_prefix,a1,a2,w_str)
				print(cmd)
				os.system(cmd)

				meta_information = out_prefix+'_filenames.txt'
				with open(meta_information) as fin:
					for line in fin:
						if line[0] == '#':
							continue
						row = line.strip().split()
						w1 = float(row[2])
						w2 = float(row[3])
						fname = row[5]

						if w1 not in scores[ev_type][a1][a2]:
							scores[ev_type][a1][a2][w1] = {}

						# compute average score for these ground truth edges.
						scores[ev_type][a1][a2][w1][w2] = 0
						num_edges = 0
						with open(fname) as efin:
							for eline in efin:
								if eline[0] == '#':
									continue
								erow = eline.strip().split()
								scores[ev_type][a1][a2][w1][w2] += float(erow[2])
								num_edges+=1
						scores[ev_type][a1][a2][w1][w2] /= num_edges

				out = open(score_file,'w')
				out.write('\t'+'\t'.join(['W2=%.3f'%w for w in W_PARAMS])+'\n')
				for w1 in W_PARAMS:
					score_str = '\t'.join([str(scores[ev_type][a1][a2][w1].get(w2,'NaN')) for w2 in W_PARAMS])
					out.write('W1=%.3f\t' % (w1)+score_str+'\n')

				out.close()
				print('wrote scores to "%s"' % (score_file))

		## make heatmap just for funsies.
		fig, ax_array = plt.subplots(nrows=len(A_PARAMS),ncols=len(A_PARAMS),figsize=(40,40))
		current_cmap = matplotlib.cm.get_cmap()
		current_cmap.set_bad(color='black')
		for i in range(len(A_PARAMS)):
			a1 = A_PARAMS[i]
			for j in range(len(A_PARAMS)):
				a2 = A_PARAMS[j]

				out_prefix = interactome_prefix+'/a1_%.3f_a2_%.3f' % (a1,a2)
				score_file = out_prefix+'_scores.txt'

				M = []
				with open(score_file) as fin:
					for line in fin:
						M.append(line.strip().split())
				M = M[1:] # strip first row
				M = [m[1:] for m in M] # strip first col
				for mi in range(len(M)):
					for mj in range(len(M)):
						if M[mi][mj] == 'NaN':
							M[mi][mj] = np.nan
						else:
							M[mi][mj] = float(M[mi][mj])
				ax = ax_array[i][j]
				im = ax.imshow(M,interpolation='none',origin='lower',vmin=0,vmax=1)
				ax.set_title('A1=%.3f A2=%.3f\n(max=%.3f)' % (a1,a2,max([max(m) for m in M])))
				ax.set_ylabel('W1')
				ax.set_yticklabels(W_PARAMS)
				ax.set_yticks(range(len(W_PARAMS)))
				ax.set_xlabel('W2')
				ax.set_xticks(range(len(W_PARAMS)))
				ax.set_xticklabels(W_PARAMS)

				fig.colorbar(im, ax=ax)
				
		plt.tight_layout()
		plt.savefig(interactome_prefix+'/score_heatmap.png')

	return


if __name__ == '__main__':
	if len(sys.argv) != 4:
		sys.exit('USAGE: python parameter_sweep.py <INTERACTOME> <ETYPE_PROBS> <DB>')
	network = sys.argv[1]
	probfile = sys.argv[2]
	ground_truth_db = sys.argv[3]
	main(network,probfile,ground_truth_db)