from __future__ import print_function
import sys
import os
## NetworkX imports
import numpy as np
import networkx as nx			   
from graph import Graph
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

A_PARAMS = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 3]

W_PARAMS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

MIN_INTERACTIONS = 50

# skip interolog mapping
# FB2015_05 is the same as flybase_ppi
# hub-spoke-model is the same as anti-tag-coip & perrimon_coapcomplex
FILTERS = ['interolog','FB2015_05','hub-spoke-model','perrimon_coapcomplex']

def main(network,weighted_network,probfile,ground_truth_db,force):
	

	print('\nReading the network from %s' %(network))
	net = Graph()
	net.read(network,strip_db=False,collapsed=True,)
	print(nx.info(net))
	print('\nResetting original network weights')
	with open(weighted_network) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			net[row[0]][row[1]]['weight'] = float(row[2])
	
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
		if ev_type != 'droid:MI:0892':
			print('skipping ',ev_type)
			continue
		# get edges with that attribute:
		edges_with_ev = [(u,v) for (u,v,attrs) in net.edges(data=True) if ev_type in attrs['types']]
		mult_ev_edges = [e for e in edges_with_ev if len(net[e[0]][e[1]]['types']) > 1]
		if len(mult_ev_edges) >= MIN_INTERACTIONS:
			ev_types_filtered[ev_type] = mult_ev_edges
			print('  %d edges with evidence type "%s"; %d have mult evidence types' % (len(edges_with_ev),ev_type,len(mult_ev_edges)))
	print('%d %s evidence types that passed filtering' % (len(ev_types_filtered),ground_truth_db))


	test_dir = 'param-sweep-'+ground_truth_db+'/'
	if not os.path.isdir(test_dir):

		print('making directory')
		print('mkdir %s' % (test_dir))
		val = os.system('mkdir %s' % (test_dir))
		if val != 0:
			sys.exit()
	
	scores = {} # scores[ev_type][a1][a2][w1][w2] = Jaccard overlap of top quartile.
	for ev_type in ev_types_filtered:
		edges_to_keep = set(ev_types_filtered[ev_type])
		scores[ev_type] = {}
		

		interactome_prefix =  test_dir+'/'+ev_type.replace(':','-')
		if not os.path.isdir(interactome_prefix):
			print('mkdir %s' % (interactome_prefix))
			val = os.system('mkdir %s' % (interactome_prefix))
			if val != 0:
				sys.exit()
			
		# Generate interactome, removing ev_type from the list.
		# TODO (long term) -- remove pubmed IDs if we can link them back to evtype?
		interactome_file = interactome_prefix+'/interactome.txt'
		if os.path.isfile(interactome_file) and not force:
			print('File "%s" exists -- not overwriting.' % (interactome_file))
		else:
			out = open(interactome_file,'w')
			out.write('#symbol1\tsymbol2\tPubMedIDs\tFlyBase1\tFlyBase2\tDBs\tEvidence\n')
			for (u,v) in edges_to_keep:
				id1 = net.node[u]['id']
				id2 = net.node[v]['id']
				attrs = net[u][v]
				remaining_evs = [ev for ev in attrs['types'] if ev != ev_type]
				if len(remaining_evs) == 0:
					continue
				remaining_dbs = list(set([e.split(':')[0] for e in remaining_evs]))
				out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % \
					(u,v,';'.join(attrs['pmids']), id1,id2,';'.join(remaining_dbs), \
						';'.join(remaining_evs)
					))
			out.close()
			print('Wrote interactome file to "%s"' % (interactome_file))

		## Generate evidence-weighted probabilities. Use --prob value.
		weighted_interactome_prefix = interactome_prefix+'/interactome-weighted'
		weighted_interactome_file = weighted_interactome_prefix+'.txt'
		if os.path.isfile(weighted_interactome_file) and not force:
			print('File "%s" exists -- not overwriting.' % (weighted_interactome_file))
		else:
			cmd = 'python3 weight-edges-by-evidence.py -n %s --probs %s -o %s' % \
				(interactome_file,probfile,weighted_interactome_prefix)
			print(cmd)
			val = os.system(cmd)
			if val != 0:
				sys.exit()
			
		# Run parameter sweeps
		w_str = ''
		for w in W_PARAMS:
			w_str += ' --w1 %.3f --w2 %.3f' % (w,w)		
		i=1
		print('\nWeighting edges...')
		for a1 in A_PARAMS:
			scores[ev_type][a1] = {}
			for a2 in A_PARAMS:

				out_prefix = interactome_prefix+'/a1_%.3f_a2_%.3f' % (a1,a2)
				meta_information = out_prefix+'_filenames.txt'
				if os.path.isfile(meta_information) and not force:
					print('File "%s" already exists -- not overwriting.' % (meta_information))
					continue

				print('\n%d of %d: A1 = %f, A2 = %f' % (i,len(A_PARAMS)**2,a1,a2))
				i+=1
				cmd = 'python3 weight-edges.py -n %s -c -e %s -o %s --a1 %.3f --a2 %.3f %s\n' % \
							(interactome_file,weighted_interactome_file,out_prefix,a1,a2,w_str)
				print(cmd)
				val = os.system(cmd)
				if val != 0:
					sys.exit()

		print('\nCalculating scores...')
		quartile = int(nx.number_of_edges(net)/4)
		print('There are %d edges in the first quartile.' % (quartile))
		for a1 in A_PARAMS:
			scores[ev_type][a1] = {}
			for a2 in A_PARAMS:
				out_prefix = interactome_prefix+'/a1_%.3f_a2_%.3f' % (a1,a2)
				score_file = out_prefix+'_scores.txt'
				if os.path.isfile(score_file) and not force:
					print('File "%s" already exists -- not overwriting.' % (score_file))
					continue

				scores[ev_type][a1][a2] = {}
				meta_information = out_prefix+'_filenames.txt'
				with open(meta_information) as fin:

					for line in fin:
						if line[0] == '#':
							continue
						row = line.strip().split()
						w1 = float(row[2])
						w2 = float(row[3])
						fname = row[5]
						print('  reading %s'  %(fname))

						if w1 not in scores[ev_type][a1][a2]:
							scores[ev_type][a1][a2][w1] = {}

						# get scores for all edges.  Start with ones that
						# are affected by ev_type disappearing.
						score_dict = {}
						with open(fname) as efin:
							for eline in efin:
								if eline[0] == '#':
									continue
								erow = eline.strip().split()
								score_dict[tuple(sorted([erow[0],erow[1]]))] = float(erow[2])
						# now add remainder of edges from original network; 
						# confirm that these are not in edges_to_keep.
						for u,v in net.edges():
							if (u,v) not in score_dict:
								if tuple(sorted([u,v])) in edges_to_keep:
									sys.exit('ERROR: (%s,%s) is a kept edge but not yet in score dictionary.' % (u,v))
								score_dict[tuple(sorted([u,v]))] = net[u][v]['weight']
								
						sorted_scores = sorted(score_dict.items(),key=lambda x: x[1], reverse=True)
						top_edges = set([e for e,s in sorted_scores[:quartile]])
						
						# compute jaccard overlap
						scores[ev_type][a1][a2][w1][w2] = float(len(top_edges.intersection(edges_to_keep)))/min(len(top_edges),len(edges_to_keep))
						#print('PARAMS:',a1,a2,w1,w2,'EDGES:',len(edges_to_keep),len(top_edges),len(top_edges.intersection(edges_to_keep)),'SCORE:',scores[ev_type][a1][a2][w1][w2])

				out = open(score_file,'w')
				out.write('\t'+'\t'.join(['W2=%.3f'%w for w in W_PARAMS])+'\n')
				for w1 in W_PARAMS:
					score_str = '\t'.join([str(scores[ev_type][a1][a2][w1].get(w2,'NaN')) for w2 in W_PARAMS])
					out.write('W1=%.3f\t' % (w1)+score_str+'\n')

				out.close()
				print('wrote scores to "%s"' % (score_file))

		## make heatmap just for funsies.
		figfile = interactome_prefix+'/%s_score_heatmap.png' % (ev_type.replace(':','-'))
		if os.path.isfile(figfile) and not force:
			print('File "%s" already exists -- not overwriting.' % (figfile))
		else:
			print('\nMaking a heatmap (just for fun)')
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
			
			fig.set_tight_layout(True)
			fig.savefig(figfile)
			print('wrote to %s' % (figfile))		
			#plt.tight_layout()
			#plt.savefig(figfile)

	return


if __name__ == '__main__':
	if len(sys.argv) != 5 and len(sys.argv) != 6:
		sys.exit('USAGE: python3 parameter_sweep.py <INTERACTOME> <WEIGHTED_INTERACTOME> <ETYPE_PROBS> <DB> <FORCE-optional>')
	network = sys.argv[1]
	weighted_network = sys.argv[2]
	probfile = sys.argv[3]
	ground_truth_db = sys.argv[4]
	if len(sys.argv)==6:
		force=True
	else:
		force=False
	main(network,weighted_network,probfile,ground_truth_db,force)