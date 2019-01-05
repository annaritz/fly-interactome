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
import glob

A_PARAMS = [0.05, 0.1, 0.5, 1.0, 1.5, 2.0]
W_PARAMS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

def main(network,weighted_network,outprefix,force):

	print('\nRead All Information from param-sweep-droid/.')
	param_dir = 'param-sweep-droid/'
	ev_file = 'param-sweep-droid-evtypes.txt'
	ev_types = {} # {ev_type: (tot, overlap)}
	max_val = 0
	with open(ev_file) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			row[0] = row[0].replace(':','-')
			ev_types[row[0]] = (int(row[1]),int(row[2]))
			max_val = max(0,float(int(row[1])/int(row[2])))
	# normalize values to be between 0 and 1:
	ev_type_weights = {} # {ev_type: weight for f}
	for e in ev_types:
		ev_type_weights[e] = (float(ev_types[e][0])/float(ev_types[e][1]))/max_val
	
	print('%d evidence types' % (len(ev_types)))
	scores = {} # scores[ev_type][a1][a2]
	for ev_type in ev_types:
		print('  reading in %s (weight = %.4f)'  %(ev_type,ev_type_weights[ev_type]))
		scores[ev_type] = {}
		for a1 in A_PARAMS:
			scores[ev_type][a1] = {}
			for a2 in A_PARAMS:
				scorefile=param_dir+'/'+ev_type+'/a1_%.3f_a2_%.3f_scores.txt' % (a1,a2)
				scores[ev_type][a1][a2] = np.loadtxt(scorefile,skiprows=1,usecols=range(1,len(W_PARAMS)+1))
				#print(scorefile)
				#print(scores[ev_type][a1][a2])
				#max_val = max(max_val,np.nanmax(scores[ev_type][a1][a2]))
	
	print('\nCompute f for every parameter combination:')
	max_val = -10000000
	min_val = 0
	tot_opts = 0
	f = {} # f[a1][a2]
	for a1 in A_PARAMS:
		f[a1] = {}
		for a2 in A_PARAMS:
			f[a1][a2] = np.zeros((len(W_PARAMS),len(W_PARAMS)))
			for ev_type in ev_types:
				f[a1][a2] = np.add(f[a1][a2],ev_type_weights[ev_type]*np.log2(scores[ev_type][a1][a2]))
				tot_opts += len(W_PARAMS)*len(W_PARAMS)/2
				#print(ev_type)
				#print(f[a1][a2])
			max_val = max(max_val,np.nanmax(f[a1][a2]))
			min_val = min(min_val,np.nanmin(f[a1][a2]))
	print('Min Value =',(min_val))
	print('Max Value =',(max_val))

	param_combos = set() # set of (a1,a2,w1,w2) parameter combos that are optimal.
	for a1 in A_PARAMS:
		for a2 in A_PARAMS:
			for i in range(len(W_PARAMS)):
				for j in range(len(W_PARAMS)):
					if f[a1][a2][i][j] == max_val:
						print('OPTIMAL PARAM COMBO FOUNDS:',a1,a2,W_PARAMS[i],W_PARAMS[j])
						param_combos.add((a1,a2,W_PARAMS[i],W_PARAMS[j]))
	print ('%d parameter combinations out of %d are optimal.' % (len(param_combos),tot_opts))

	out_dir = outprefix+"-full-weights/"
	if not os.path.isdir(out_dir):
		val = os.system('mkdir %s'  %(out_dir))
		if val != 0:
			sys.exit()
	
	print('\nCalculate weights and compute the interquartile range for every combo')
	out = open(outprefix+'iqr.txt','w')
	out.write('#a1\ta2\tw1\tw2\tiqr\n')
	best_combo = None
	best_val = 0
	i = 0
	for a1,a2,w1,w2 in param_combos:
		print('i=%d:' % (i),a1,a2,w1,w2)
		i+=1
		cmd = 'python3 weight-edges.py -n %s -c -e %s -o %s --a1 %.3f --a2 %.3f --w1 %.3f --w2 %.3f > out.log' % \
		(network,weighted_network,out_dir+'run',a1,a2,w1,w2)
		#print(cmd)
		val = os.system(cmd)
		if val != 0:
			sys.exit()
		weight_file = out_dir+'run_w1_%.3f_w2_%.3f.txt' % (w1,w2)
		these_weights = np.loadtxt(weight_file,skiprows=0,usecols=(2,))
		q75, q25 = np.percentile(these_weights, [75 ,25])
		iqr = q75 - q25
		out.write('%.3f\t%.3f\t%.3f\t%.3f\t%e\n' % (a1,a2,w1,w2,iqr))
		print('    q75=%f, q25=%f, IQR=%f' % (q75,q25,iqr))
		if iqr > best_val:
			best_val = iqr
			best_combo = (a1,a2,w1,w2)
			print('    BEST SO FAR')
	out.close()
	print("BEST IQR IS "+best_val)
	print('BEST COMBO IS %.3f\t%.3f\t%.3f\t%.3f' % (a1,a2,w1,w2))
	print('wrote to '+outprefix+'iqr.txt')

	return


if __name__ == '__main__':
	if len(sys.argv) != 4 and len(sys.argv) != 5:
		sys.exit('USAGE: python3 parameter_sweep.py <INTERACTOME> <WEIGHTED_INTERACTOME> <OUTPREFIX> <FORCE-optional>')
	network = sys.argv[1]
	weighted_network = sys.argv[2]
	outprefix = sys.argv[3]
	if len(sys.argv)==5:
		force=True
	else:
		force=False
	main(network,weighted_network,outprefix,force)