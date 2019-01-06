# import the GraphSpace and GSGraph classes from graphspace_python
# learn more: https://python.org/pypi/graphspace_python
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
import sys
sys.path.append('/Users/aritz/git/fly-interactome/interactome/weighted-interactome/')
import networkx as nx			   
from graph import Graph
import re

NUM_TO_PLOT = 100
GOterms = {'GO:0003383':'apical constriction',
	'GO:0008360':'regulation of cell shape',
	'GO:0030048':'actin filament based movement',
	'GO:0032956':'regulation of actin cytoskeleton organization',
	'GO:0017022':'myosin binding',
	'GO:0016459':'myosin complex',
	'GO:0003779': 'actin binding'}
GOfile = '/Users/aritz/git/fly-interactome/interactome/weighted-interactome/gene_association.gmt'
MIfile = '/Users/aritz/git/fly-interactome/interactome/mi.owl'

def main(network,weighted_network,weighted_nodes,q_val):
	
	print('\nReading the network from %s' %(network))
	net = Graph()
	net.read(network,True) # collapsed=True
	print('\nResetting the network weights')
	with open(weighted_network) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			net[row[0]][row[1]]['weight'] = float(row[2])
	print('\nWeighting nodes')	
	rank = 1
	with open(weighted_nodes) as fin:
		for line in fin:
			if line[0] == '#':
				continue
			row = line.strip().split()
			assert net.has_node(row[0])
			net.add_node(row[0],rank=rank,weight=float(row[1]))

			rank+=1
	print(nx.info(net))

	print('\nReading gene associations from %s' % (GOfile))
	GOgenes = {}
	with open(GOfile) as fin:
		for line in fin:
			row = line.strip().split()
			if row[0] in GOterms:
				GOgenes[row[0]] = row[2:]
	for term in GOterms:
		if term in GOgenes:
			print('%s (%s): %d genes' % (term,GOterms[term],len(GOgenes[term])))
		else: 
			print('%s (%s): not found' % (term,GOterms[term]))

	print('\nReading Molecular Interactio Names from %s' % (MIfile))
	term2name, name2term = read_owl(MIfile)

	print('\nGet Induced Subgraph for plotting')
	node_subset = [n for n,attrs in net.nodes(data=True) if attrs['rank'] <= NUM_TO_PLOT]
	print('  plotting %d nodes' % (len(node_subset)))
	subset_net = net.subgraph(node_subset)
	nx.info(subset_net)


	print('\nMake GraphSpace graph')
	G = GSGraph()
	large_size = 60
	small_size=40
	thick_width=3
	thin_width=1
	for n in subset_net.nodes():
		this_rank = net.node[n]['rank']
		step_val = float(this_rank)/(NUM_TO_PLOT*2)
		this_color = rgb_to_hex(step_val*2,.8,1-(step_val*2))
		popup,inGO=get_node_popup(n,net,GOgenes)
		if inGO:
			this_width=thick_width
		else:
			this_width=thin_width
		print(n,this_rank,this_color,net.degree(n),inGO)
		G.add_node(n,label='%s\n#%d' % (n,this_rank),popup=popup,k=this_rank)
		G.add_node_style(n,shape='ellipse',color=this_color,height=large_size,width=large_size,border_width=this_width)
	for n1,n2 in subset_net.edges():
		this_rank = max(net.node[n1]['rank'],net.node[n2]['rank'])
		G.add_edge(n1,n2,directed=False,popup=get_edge_popup(n1,n2,net,term2name),k=this_rank)
		G.add_edge_style(n1,n2,directed=False,width=2*net.edge[n1][n2]['weight'])
	G.set_name('Top %d Ranked Nodes (teleprob %s)' % (len(node_subset),q_val))
	desc = 'Hey Derek!</br>Nodes are ranked by <b>random walk</b> starting from Node #1 (teleporting back to that node with probability %s.<br>' % (q_val)
	desc += '<b>Bolded Nodes</b> Appear in One of the GO Terms:<br>'
	for term in GOterms:
		desc+= '<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank"><font color="blue">%s</font></a> (%s)<br>' % (term,term,GOterms[term])
	desc+='<br>'
	G.set_data({'description':desc})
	username = 'aritz@reed.edu'
	password = 'platypus'
	graphid = post(G,username,password)
	print('DONE WITH GRAPH',graphid)
	return

def get_node_popup(node,net,GOgenes):
	node_id = net.node[node]['id']
	node_weight = net.node[node]['weight']
	node_rank = net.node[node]['rank']
	popup = '<b>%s</b> (<a href="http://flybase.org/reports/%s.html" target="_blank"><font color="blue">%s</font></a>)<br>Rank=%d<br>weight=%f<br>%d Neighbors<br>' % (node,node_id,node_id,node_rank,node_weight,net.degree(node))
	inGO = False
	if any([node in GOgenes[term] for term in GOgenes]):
		inGO = True
		popup += '<br><b>GO Terms</b>(<a href="http://amigo.geneontology.org/amigo/gene_product/FB:%s" target="_blank"><font color="blue">AmiGO link</font></a>):<br>' % (node_id)
		for term in GOgenes:
			if node in GOgenes[term]:
				popup += '<a href="http://amigo.geneontology.org/amigo/term/%s" target="_blank"><font color="blue">%s</font></a>: %s<br>' % (term,term,GOterms[term])
	return popup,inGO

def get_edge_popup(node1,node2,net,term2name):
	node1_id = net.node[node1]['id']
	node2_id = net.node[node2]['id']
	edge_weight = net.edge[node1][node2]['weight']
	pmids = net.edge[node1][node2]['pmids']
	dbs = net.edge[node1][node2]['dbs']
	ev_types = net.edge[node1][node2]['types']
	popup = '%s-%s<br>%s-%s<br>Weight=%f<br>' % (node1,node2,node1_id,node2_id,edge_weight)
	popup += '<br><b>%d PubMed IDs</b>:<br>' % (len(pmids))
	for p in pmids:
		if p == 'unpublished':
			popup += 'unpublished<br>'
		else:
			popup += '<a href="https://www.ncbi.nlm.nih.gov/pubmed/%s" target="_blank"><font color="blue">%s</font></a><br>' % (p,p)
	popup += '<br><b>%d Evidence Types</b>:<br>' % (len(ev_types))
	for t in ev_types:
		if 'MI:' in t:
			popup += '<a href="http://purl.obolibrary.org/obo/%s" target="_blank"><font color="blue">%s</font></a>: %s<br>' % (t.replace(':','_'),t,term2name[t])
		else:
			popup += '%s<br>' % (t)
	popup += '<br><b>%d Databases:</b>:<br>' % len(dbs)
	for db in dbs:
		popup += '%s<br>' % (db)
	return popup


def rgb_to_hex(r,g,b,to255=True):
	return '#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255))

def read_owl(infile):
	filestring = open(infile).read()
	filelist = filestring.split('[Term]')
	filelist = filelist[1:] # skip first one

	pattern = re.compile('id: (MI:\d+).*\nname: (.*?)\n')
	term2name = {}
	name2term = {}
	for item in filelist:
		m = re.search(pattern,item)
		term = m.group(1)
		name = m.group(2).replace(' ','-')
		term2name[term] = name
		name2term[name] = term
		

	return term2name, name2term


# Function to post a graph to graphspace with a specific
# username and password.  
# Inputs: Graph (GSGraph Object), username (string), password (string)
# Output: the graph ID (int)
def post(G,username,password):
  # connect to GraphSpace with the username and password.
  gs = GraphSpace(username,password)
  try:
    # try updating the graph. If the graph does not exist, 
    # this will throw an error.
    graph = gs.update_graph(G)
  except:
    # catch the error and try posting a new graph.
    graph = gs.post_graph(G)
  return graph.id

if __name__ == '__main__':
	if len(sys.argv) != 5:
		sys.exit('USAGE: python viz-network.py <network> <weighted_network> <weighted_nodes> <q_val>')
	network = sys.argv[1]
	weighted_network = sys.argv[2]
	weighted_nodes = sys.argv[3]
	q_val = sys.argv[4]
	main(network,weighted_network,weighted_nodes,q_val)
