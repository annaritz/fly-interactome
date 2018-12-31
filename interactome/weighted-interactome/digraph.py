'''
Subclasses NetworkX DiGraph
'''
__author__ = '''Christopher L. Poirel (chris.poirel@gmail.com)'''

import sys
import networkx as nx
from collections import deque


class DiGraph(nx.DiGraph):
    '''
    Subclass of the networkx DiGraph class
    '''
    def __init__(self, name=''):
        nx.DiGraph.__init__(self, name=name)

    def read(self, filename, typesToIgnore=[]):
        infile = open(filename, 'r')
        for line in infile:
            items = [x.strip() for x in line.rstrip().split('\t')]
            if line=='':
                continue
            if line[0]=='#':
                continue
            if len(items) < 2:
                continue

            id1 = items[0]
            id2 = items[1]
            eweight = 1
            etypes = []
            
            if id1==id2:
                continue
            
            if len(items)>2:
                eweight = float(items[2])
            if len(items)>3:
                newType = items[3]
                if newType in typesToIgnore:
                    continue
                # if the edge exists, record known types
                if self.has_edge(id1, id2):
                    edata = self.get_edge_data(id1, id2)
                    if 'types' in edata:
                        etypes = edata['types']
                    self.remove_edge(id1, id2)
                etypes.append(newType)

            etypes = list(set(etypes)) # remove duplicates
            self.add_edge(id1, id2, weight=eweight, types=etypes)

    def add_edge(self, u, v, attr_dict=None, **attr):
        nx.DiGraph.add_edge(self,u,v,attr_dict=attr_dict,**attr)
        if 'weight' not in self[u][v]:
            self[u][v]['weight'] = 1
        if 'types' not in self[u][v]:
            self[u][v]['types'] = []

    def add_edges_from(self, ebunch, attr_dict=None, **attr):  
        for e in ebunch:
            u,v=e[0:2]
            if len(e)>2:
                attr_dict = e[2]
            self.add_edge(u,v,attr_dict=attr_dict,**attr)

    def printGraph(self, outfile=None, weight='weight'):
        ostream = sys.stdout
        if outfile!=None:
            ostream = open(outfile, 'w')
            
        ostream.write('#tail\thead\tweight\ttype\n')
        for tail, head, edata in self.edges_iter(data=True):
            if len(edata['types'])==0:
                ostream.write('%s\t%s\t%0.5e\t%s\n' %(tail, head, edata.get(weight, 1.0), ''))
                continue
            for t in edata['types']:
                ostream.write('%s\t%s\t%0.5e\t%s\n' %(tail, head, edata.get(weight, 1.0), t))
        if ostream!=sys.stdout:
            ostream.close()
    def getEdgeTypes(self):
        etypesDict = {}
        for (t,h,data) in self.edges_iter(data=True):
            for etype in data['types']:
                if etype not in etypesDict:
                    etypesDict[etype] = set()
                etypesDict[etype].add((t,h))
        return etypesDict
    
    def reachable(self, nodes):
        visited = set()
        toVisit = deque(nodes)
        reachable = set()
        while len(toVisit)>0:
            t = toVisit.pop()
            visited.add(t)
            reachable.add(t)
            toVisit.extend( set([h for h in self.successors(t) if h not in visited and h not in toVisit]) )
        return reachable

    def count_uni_bi_directed_edges(self):
        if self.number_of_selfloops() > 0:
            print "There are %d selfloops in the graph. Excluding them from counts below" %self.number_of_selfloops()
        else:
            print "There are no selfloops in the graph."
        G = self.copy()
        G.remove_edges_from(G.selfloop_edges())
        counter_unidirected = 0
        counter_bidirected = 0
        for e in G.edges_iter():
            if G.has_edge(e[1], e[0]):
                counter_bidirected += 1
            else:
                counter_unidirected += 1
        return "Number of uni-directed edges: %d\nNumber of bi-directed edges: %d" %(counter_unidirected , (counter_bidirected/2))
        
        
