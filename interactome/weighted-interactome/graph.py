# uncompyle6 version 3.2.5
# Python bytecode 2.7 (62211)
# Decompiled from: Python 3.5.1 (default, Dec 24 2015, 12:10:56) 
# [GCC 4.2.1 Compatible Apple LLVM 7.0.2 (clang-700.1.81)]
# Embedded file name: /Users/aritz/git/fly-interactome/interactome/weighted-interactome/graph.py
# Compiled at: 2018-12-31 09:01:11
"""
Subclasses NetworkX Graph
"""
__author__ = 'Anna Ritz (modified from Christopher L. Poirel (chris.poirel@gmail.com))'
import sys, networkx as nx
from collections import deque

class Graph(nx.Graph):
    """
    Subclass of the networkx DGraph class
    """

    def __init__(self, name=''):
        nx.Graph.__init__(self, name=name)

    def read(self, filename, collapsed=True):
        infile = open(filename, 'r')
        eweight = 1
        for line in infile:
            items = [ x.strip() for x in line.rstrip().split('\t') ]
            if line == '' or line[0] == '#' or len(items) < 2:
                continue
            name1 = items[0]
            name2 = items[1]
            if name1 == name2:
                continue
            if collapsed:
                if len(items) != 7:
                    print (
                     'SKIPPING LINE', line)
                    continue
                pubmedids = items[2].split(';')
                id1 = items[3]
                id2 = items[4]
                databases = items[5].split(';')
                evidence = set()
                for item in items[6].split(';'):
                    e = item.split(':')
                    if len(e) == 1:
                        continue
                    if len(e) == 2:
                        evidence.add(e[1])
                    else:
                        evidence.add((':').join(e[1:]))

                evidence = list(evidence)
            else:
                sys.exit('IMPLEMENT!!')
            if not self.has_node(name1):
                self.add_node(name1, id=id1)
            if not self.has_node(name2):
                self.add_node(name2, id=id2)
            self.add_edge(name1, name2, weight=eweight, types=evidence, dbs=databases, pmids=pubmedids)

    def add_edge(self, u, v, attr_dict=None, **attr):
        nx.Graph.add_edge(self, u, v, attr_dict=attr_dict, **attr)
        if 'weight' not in self[u][v]:
            self[u][v]['weight'] = 1
        for attr in ['types', 'dbs', 'pmids']:
            if attr not in self[u][v]:
                self[u][v][attr] = []

    def add_edges_from(self, ebunch, attr_dict=None, **attr):
        for e in ebunch:
            u, v = e[0:2]
            if len(e) > 2:
                attr_dict = e[2]
            self.add_edge(u, v, attr_dict=attr_dict, **attr)

    def printGraph(self, outfile=None, weight='weight'):
        ostream = sys.stdout
        if outfile != None:
            ostream = open(outfile, 'w')
        ostream.write('#tail\thead\tweight\ttype\n')
        for tail, head, edata in self.edges_iter(data=True):
            if len(edata['types']) == 0:
                ostream.write('%s\t%s\t%0.5e\t%s\n' % (tail, head, edata.get(weight, 1.0), ''))
                continue
            for t in edata['types']:
                ostream.write('%s\t%s\t%0.5e\t%s\n' % (tail, head, edata.get(weight, 1.0), t))

        if ostream != sys.stdout:
            ostream.close()
        return

    def getEdgeTypes(self):
        etypesDict = {}
        for t, h, data in self.edges_iter(data=True):
            for etype in data['types']:
                if etype not in etypesDict:
                    etypesDict[etype] = set()
                etypesDict[etype].add((t, h))

        return etypesDict

    def reachable(self, nodes):
        visited = set()
        toVisit = deque(nodes)
        reachable = set()
        while len(toVisit) > 0:
            t = toVisit.pop()
            visited.add(t)
            reachable.add(t)
            toVisit.extend(set([ h for h in self.successors(t) if h not in visited and h not in toVisit ]))

        return reachable

    def count_uni_bi_directed_edges(self):
        if self.number_of_selfloops() > 0:
            print 'There are %d selfloops in the graph. Excluding them from counts below' % self.number_of_selfloops()
        else:
            print 'There are no selfloops in the graph.'
        G = self.copy()
        G.remove_edges_from(G.selfloop_edges())
        counter_unidirected = 0
        counter_bidirected = 0
        for e in G.edges_iter():
            if G.has_edge(e[1], e[0]):
                counter_bidirected += 1
            else:
                counter_unidirected += 1

        return 'Number of uni-directed edges: %d\nNumber of bi-directed edges: %d' % (counter_unidirected, counter_bidirected / 2)
# okay decompiling graph.pyc
