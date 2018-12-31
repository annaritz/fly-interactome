'''
Annotations and GOdag classes
'''
__author__ = '''Christopher L. Poirel (chris.poirel@gmail.com)'''

import sys
import networkx as nx

class Annotations():
    
    def __init__(self):
        self.functions = set()
        self.descriptions = {}
        self.genesets = {}
        self.annotatedGenes = set()
    
    def readGMT(self, infile, mapper={}):
        '''
        Read annotations from infile, which should be in the GMT format.
        The GMT format is a tab delimited file where each row represents a
        gene set. The first item in the row is the name of the gene set, the
        second item is a string description of the gene set (or "na"), and
        additional items are the tab-delimited list of genes in that gene set.
        '''
        f = open(infile, 'r')
        for line in f:
            if line[0]=='#': #ignore comments
                continue
            
            items = line.rstrip().split('\t')
            if len(items)<3:
                continue
            
            func = items[0]
            descr = items[1]
            geneset = set([mapper.get(i, i) for i in items[2:]])
            
            self.functions.add(func)
            self.descriptions[func] = descr
            self.genesets[func] = geneset
            self.annotatedGenes.update(geneset)
        f.close()

    def readGMT_FuncAssociate(self, infile, mapper={}):
        '''
        Read annotations from a downloaded association file from FuncAssociate.
        This is the same as the GMT format except that the genes in the gene set 
        are delimited by a space instead of a tab.  See readGMT().
        '''
        f = open(infile, 'r')
        for line in f:
            if line[0]=='#': #ignore comments
                continue
            
            items = line.rstrip().split('\t')
            if len(items)<3:
                continue
            
            func = items[0]
            descr = items[1]
            ## this is the only line that changed from readGMT()
            geneset = set(items[2].split(' '))
            
            self.functions.add(func)
            self.descriptions[func] = descr
            self.genesets[func] = geneset
            self.annotatedGenes.update(geneset)
        f.close()

    def readTSV(self, infile, mapper={}):
        '''
        Read annotations from infile, which should be a tab delimited file where
        each row represents a single gene-function pair. There should be two columns
        with column headers "function" and "gene". The header line can optionally
        begin with a "#" symbol. All other lines beginning with "#" will be ignored.
        '''
        f = open(infile, 'r')
        headeritems = f.readline().lstrip('#').rstrip().split('\t')
        
        if ('gene' not in headeritems) or ('function' not in headeritems):
            print 'Incorrect header in the annotations file. Required: "function", "gene"'
            return
        
        functionCol = headeritems.index('function')
        geneCol = headeritems.index('gene')
        
        for line in f:
            if line[0]=='#':
                continue
            
            items = line.rstrip().split('\t')
            if len(items)<max(functionCol, geneCol)+1:
                continue
            
            gene = mapper.get(items[geneCol], items[geneCol])
            func = items[functionCol]
            
            if func in self.functions:
                self.genesets[func].add(gene)
                self.annotatedGenes.add(gene)
            else:
                self.functions.add(func)
                self.descriptions[func] = 'NA'
                self.genesets[func] = set([gene])
                self.annotatedGenes.add(gene)
    
        f.close()
        
    def getFunctions(self):
        return self.functions
        
    def getAnnotatedGenes(self, func=None):
        if func==None:
            return self.annotatedGenes
        else:
            return self.genesets.get(func, None)
        
    def getDescription(self, func):
        return self.descriptions.get(func, None)
    
    def printAnnotations(self, outfile=None):
        ostream = sys.stdout
        if outfile!=None:
            ostream = open(outfile, 'w')
        for f in self.functions:
            ostream.write('%s\t%s\t%s\n' %(f, self.descriptions[f], '\t'.join(self.genesets[f])))

        if ostream!=sys.stdout:
            ostream.close()

    def keepAnnotationsForGenes(self, genes):
        '''Keep only annotations for the given genes'''
        
        self.annotatedGenes.intersection_update(genes)
        
        emptyFuncs = []
        for func, genesets in self.genesets.iteritems():
            genesets.intersection_update(genes)
            if len(genesets)==0:
                emptyFuncs.append(func)
                
        for func in emptyFuncs:
            self.descriptions.pop(func)
            self.functions.remove(func)
            self.genesets.pop(func)

# some GO specific functions below 
    def applyTruePathRule(self, godag, typesToApply=set(['is_a', 'part_of'])):
        '''Apply the true path rule to the current annotations given the GO DAG.'''
        
        print 'Applying the true path rule up the GO DAG.'
        topo = nx.topological_sort(godag)
        
        for func in topo:
            # nothing to do if this function annotates 0 genes
            if func not in self.functions:
                continue
            for neighbor in godag.neighbors(func):
                # only apply true-path rule up the specified link types
                if len(godag.edge[func][neighbor]['link_types'].intersection(typesToApply))==0:
                    continue
                if neighbor not in self.functions:
                    self.functions.add(neighbor)
                    self.descriptions[neighbor] = 'NA'
                    self.genesets[neighbor] = set()
                self.genesets[neighbor].update(self.genesets[func])

    def addFunctionDescriptions(self, godag):
        for f in self.descriptions.keys():
            if f in godag:
                self.descriptions[f] = godag.names[f]
                
    def keepAnnotationsForNamespace(self, godag, namespace):
        '''
        namespace can be: biological_process, b, molecular_function, m, cellular_component, c
        '''
        validNamespaces = {'biological_process' : 'biological_process', \
                           'molecular_function' : 'molecular_function', \
                           'cellular_component' : 'cellular_component', \
                           }
        if namespace not in validNamespaces.keys():
            print 'WARNING: illegal namespace given.'
            return
        namespace = validNamespaces[namespace]
        
        for func in self.functions.copy():
            if godag.getNamespace(func)==namespace:
                continue
            self.descriptions.pop(func)
            self.functions.remove(func)
            self.genesets.pop(func)
        
        # recompute self.annotatedGenes, since we potentially removed several functions.
        self.annotatedGenes.clear()
        for s in self.genesets.values():
            self.annotatedGenes.update(s)     
                
class GOdag(nx.DiGraph):
    
    def __init__(self, obofile=None):
        nx.DiGraph.__init__(self)
        self.names = {}
        self.namespaces = {}
        self.altIds = {}
        self.synonyms = {}

        # for obsolete terms
        self.replacedBy = {}
        self.consider = {}
        
        if obofile!=None:
            print 'Reading GO DAG from "%s".' %(obofile)
            self.readOBO(obofile)
        
    def readOBO(self, obofile):
        linkTypes = {'is_a'                 : 0, \
                     'part_of'              : 0, \
                     'regulates'            : 0, \
                     'positively_regulates' : 0, \
                     'negatively_regulates' : 0, \
                     'develops_from'        : 0, \
                     'related_to'           : 0, \
                    }
    
        infile = open(obofile, 'r')
        
        # items we should identify for each term we read
        goid = None
        name = None
        namespace = None
        
        isObsolete = False
        curr_replaced = []
        curr_consider = []
        curr_alt = []
        curr_syn = []
        
        readingTerm = False
        linenum=0
        for line in infile:
            currstring = line.rstrip()
            linenum += 1
            
            # we must handle blank lines first
            if currstring=='' and readingTerm:
                readingTerm = False
                # now that I'm finished, update some information about the GO term I just read
                if isObsolete:
                    if len(curr_replaced)>0:
                        self.replacedBy[goid] = list(set(curr_replaced))
                    if len(curr_consider)>0:
                        self.consider[goid] = list(set(curr_consider))
                else:
                    self.names[goid] = name
                    self.namespaces[goid] = namespace
                    self.synonyms[goid] = curr_syn
                    for a in curr_alt:
                        self.altIds[a] = goid
                
                # and set those values back to None
                goid = None
                name = None
                namespace = None
                isObsolete = False
                curr_replaced = []
                curr_consider= []
                curr_alt = []
                curr_syn = []
                continue
            elif currstring=='' and not readingTerm:
                continue
            
            # ignore comments
            if currstring[0]=='!':
                continue
                    
            # we may be starting a new term
            if currstring=='[Term]':
                readingTerm = True
                continue
            
            # if we aren't reading a term, ignore the line
            if not readingTerm:
                continue
            
            # now, we need to handle what can occur inside a term
            tag, sep, value = currstring.partition(': ')
            if tag == 'is_obsolete':
                isObsolete = True
                continue 
            if tag=='id':
                goid = value
                continue
            elif tag=='name':
                name = value
                continue
            elif tag=='namespace':
                namespace = value
                continue
            elif tag=='alt_id':
                curr_alt.append(value.strip())
            elif tag in linkTypes:
                linkTypes[tag] += 1
                target = value.split()[0]
                self.add_edge(goid, target)
                self.add_link_type(goid, target, tag)
                continue
            elif tag=='relationship':
                valueItems = value.split()
                if len(valueItems)<2:
                    continue
                ltype = valueItems[0]
                target = valueItems[1]
                if ltype in linkTypes:
                    linkTypes[ltype] += 1
                    self.add_edge(goid, target)
                    self.add_link_type(goid, target, ltype)
                continue
            elif tag=='replaced_by':
                curr_replaced.append(value.strip())
            elif tag=='consider':
                curr_consider.append(value.strip())
            elif tag=='synonym':
                curr_syn.append(value.split('"')[1]) # take text in quotes
            
                
        infile.close()
        
        for ltype, count in linkTypes.iteritems():
            print '\tFound %d "%s" links.' %(count, ltype)

    def add_link_type(self, u, v, ltype):
        if 'link_types' in self.edge[u][v]:
            self.edge[u][v]['link_types'].add(ltype)
        else:
            self.edge[u][v]['link_types'] = set([ltype])

    def getDescendants(self, term, linkTypes=set(['is_a', 'part_of'])):
        '''Return all descendants of the given term'''
        children = set()
        tovisit = set([term])
        while len(tovisit)!=0:
            curr = tovisit.pop()
            children.add(curr)
            tovisit.update( set([p for p in self.predecessors(curr) if len(self.edge[p][curr]['link_types'].intersection(linkTypes))>0]) )
            tovisit.difference_update(children)
        return children
    
    def getAncestors(self, term, linkTypes=set(['is_a', 'part_of'])):
        '''Return all ancestors of the given term'''
        anc = set()
        tovisit = set([term])
        while len(tovisit)!=0:
            curr = tovisit.pop()
            anc.add(curr)
            tovisit.update( set([p for p in self.successors(curr) if len(self.edge[curr][p]['link_types'].intersection(linkTypes))>0]) )
            tovisit.difference_update(anc)
        return anc
    
    def commonAncestors(self, term1, term2, linkTypes=set(['is_a', 'part_of'])):
        '''Return all common ancestors shared by term1 and term2.'''
        anc1 = self.getAncestors(term1, linkTypes)
        anc2 = self.getAncestors(term2, linkTypes)
        return anc1.intersection(anc2) 
    
    def getNamespace(self, term):
        '''Return the namespace of the given term.'''
        return self.namespaces.get(term, '')
