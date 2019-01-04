from __future__ import print_function
#!/usr/bin/python

## Bayesian weighting Implementation for Fly Interactome

## History:
## Written by Chris Poirel 
## Updated by Anna Ritz Septembre 2014
##
## Implements ResponseNet's interactome weighting scheme.
## Detailed in Chris's paper (supplementary section):
## Top-down network analysis to drive bottom-up modeling of physiological processes.
## Poirel CL, Rodrigues RR, Chen KC, Tyson JJ, Murali TM.
## J Comput Biol. 2013 May;20(5):409-18. doi: 10.1089/cmb.2012.0274.
## http://www.ncbi.nlm.nih.gov/pubmed/23641868

## General python imports
from optparse import OptionParser   # for parsing options
import sys                          # for sys.exit() calls when there's an error
import random                       # random number generator
from math import log

## NetworkX imports
import networkx as nx               
from graph import Graph

## From SVN: src/python/scripts/trunk/
from utilsPoirel import *

## From SVN: src/python/CellCycle (same location as this file)
from annotations import Annotations
from annotations import GOdag

################################################
def main(args):

    # parse arguments
    opts = parse_arguments(args)

    # get the undirected network
    print('\nReading the network from %s' %(opts.network))
    net = Graph()
    net.read(opts.network)
    print(nx.info(net))
    # get edge types from the network
    etypes = net.getEdgeTypes()
    
    if opts.probs:
        print('Reading all required information from file "%s"' % (opts.probs))
        etypeProbs,num_pos,num_neg = read_etype_probs(opts,net,etypes)
    else:
        # merge edges from evidences that annotate <25 edges.
        print('Get miscellaneous edge types.')
        output = open('%s-miscellaneous_edge_types.txt' %(opts.outprefix), 'w')
        output.write('#EdgeType\tNumEdges\n')
        to_merge = [e for e in etypes if len(etypes[e]) < 25]
        etypes['miscellaneous'] = set()
        for e in to_merge:
            output.write('%s\t%d\n' % (e,len(etypes[e])))
            etypes['miscellaneous'].update(etypes[e])
            del etypes[e] # delete the small type
        output.close()

        print('Collecting GO annotations, positive functions, and sampling pos and neg sets.')
        # read mapping file.
        if opts.mapper != None:
            mapper = readDict(opts.mapper,opts.fromcol,opts.tocol)
        else:
            mapper = {}

        # get GO annotations
        ann = get_annotations(opts,net,mapper)
       
        # get and filter positive functions
        posFuncs = get_pos_funcs(opts,ann)

        # compute the weight of individual experiments
        etypeProbs,num_pos,num_neg = compute_etype_probs(opts,net,etypes,ann,posFuncs)
        
        # compute the weight of edge with only experiment k
        compute_single_evidence_edge_weights(opts,etypes,etypeProbs,num_pos,num_neg)
    
    # Finally, we can weight the edges.
    # Note, priors are calculated previously, but they are so simmple
    # recalculate them again here.
    pr_I1 = float(num_pos)/(num_pos + num_neg) # Pr(I=1)
    pr_I0 = 1.0 - pr_I1                        # Pr(I=0)
    weight_edges(opts,net,etypes,etypeProbs,pr_I1,pr_I0)

    print('\nWrote to positive GO terms to %s-positive-GO-terms.txt' %(opts.outprefix))
    print('Wrote edge types that are merged to "miscellaneous" to %s-miscellaneous_edge_types.txt' %(opts.outprefix) )
    print('Wrote edge type probabilities to %s-edge_type_probs.txt' %(opts.outprefix))
    print('Wrote edge type weights to %s-edge_type_weights.txt' %(opts.outprefix))
    print('Wrote weighted interactome to %s.txt' %(opts.outprefix))
    print('DONE.')
    return

################################################
def get_annotations(opts,net,mapper):
    # get the annotations
    print('\nReading the annotations from %s' %(opts.annotations))
    ann = Annotations()
    ann.readGMT(opts.annotations, mapper)
    print('\t%d functions' %(len(ann.getFunctions())))
    print('\t%d annotated genes' %(len(ann.getAnnotatedGenes())))
    
    # remove annotations for genes that are not in the network    
    print('\nKeeping annotations for only genes in the network')
    ann.keepAnnotationsForGenes(net.nodes())
    print('\t%d functions' %(len(ann.getFunctions())))
    print('\t%d annotated genes' %(len(ann.getAnnotatedGenes())))

    return ann

################################################
def get_pos_funcs(opts,ann):
    # get the GO DAG
    print('')
    godag = GOdag(opts.ontology)
    ann.applyTruePathRule(godag)
    ann.keepAnnotationsForNamespace(godag, 'biological_process')
    print('\nKeeping only biological processes.')
    print('\t%d functions' %(len(ann.getFunctions())))
    print('\t%d annotated genes' %(len(ann.getAnnotatedGenes())))
    
    # read functions
    posFuncs = readItemSet(opts.functions,1)
    posSubgraph = godag.subgraph(posFuncs)
    posFuncs = set(posSubgraph.nodes())

    print('\nRemoving functions: (1) not in GOdag, (2) too large, (3) descendant of another function, and (4) too small, in that order.')

    ## remove any functions that aren't in the GOdag
    posFuncs = set([f for f in posFuncs if f in ann.genesets])
    print('%d GO terms after removing terms not in GOdag' % (len(posFuncs)))
    if len(posFuncs)<50:
        for p in posFuncs:
            print('  %s: %d genes' % (p,len(ann.genesets[p])))

    ## remove functions that are too large
    posFuncs = set([f for f in posFuncs if len(ann.genesets[f])<=opts.maxsetsize])
    print('%d GO terms after removing terms with >%d genes' % \
        (len(posFuncs),opts.maxsetsize))

    ## remove any functions with an ancestor in the list
    # for f in posFuncs:
    #     print f,ann.descriptions[f],len(ann.genesets[f])
    #     A = godag.getAncestors(f)
    #     for a in A.intersection(posFuncs):
    #         print ' ancestor',a,ann.descriptions[a],len(ann.genesets[a])      
    #     sys.exit()
    posFuncs = set([f for f in posFuncs if len(godag.getAncestors(f).intersection(posFuncs))==1])
    print('%d GO terms after removing descendants' % (len(posFuncs)))

    ## remove functions that are too small
    posFuncs = set([f for f in posFuncs if len(ann.genesets[f])>=opts.minsetsize])
    print('%d GO terms after removing terms with <%d genes' % \
        (len(posFuncs),opts.minsetsize))

    if len(posFuncs) == 0: # if we've removed all functions, exit.
        sys.exit('ERROR: removed all GO terms.')

    # if there are fewer than 50, print them.
    if len(posFuncs)<50:
        for p in posFuncs:
            print('  %s: %d genes' % (p,len(ann.genesets[p])))

    ## write selected terms to file.
    output = open('%s-positive-GO-terms.txt' %(opts.outprefix), 'w')
    output.write('#numgenes\tterm\tdescription\n')
    for term in posFuncs:
        output.write('%d\t%s\t%s\n' % (len(ann.genesets[term]),term,ann.descriptions[term]))
    output.close()

    return posFuncs

################################################
def get_pos_neg_edges(opts,net,ann,posFuncs):

    # major difference from Poirel et al weighting; all pairs are all INTERACTING
    # pairs (that is, all nodes). It used to be all n choose 2 possibilities. This
    # is commented out to illustrate difference.
    allPairs = set(net.edges())
    #allPairs = set( [(u,v) for u in net.nodes_iter() for v in net.nodes_iter() if u<v] )
    

    pos = set()
    neg = set()
    print('\nGenerating positive and negative pairs.' )
    print('Using the union of %d functions to generate positives' %(len(posFuncs)))
    print('\t%d total pairs (number of edges)' % (len(allPairs)))
    # get all pairs
    for f in posFuncs:
        funcgenes = ann.genesets[f]
        pos.update([(u,v) for u in funcgenes for v in funcgenes])
    # retain those that are edges
    pos = pos.intersection(allPairs)
    print('\t%d positive pairs' %(len(pos)))
    
    # randomly sample non-positive pairs for negative examples
    if opts.samplesize*len(pos) < len(allPairs):
        neg = set( random.sample(allPairs.difference(pos), opts.samplesize*len(pos)) )
    else:
        print('\tWARNING: not sampling: proportion of negatives is %.4f times the size of positives' % \
            ((len(allPairs)-len(pos))/float(len(pos))))
        neg = allPairs.difference(pos)
    print('\t%d negative pairs' %(len(neg)))
    return pos,neg

################################################
def compute_etype_probs(opts,net,etypes,ann,posFuncs):

    # get positive and negative edges
    pos,neg = get_pos_neg_edges(opts,net,ann,posFuncs)
    num_pos = len(pos)
    num_neg = len(neg)

    print('\nComputing the weight of individual experiments Pr(I=1|E=1).')
    output = open('%s-edge_type_probs.txt' %(opts.outprefix), 'w')
    output.write('# %d positives\n' % (len(pos)))
    output.write('# %d negatives\n' % (len(neg)))
    output.write('#edge_type\tnumedges\tnumpos\tnumneg\n')
    etypeProbs = {}
    for et, edges in etypes.iteritems():
        
        pr_E1_I1 = float(len(pos.intersection(edges))) / len(pos) # Pr(E=1|I=1)
        pr_E0_I1 = 1.0 - pr_E1_I1                                 # Pr(E=0|I=1)
        
        pr_E1_I0 = float(len(neg.intersection(edges))) / len(neg) # Pr(E=1|I=0)  
        pr_E0_I0 = 1.0 - pr_E1_I0                                 # Pr(E=0|I=0)
        etypeProbs[et] = {'pr_E1_I1' : pr_E1_I1, \
                          'pr_E1_I0' : pr_E1_I0, \
                          'pr_E0_I1' : pr_E0_I1, \
                          'pr_E0_I0' : pr_E0_I0, \
                         }        
        p = float(len(pos.intersection(edges)))
        n = float(len(neg.intersection(edges)))
        output.write('%s\t%d\t%d\t%d\n' % (et,len(edges),p,n))
    output.close()

    return etypeProbs,num_pos,num_neg

def read_etype_probs(opts,net,etypes):
    ## like above compute_etype_probs but reads all information from the 
    ## opts.prob file instead.

    fin = open(opts.probs)
    # first two lines are number positives and number negatives.
    num_pos = int(fin.readline().split()[1])
    num_neg = int(fin.readline().split()[1])
    # next line is a header
    fin.readline()

    # for remainder, populated etypeProbs
    etypeProbs = {}
    for line in fin:

        #edge_type  numedges    numpos  numneg
        row = line.strip().split()
        et = row[0]
        numedges = int(row[1])
        p = int(row[2])
        n = int(row[3])
        
        pr_E1_I1 = float(p) / num_pos # Pr(E=1|I=1)
        pr_E0_I1 = 1.0 - pr_E1_I1     # Pr(E=0|I=1)
        
        pr_E1_I0 = float(n) / num_neg # Pr(E=1|I=0)  
        pr_E0_I0 = 1.0 - pr_E1_I0     # Pr(E=0|I=0)
        etypeProbs[et] = {'pr_E1_I1' : pr_E1_I1, \
                          'pr_E1_I0' : pr_E1_I0, \
                          'pr_E0_I1' : pr_E0_I1, \
                          'pr_E0_I0' : pr_E0_I0, \
                         }        
    return etypeProbs,num_pos,num_neg


################################################
def weight_edges(opts,net,etypes,etypeProbs,pr_I1,pr_I0):
    print('\nWeighting each edge with the new edge weight.')
    output = open('%s.txt' %(opts.outprefix), 'w')
    output.write('#tail\thead\tedge_weight\tedge_type\n')
    for t,h in net.edges():
        on = set(net[t][h]['types'])
        
        off = set(etypes.keys())
        off.difference_update(on)
        
        num_I1 = pr_I1
        num_I0 = pr_I0
        
        for e in on:
            num_I1 *= etypeProbs.get(e, etypeProbs['miscellaneous'])['pr_E1_I1']
            num_I0 *= etypeProbs.get(e, etypeProbs['miscellaneous'])['pr_E1_I0']

        for e in off:
            num_I1 *= etypeProbs.get(e, etypeProbs['miscellaneous'])['pr_E0_I1']
            num_I0 *= etypeProbs.get(e, etypeProbs['miscellaneous'])['pr_E0_I0']
        
        if (num_I1+num_I0)==0:
            weight = 0
        else:
            weight = num_I1/(num_I1+num_I0)

        output.write('%s\t%s\t%0.5e\t%s\n' % (t,h,weight,'|'.join(on)))
    output.flush()
    output.close()
    return

################################################
def compute_single_evidence_edge_weights(opts,etypes,etypeProbs,num_pos,num_neg):
    
    pr_I1 = float(num_pos)/(num_pos + num_neg) # Pr(I=1)
    pr_I0 = 1.0 - pr_I1                        # Pr(I=0)

    print('\nComputing the weight of edge with only experiment k Pr(I=1|E), where E_k=1.')
    etypeWeights = {}
    for currType in etypes.keys():
        on = [currType]
        off = etypes.keys()
        off.remove(currType)
        
        num_I1 = pr_I1
        num_I0 = pr_I0
        
        for e in on:
            num_I1 *= etypeProbs[e]['pr_E1_I1']
            num_I0 *= etypeProbs[e]['pr_E1_I0']

        for e in off:
            num_I1 *= etypeProbs[e]['pr_E0_I1']
            num_I0 *= etypeProbs[e]['pr_E0_I0']
            
        etypeWeights[currType] = {}
        etypeWeights[currType]['num_I1'] = float(num_I1)
        etypeWeights[currType]['num_I0'] = float(num_I0)
        #etypeWeights[currType]['ratio'] = float(num_I1/num_I0)
        if (num_I1+num_I0)>0:
            etypeWeights[currType]['weight'] = float(num_I1)/(num_I1+num_I0)
        else:
            etypeWeights[currType]['weight'] = 0.0
    
    output = open('%s-edge_type_weights.txt' %(opts.outprefix), 'w')
    output.write('#total\tPr(I=1|E)\tPr(I=0|E)\tprob\tedge_type\n')
    for e, w in sorted(etypeWeights.iteritems(), key=lambda x: x[1]['weight'], reverse=True):
        output.write('%d\t%.5e\t%.5e\t%.5e\t%s\n' %(len(etypes[e]), w['num_I1'], w['num_I0'], w['weight'], e))
    output.close()
    return

################################################
def parse_arguments(args):
    
    usage = '''
    python weight-edges-by-evidence.py [OPTIONS]
    '''

    parser = OptionParser(usage=usage)
    parser.add_option('-n','--network',type='string',metavar='STR',\
                          help='protein-protein interactome (tab-delimited). First two columns are the tail and head nodes, the last column is the evidence type. Required.')
    parser.add_option('-a','--annotations',type='string',metavar='STR',\
                          help='GO annotation file.  Multi-column tab delimited with "GO term" "description" "tab-delimited list of common names". Required unless --probs specified.')
    parser.add_option('-t','--ontology',type='string',metavar='STR',\
                          help='GO Ontology (OBO) file.  Downloaded from the Gene Ontology database. Required.')
    parser.add_option('-f','--functions',type='string',metavar='STR',\
                          help='File of GO terms that are "positives".  Originally, a 2-column tab-delimited with "GO term" "description". Now, only first column is needed. Required unless --probs specified.')
    parser.add_option('-o','--outprefix',type='string',metavar='STR',\
                          help='output file prefix. Required.')
    parser.add_option('','--probs',type='string',metavar='STR',\
                        help='Bypass calculation of evidence type probabilities - pass it a *edge_type_probs.txt file from a previous run. Optional')
    parser.add_option('-m','--mapper',type='string',metavar='STR',\
                          help='ID Mapping File. Tab-delimited. Optional')
    parser.add_option('','--fromcol',type='int',metavar='INT',\
                          help='Mapping from column x (systematic name), indexed starting from 1. Optional.')
    parser.add_option('','--tocol',type='int',metavar='INT',\
                          help='Mapping to column y (systematic name), indexed starting from 1. Optional.')
    parser.add_option('','--randseed',type='int',metavar='INT',default=1234567,\
                      help='Seed for random number generator Optional; default = 1234567.')
    parser.add_option('','--minsetsize',type='int',metavar='INT',default=2,\
                          help='Minimum number of genes annotated to a term. Optional; default = 2')
    parser.add_option('','--maxsetsize',type='int',metavar='INT',default=400,\
                          help='Maximum number of genes annotated to a term. Optional; default = 400')
    parser.add_option('','--samplesize',type='int',metavar='INT',default=10,\
                          help='For the negative set, sample X times the number of positives. Optional; default = 10')
    
    # General Options
    (opts, args) = parser.parse_args()
    print('ARGUMENTS ARE:',opts)
    
    if opts.network == None:
        parser.print_help()
        sys.exit('\nERROR: network file required')

    if opts.probs == None:
        if opts.annotations == None:
            parser.print_help()
            sys.exit('\nERROR: annotations file required')

        if opts.ontology == None:
            parser.print_help()
            sys.exit('\nERROR: ontology (OBO) file required')

        if opts.functions == None:
            parser.print_help()
            sys.exit('\nERROR: functions file required')
    
    if opts.outprefix == None:
        parser.print_help()
        sys.exit('\nERROR: output prefix required')  

    # seed the random number generator
    random.seed(opts.randseed)

    return opts


if __name__=='__main__':
    main(sys.argv)
