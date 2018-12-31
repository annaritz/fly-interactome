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
    print '\nReading the network from %s' %(opts.network)
    net = Graph()
    net.read(opts.network)
    print nx.info(net)
    
    # read mapping file.
    if MAPFILE != None:
        mapper = readDict(MAPFILE,FROM_COL, TO_COL)
    else:
        mapper = {}

    ## from yeast remove values '!?' (synergizer), '---' (AFFY), 'Input not found'
    ## from human remove keys/values '-','--' (Uniprot)
    tofilter = ['!?','Input not found','-','--','---']
    mapper = {k:mapper[k] for k in mapper if k not in tofilter and mapper[k] not in tofilter}
        
    # get the annotations
    print '\nReading the annotations from %s' %(ANNOTATIONS_FILE)
    ann = Annotations()
    if FUNCASSOC:
        ann.readGMT_FuncAssociate(ANNOTATIONS_FILE, mapper)
    else:
        ann.readGMT(ANNOTATIONS_FILE, mapper)
    print '\t%d functions' %(len(ann.getFunctions()))
    print '\t%d annotated genes' %(len(ann.getAnnotatedGenes()))
    
    # remove annotations for genes that are not in the network    
    print '\nKeeping annotations for only genes in the network'
    ann.keepAnnotationsForGenes(net.nodes())
    print '\t%d functions' %(len(ann.getFunctions()))
    print '\t%d annotated genes' %(len(ann.getAnnotatedGenes()))
   
    # get the GO DAG
    godag = GOdag(OBO_FILE)
    ann.applyTruePathRule(godag)
    ann.keepAnnotationsForNamespace(godag, 'biological_process')
    print 'Keeping only biological processes.'
    print '\t%d functions' %(len(ann.getFunctions()))
    print '\t%d annotated genes' %(len(ann.getAnnotatedGenes()))
    
    # get the functions that define positives
    if FUNCASSOC:
        lines = readColumns(FUNCS_FILE,3,5,6)
        posFuncs = set([goterm for logodds,pval,goterm in lines if \
                            float(logodds)>FUNCASSOC_LOGODDSTHRES and \
                            ( (FUNCASSOC_PVALTHRES >= 0.001 and pval == '<0.001')\
                                  or float(pval)<FUNCASSOC_PVALTHRES)])
        print '%d GO terms with adjusted pval < %f and log odds ratio > %f.' % \
            (len(posFuncs),FUNCASSOC_PVALTHRES,FUNCASSOC_LOGODDSTHRES)
    else:
        posFuncs = readItemSet(FUNCS_FILE,1)

    posSubgraph = godag.subgraph(posFuncs)
    posFuncs = set(posSubgraph.nodes())

    print '\nRemoving functions: (1) not in GOdag, (2) too large, (3) descendant of another function, and (4) too small, in that order.'
    print ann.genesets.keys()

    ## remove any functions that aren't in the GOdag
    posFuncs = set([f for f in posFuncs if f in ann.genesets])
    print '%d GO terms after removing terms not in GOdag' % (len(posFuncs))

    ## remove functions that are too large
    posFuncs = set([f for f in posFuncs if len(ann.genesets[f])<=MAXSETSIZE])
    print '%d GO terms after removing terms with >%d genes' % \
        (len(posFuncs),MAXSETSIZE)

    ## remove any functions with an ancestor in the list
    # for f in posFuncs:
    #     print f,ann.descriptions[f],len(ann.genesets[f])
    #     A = godag.getAncestors(f)
    #     for a in A.intersection(posFuncs):
    #         print ' ancestor',a,ann.descriptions[a],len(ann.genesets[a])      
    #     sys.exit()
    posFuncs = set([f for f in posFuncs if len(godag.getAncestors(f).intersection(posFuncs))==1])
    print '%d GO terms after removing descendants' % (len(posFuncs))

    ## remove functions that are too small
    posFuncs = set([f for f in posFuncs if len(ann.genesets[f])>=MINSETSIZE])
    print '%d GO terms after removing terms with <%d genes' % \
        (len(posFuncs),MINSETSIZE)

    # if there are fewer than 50, print them.
    if len(posFuncs)<50:
        print posFuncs

    if POSITIVE_GENEFILE != None:
        ## check coverage of all genes.
        print '\nChecking coverage of genes in positive gene file with the selected GO terms...'
        studyset = readItemSet(POSITIVE_GENEFILE,1)
        print '%d genes in studyset' % (len(studyset))
        coveredgenes = set()
        for term in posFuncs:
            interm = ann.genesets[term].intersection(studyset)
            #print '%d in term %s (%d tot)' % (len(interm),term,len(ann.genesets[term]))
            coveredgenes.update(interm)
        print '%d genes (%.2f) covered by selected GO functions' % (len(coveredgenes),len(coveredgenes)/float(len(studyset)))


    if FUNCASSOC_GENEFILE != None:
        ## check coverage of all genes.
        print '\nChecking coverage of genes in FuncAssociate study set with the selected GO terms...'
        studyset = readItemSet(FUNCASSOC_GENEFILE,1)
        print '%d genes in studyset' % (len(studyset))
        coveredgenes = set()
        for term in posFuncs:
            interm = ann.genesets[term].intersection(studyset)
            #print '%d in term %s (%d tot)' % (len(interm),term,len(ann.genesets[term]))
            coveredgenes.update(interm)
        print '%d genes (%.2f) covered by selected GO functions' % (len(coveredgenes),len(coveredgenes)/float(len(studyset)))

    ## write selected terms to file.
    output = open('%s-positive-GO-terms.txt' %(OUTPREFIX), 'w')
    if FUNCASSOC_GENEFILE != None:
        output.write('#numinstudyset\tnumgenes\tterm\tdescription\n')
        for term in posFuncs:
            output.write('%d\t%d\t%s\t%s\n' % (len(ann.genesets[term].intersection(studyset)),len(ann.genesets[term]),\
                                               term,ann.descriptions[term]))
    else:
        output.write('#numgenes\tterm\tdescription\n')
        for term in posFuncs:
            output.write('%d\t%s\t%s\n' % (len(ann.genesets[term]),term,ann.descriptions[term]))
    output.close()

    allPairs = set( [(u,v) for u in net.nodes_iter() for v in net.nodes_iter() if u<v] )
    pos = set()
    neg = set()
    pos_sets = []
    print '\nGenerating positive and negative pairs.' 
    for f in posFuncs:
        currpos = set()
        funcgenes = ann.genesets[f]
        for p in [(u,v) for u in funcgenes for v in funcgenes if u<v]:
            currpos.add(p)
        pos_sets.append(currpos)
    if INTERSECTTERMS is False:
        print 'Using the union of %d functions to generate positives' %(len(posFuncs))
        pos = set.union(*pos_sets)
    else:
        print 'Using the intersection of %d functions to generate positives' %(len(posFuncs))
        pos = set.intersection(*pos_sets)
    print '\t%d positive pairs' %(len(pos))
    
    # randomly sample non-positive pairs for negative examples
    if SAMPLESIZE*len(pos) < len(allPairs):
        neg = set( random.sample(allPairs.difference(pos), SAMPLESIZE*len(pos)) )
    else:
        print '\tWARNING: not sampling: proportion of negatives is %.4f times the size of positives' % \
            ((len(allPairs)-len(pos))/float(len(pos)))
        neg = allPairs.difference(pos)
    print '\t%d negative pairs' %(len(neg))
    
    di_etypes = di_net.getEdgeTypes()
    etypes = {}
    # convert etypes to undirected values
    for et, edges in di_etypes.iteritems():
        etypes[et] = set([(u,v) if u<v else (v,u) for (u,v) in edges])
    
    # merge edges from evidences that annotate <25 edges.
    etypes['miscellaneous'] = set()
          
    print '\nComputing the weight of individual experiments Pr(I=1|E=1).'
    etypeProbs = {}
    for et, edges in etypes.iteritems():
        
        pr_E1_I1 = float(len(pos.intersection(edges))) / len(pos) # Pr(E=1|I=1)
        pr_E0_I1 = 1.0 - pr_E1_I1                              # Pr(E=0|I=1)
        
        pr_E1_I0 = float(len(neg.intersection(edges))) / len(neg) # Pr(E=1|I=0)  
        pr_E0_I0 = 1.0 - pr_E1_I0                              # Pr(E=0|I=0)
        etypeProbs[et] = {'pr_E1_I1' : pr_E1_I1, \
                          'pr_E1_I0' : pr_E1_I0, \
                          'pr_E0_I1' : pr_E0_I1, \
                          'pr_E0_I0' : pr_E0_I0, \
                         }        
        p = float(len(pos.intersection(edges)))
        n = float(len(neg.intersection(edges)))
                    
    pr_I1 = float(len(pos))/(len(pos) + len(neg)) # Pr(I=1)
    pr_I0 = 1.0 - pr_I1                        # Pr(I=0)
    
    print '\nComputing the weight of edge with only experiment k Pr(I=1|E), where E_k=1.'
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
    
    output = open('%s-edge_type_weights.txt' %(OUTPREFIX), 'w')
    output.write('#total\tPr(I=1|E)\tPr(I=0|E)\tprob\tedge_type\n')
    for e, w in sorted(etypeWeights.iteritems(), key=lambda x: x[1]['weight'], reverse=True):
        output.write('%d\t%.5e\t%.5e\t%.5e\t%s\n' %(len(etypes[e]), w['num_I1'], w['num_I0'], w['weight'], e))
    output.close()
    
    print '\nWeighting each edge with the new edge weight.'
    output = open('%s.txt' %(OUTPREFIX), 'w')
    output.write('#tail\thead\tedge_weight\tedge_type\n')
    if WEIGHTCAP:
        cap_output = open('%s-cap%s.txt' %(OUTPREFIX, str(WEIGHTCAP).replace('.','_')), 'w')
        cap_output.write('#tail\thead\tedge_weight\tedge_type\n')
    for t,h in di_net.edges():
        on = set(di_net.edge[t][h]['types'])
        
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
        if COMPRESSED_EVIDENCE: # write one line; joine evidence types by pipe.
            output.write('%s\t%s\t%0.5e\t%s\n' % (t,h,weight,'|'.join(on)))
            ## Allow setting a threshold for the weight. 
            ## for example: weight = min(0.75, num_I1/(num_I1+num_I0))
            ## JEFF EDITED JAN 2017
            if WEIGHTCAP:
                weight = min(WEIGHTCAP, num_I1/(num_I1+num_I0)) if weight != 0 else weight
                cap_output.write('%s\t%s\t%0.5e\t%s\n' % (t,h,weight,'|'.join(on)))
        else: # write each evidence type on separate lines.
            # comment leftover from Chris:
            # I Shouldn't use "for etype in on:" here. I only want the etypes for this directed edge in di_net
            for etype in on:
                output.write('%s\t%s\t%0.5e\t%s\n' %(t, h, weight, etype))
            if WEIGHTCAP:
                weight = min(WEIGHTCAP, num_I1/(num_I1+num_I0)) if weight != 0 else weight
                cap_output.write('%s\t%s\t%0.5e\t%s\n' %(t, h, weight, etype))
    output.close()

    print '\nWrote to positive GO terms to %s-positive-GO-terms.txt' %(OUTPREFIX)
    print 'Wrote edge type probabilities to %s-edge_type_weights.txt' %(OUTPREFIX)
    print 'Wrote weighted interactome to %s.txt' %(OUTPREFIX)
    print 'DONE.'

################################################
def parse_arguments(args):
    
    usage = '''
    python weight-edges.py [OPTIONS]
    '''

    parser = OptionParser(usage=usage)
    parser.add_option('-n','--network',type='string',metavar='STR',\
                          help='protein-protein interactome (tab-delimited). First two columns are the tail and head nodes, the last column is the evidence type. Required.')
    parser.add_option('-a','--annotations',type='string',metavar='STR',\
                          help='GO annotation file.  Multi-column tab delimited with "GO term" "description" "tab-delimited list of common names". Required.')
    parser.add_option('-t','--ontology',type='string',metavar='STR',\
                          help='GO Ontology (OBO) file.  Downloaded from the Gene Ontology database. Required.')
    parser.add_option('-f','--functions',type='string',metavar='STR',\
                          help='File of GO terms that are "positives".  Originally, a 2-column tab-delimited with "GO term" "description". Now, only first column is needed. Required.')
    parser.add_option('-m','--mapper',type='string',metavar='STR',\
                          help='ID Mapping File. Tab-delimited. Optional')
    parser.add_option('','--fromcol',type='int',metavar='INT',\
                          help='Mapping from column x (systematic name), indexed starting from 1. Optional.')
    parser.add_option('','--tocol',type='int',metavar='INT',\
                          help='Mapping to column y (systematic name), indexed starting from 1. Optional.')
    parser.add_option('-o','--outprefix',type='string',metavar='STR',\
                          help='output file prefix. Required.')
    parser.add_option('','--randseed',type='int',metavar='INT',default=1234567,\
                      help='Seed for random number generator Optional; default = 1234567.')
    parser.add_option('','--minsetsize',type='int',metavar='INT',default=2,\
                          help='Minimum number of genes annotated to a term. Optional; default = 2')
    parser.add_option('','--maxsetsize',type='int',metavar='INT',default=400,\
                          help='Maximum number of genes annotated to a term. Optional; default = 400')
    parser.add_option('','--samplesize',type='int',metavar='INT',default=10,\
                          help='For the negative set, sample X times the number of positives. Optional; default = 10')
    parser.add_option('','--weightcap',type='float',metavar='FLOAT',\
                          help='Cap the edge weights at the specified value (will write results to another file with -cap0_XX appended to output prefix. Usual value is 0.75')
    
    # General Options
    (opts, args) = parser.parse_args()
    print 'ARGUMENTS ARE:',opts
    
    if opts.network == None:
        parser.print_help()
        sys.exit('\nERROR: network file required')

    if opts.annotations == None:
        parser.print_help()
        sys.exit('\nERROR: annotations file required')

    if opts.ontology == None:
        parser.print_help()
        sys.exit('\nERROR: ontology (OBO) file required')
    
    if opts.functions == None:
        parser.print_help()
        sys.exit('\nERROR: functions file required')

    #if opts.mapper == None:
    #    parser.print_help()
    #    sys.exit('\nERROR: mapping file required')  
  
    #if opts.fromcol == None:
    #    parser.print_help()
    #    sys.exit('\nERROR: from column required')  

    #if opts.tocol == None:
    #    parser.print_help()
    #    sys.exit('\nERROR: to column required')  
    
    if opts.outprefix == None:
        parser.print_help()
        sys.exit('\nERROR: output prefix required')  

    # seed the random number generator
    random.seed(opts.randseed)

    return opts


if __name__=='__main__':
    main(sys.argv)
