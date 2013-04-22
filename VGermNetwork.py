# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 11:55:56 2013

@author: hok1
"""

import networkx as nx
import numpy as np
import csv

class GeneNetwork:
    def __init__(self, filename):
        geneDists = np.loadtxt(filename, dtype={'names': ('germ1', 'germ2', 'dist'),
                                                'formats': ('S30', 'S30', np.integer)},
                               skiprows=1, delimiter='\t')
        for geneDist in geneDists:
            geneDist[0] = geneDist[0].strip()
            geneDist[1] = geneDist[1].strip()
        self.geneNet = nx.Graph()
        self.geneNet.add_weighted_edges_from(geneDists)

    def island(self, threshold=50):
        g2 = nx.Graph()
        for germ1, germ2, edata in self.geneNet.edges(data=True):
            if edata['weight'] < threshold:
                g2.add_edge(germ1, germ2, edata)
        return g2
        
    def subgraphs(self, threshold=50):
        return nx.connected_component_subgraphs(self.island(threshold=threshold))
        
def writeSubGraphsToFile(inputfile, outputfile, threshold=50):
    gn = GeneNetwork(inputfile)
    g2 = gn.subgraphs(threshold=threshold)
    
    out = open(outputfile, 'wb')
    writer = csv.writer(out)
    germs = []
    for sub in g2:
        writer.writerow(sub.nodes())
        for node in sub.nodes():
            if not (node in germs):
                germs.append(node)
    for germ in gn.geneNet.nodes():
        if not (germ in germs):
            germs.append(germ)
            writer.writerow([germ])
        
    out.close()
