#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 09:59:32 2024

@author: rubenmitchell
"""
#%% Import libraries
import numpy as np 
rnd = np.random
import tsplib95
from pathlib import Path

#%% Random Data Class:

class rnd_data: 
    def __init__(self, n, width, seed):
        rnd.seed(seed)
        self.n = n                                                          # Quantity of vertices
        self.width = width     
        self.height = width                                                 # Dimensions of square plotting space
        self.V = range(0, n)                                                # Variable for each vertex
        self.A = [(i,j) for i in self.V for j in range (i+1, n)]            # Variable for every possible edge defined as 2D list 
        self.loc = {i:(rnd.random()*self.width,rnd.random()*self.width) for i in self.V}
        self.dists = {(i,j): np.linalg.norm(np.array(self.loc[i])- np.array(self.loc[j])) for (i,j) in self.A}



class tsplib:
    def __init__(self, filename):
        
        tsp = tsplib95.load(Path.cwd() / "TSPlib_files" / filename)
        self.n = tsp.dimension
        self.V = range(0,self.n)
        self.loc = {i:(tsp.node_coords[i+1]) for i in self.V}
        self.width = max(self.loc.values())[0]*1.05
        self.height = max(self.loc.values(), key=lambda x: x[1])[1]*1.05
        self.A = [(i,j) for i in self.V for j in range(i+1,self.n)]
        self.dists = {(i,j): tsp.get_weight(i+1,j+1) for (i,j) in self.A}