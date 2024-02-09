#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 09:46:15 2024

@author: rubenmitchell
"""

import numpy as np
rnd = np.random

class heuristics: 
    def __init__(self, model):
        self.mdl = model
        self.current_tour = []
    
    def calculate_tour_cost(self):
        cost = 0
        for i in self.current_tour:
            #print(i)
            cost += self.mdl.data.dists[i]
        return cost
        
    def greedy_tsp_with_partial_edges(self, partial_edges, connected_comps):
        unsaturated = []
        self.current_tour = partial_edges
        visited_count = [0 for i in self.mdl.data.V]
        #print("partial edges = "+str(partial_edges))
        for i in partial_edges: 
            visited_count[i[0]] += 1
            visited_count[i[1]] += 1
        for i in range(len(visited_count)):
            if visited_count[i] == 1:
                unsaturated.append(i)
        while len(unsaturated) > 2:
            cur_node = unsaturated[0]
            unsaturated.pop(0)
            options = []
            current_comp = None
            for i in connected_comps:
                if cur_node in i:
                    current_comp = i
                    break
            for i in unsaturated:
                if i not in current_comp:
                    options.append(i)
            
            nearest_node = min(options, key=lambda node: self.mdl.data.dists[min(cur_node, node), max(cur_node, node)])
            unsaturated.remove(nearest_node)
            #logic for merging connected comps:
            joining_comp = []
            for i in connected_comps:
                if nearest_node in i:
                    joining_comp = i
                    break
            connected_comps.append(current_comp+joining_comp)
            connected_comps.remove(joining_comp)
            connected_comps.remove(current_comp)
            self.current_tour.append((min(cur_node,nearest_node), max(cur_node, nearest_node)))
           
        x = unsaturated[0]
        y = unsaturated[1]
        self.current_tour.append((min(x,y),max(x,y)))
        
    
    def two_opt_calc(self, A, B, C, D):
        return self.mdl.data.dists[min(A,C),max(A,C)] + self.mdl.data.dists[min(B,D),max(B,D)] - self.mdl.data.dists[min(A,B),max(A,B)] - self.mdl.data.dists[min(C,D),max(C,D)]    

    def two_opt(self):
        tour_nodes = [self.current_tour[0][0],self.current_tour[0][1]]
        while len(tour_nodes)<len(self.current_tour):
            for i in self.current_tour:
                if tour_nodes[-1] in i and tour_nodes[-2] not in i:
                    if i[0] == tour_nodes[-1]:
                        tour_nodes.append(i[1])
                    else:
                        tour_nodes.append(i[0])
        if tour_nodes[-1] == tour_nodes[0]:
            tour_nodes.pop(-1)
        tour_len = len(tour_nodes)
        improved = True
        while improved:
            improved = False
            for i in range(0, tour_len-3):
                for j in range(i + 2, tour_len-1):
                    x = self.two_opt_calc(tour_nodes[i - 1], tour_nodes[i], tour_nodes[j -1], tour_nodes[j])
                    if x < 0:
                        tour_nodes[i:j] = reversed(tour_nodes[i:j])
                        improved = True
        self.current_tour = []
        tour_nodes.append(tour_nodes[0])
        for i in range(len(tour_nodes)-1):
            x = tour_nodes[i]
            y = tour_nodes[i+1]
            self.current_tour.append((min(x,y),max(x,y)))
    
    
        
            
        