#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 11:16:03 2024

@author: rubenmitchell
"""
#%% import libraries
from cplex.callbacks import LazyConstraintCallback, UserCutCallback, HeuristicCallback
from docplex.mp.callbacks.cb_mixin import ConstraintCallbackMixin
import numpy as np

from heuristic_functions import heuristics

#%% Incumbent Callback:

class IncumbentCallback(ConstraintCallbackMixin,LazyConstraintCallback):
    def __init__(self, env):
        LazyConstraintCallback .__init__(self, env)
        ConstraintCallbackMixin.__init__(self)
        
        
    def __call__(self):
        edges = []
        current_sol = self.make_solution_from_vars(self.mdl.x.values())
        for i in self.mdl.x:
            if current_sol.get_value(self.mdl.x[i]) > 0:
                edges.append(i)
        subtours = get_subtours_from_int_sol(edges)
        if len(subtours) > 0:
            for j in subtours:
                g = (sum(self.mdl.x[i] for i in j) <= len(j)-1)
                cpx_ct = self.linear_ct_to_cplex(g)
                self.add(cpx_ct[0],cpx_ct[1],cpx_ct[2])
                
                
#%% Custom Cut Callbacks 

class CustomCutCallback(ConstraintCallbackMixin,UserCutCallback):
     def __init__(self, env):
         UserCutCallback .__init__(self, env)
         ConstraintCallbackMixin.__init__(self)  
         self.do_mincut = False
         self.do_comb = False 
         
     def __call__(self):
         if not self.do_mincut:
             adj = [[] for i in self.mdl.data.V]
             solutions = self.get_values()
             j = 0
             for i in self.mdl.x:
                 if solutions[j] > 0:
                     adj[i[0]].append(i[1])
                     adj[i[1]].append(i[0])
                 j = j+1 # had to use this slight hack as CPLEX returns solutions as single list
             connected_comps = get_connected(adj)
             if len(connected_comps) > 1:
                 connected_comps.pop(0)
                 for k in connected_comps: # adds constraints
                     k.sort()
                     g = (sum(self.mdl.x[k[i], k[j]] for i in range(len(k)) for j in range(i+1,len(k))) <= len(k)-1)
                     cpx_ct = self.linear_ct_to_cplex(g)
                     self.add(cpx_ct[0],cpx_ct[1],cpx_ct[2]) 
         else:
             adj = [[] for i in self.mdl.data.V]
             solutions = self.get_values()
             j = 0
             for i in self.mdl.x:
                 if solutions[j] > 0:
                     adj[i[0]].append(i[1])
                     adj[i[1]].append(i[0])
                 j = j+1 # had to use this slight hack as CPLEX returns solutions as single list
             connected_comps = get_connected(adj)
             if len(connected_comps) == 1:
                  n = self.mdl.data.n
                  p = 0
                  compatible_array = np.zeros((n,n))
                  for (i,j) in self.mdl.data.A:
                      compatible_array[i][j] = compatible_array[j][i] = solutions[p]
                      p = p + 1
                  w, cut = globalMinCut(compatible_array)
                  if w < 2.0 - 1e-9: # ask about this? 
                      tuc = list(filter(lambda idx: (idx not in cut), range(0, n)))
                      expr = sum(self.mdl.x[min(i,j), max(i,j)] for i in cut for j in tuc) >= 2
                      cpx_ct = self.linear_ct_to_cplex(expr)
                      self.add(cpx_ct[0],cpx_ct[1],cpx_ct[2])   
             else:
                 connected_comps.pop(0)
                 for k in connected_comps: # adds constraints
                     k.sort()
                     g = (sum(self.mdl.x[k[i], k[j]] for i in range(len(k)) for j in range(i+1,len(k))) <= len(k)-1)
                     cpx_ct = self.linear_ct_to_cplex(g)
                     self.add(cpx_ct[0],cpx_ct[1],cpx_ct[2]) 
         if self.do_comb:
             #print("entering comb")
             add_comb_inequality(self)
             #add_comb_inequality2(self)
             #print("exited comb")


#%% Heuristic callback 
class HeuristicsCallback(ConstraintCallbackMixin, HeuristicCallback):
    
    def __init__(self, env):
        HeuristicCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)
        self.counter = 0
    
    
    def __call__(self):
        self.counter += 1
        #print("running")
        if self.counter == 1:
            self.counter = 0
            x_sol = self.make_solution_from_vars(self.mdl.x.values())
            adj = [[] for i in self.mdl.data.V]
            arcs = []
            for (i, j) in self.mdl.data.A:
                val = x_sol.get_value(self.mdl.x[i,j])
                if val > -0:
                    arcs.append(((i,j),val))
            #print("raw = "+ str(arcs))
            sorted_arcs = [i[0] for i in sorted(arcs, key = lambda x:x[1], reverse = True)]
            path = clean(sorted_arcs, self.mdl.data.n)
            #print("cleaned = "+ str(path))
            heu = heuristics(self.mdl)
            for i in path:
                adj[i[1]].append(i[0])
                adj[i[0]].append(i[1])
            heu.greedy_tsp_with_partial_edges(path, get_connected(adj))
             
            heu.two_opt()
            #print("big out")
            cost = heu.calculate_tour_cost()
               
            if cost < self.get_incumbent_objective_value():
                tb = heu.current_tour
                weights = []
                c = []
                for (i,j) in self.mdl.x:
                    c.append(self.mdl.x[i,j].name)
                    if (i,j) in tb:
                        weights.append(1)
                    else:
                        weights.append(0)
                self.set_solution([c, weights])#, cost)
                print("Incumbent updated with cost = "+ str(cost))# + " and bi = " + str(self.get_incumbent_objective_value()))
# =============================================================================
#             else:
#                 print("failed")
# =============================================================================

#%% Functions to find properties of current solutions 

def get_subtours_from_int_sol(arcs):
    s = []
    while arcs != []: 
        x = arcs[0][0]
        y = arcs[0][1]
        tour = [arcs[0]]
        arcs.pop(0)
        while x != y:
            for i in arcs:
                if x in i:
                    tour.insert(0,i)
                    x = i[0] if i[0] != x else i[1]
                    arcs.remove(i)
                    break
                elif y in i:
                    tour.append(i)
                    y = i[0] if i[0] != y else i[1] 
                    arcs.remove(i)
                    break
        s.append(tour)
    s.pop(0)
    return s

def get_connected(adj):
    
    def DFS(temp, v, visited, adj):
     
            # Mark the current vertex as visited
            visited[v] = True
            # Store the vertex to list
            temp.append(v)
            # Repeat for all vertices adjacent
            # to this vertex v
            for i in adj[v]:
                if visited[i] == False:
                    # Update the list
                    temp = DFS(temp, i, visited, adj)
            return temp
    visited = []
    cc = []
    for i in range(len(adj)):
        visited.append(False)
    for v in range(len(adj)):
        if visited[v] == False:
            temp = []
            cc.append(DFS(temp, v, visited, adj))
    return cc

def globalMinCut(A):
    A = np.copy(A)
    inf = 1e9
    n = A.shape[0]
    co = [[i] for i in range(n)]
    opt_w = inf
    opt_cut = []
    for ph in range(1, n):
        w = np.copy(A[0])
        s = t = 0
        for it in range(n - ph):
            w[t] = -inf
            s = t
            t = w.argmax()
            w += A[t]
        if w[t] - A[t, t] < opt_w:
            opt_w = w[t] - A[t, t]
            opt_cut = co[t]
        co[s].extend(co[t])
        A[s, :] += A[t, :]
        A[:, s] = A[s, :]
        A[0, t] = -inf
    return opt_w, opt_cut

def clean(edges, n):
    adj = [[] for i in range(n)]
    clean = []
    for edge in edges: 
        if len(adj[edge[0]]) < 2 and len(adj[edge[1]]) < 2:
            valid = True
            for i in get_connected(adj): 
                if edge[0] in i and edge[1] in i:
                    valid = False  
                    break
            if valid:
                adj[edge[0]].append(edge[1])
                adj[edge[1]].append(edge[0])
                clean.append(edge)
    return clean   
 
def add_comb_inequality(self):
    solutions = self.get_values()
    values = {(i,j): 0 for (i,j) in self.mdl.x}
    frac_edges = []
    frac_adj = [[] for i in self.mdl.data.V]
    teeth_potentials = []
    j = 0
    for i in self.mdl.x:
        values[i] = solutions[j]
        if 0.999 > solutions[j] > -0:
            frac_edges.append(i)
            frac_adj[i[0]].append(i[1])
            frac_adj[i[1]].append(i[0])
        elif solutions[j] > 0.999:
            teeth_potentials.append(i)
        j = j+1
    comps = [x for x in get_connected(frac_adj) if len(x)>2]
    comp_edge = []
    for i in comps:
        for j in frac_edges:
            if any(x in i for x in j):
                comp_edge.append(j)
                break
    
    
    def dfs_cycle(u, p, visited, parent, cycles):
        if visited[u] == 2:
            return
        if visited[u] == 1:
            v = []
            cur = p
            v.append(cur)
            while cur != u:
                cur = parent[cur]
                v.append(cur)
            cycles.append(v)
            return
        parent[u] = p
        visited[u] = 1
        for v in frac_adj[u]:
            if v == parent[u]:
                continue
            dfs_cycle(v, u, visited, parent, cycles)
        visited[u] = 2
    
    cycles = []
    
    for i in comp_edge:
        holder = []
        dfs_cycle(i[0], i[1], [0 for i in self.mdl.data.V], [0 for i in self.mdl.data.V], holder)
        sorted_cycles = sorted(holder, key=len, reverse=True)
        cycles.append(sorted_cycles)
            
    def find_teeth(cycle_comps, teeth):
        used = []
        found_teeth = []
        for i in cycle_comps: 
            for j in teeth:
                if i == j[0] and j[1] not in used:
                    used.append(j[1])
                    if j not in found_teeth:
                        found_teeth.append(j) 
                    break
                elif i == j[1] and j[0] not in used:
                    used.append(j[0])
                    if j not in found_teeth:
                        found_teeth.append(j) 
                    break 
        return found_teeth
    
    for i in range(len(cycles)):
       
        for j in range(len(cycles[i])):
            
                teeth = find_teeth(cycles[i][j], teeth_potentials)
                
                if len(teeth)%2 == 1 and len(teeth)> 2:
                    lhs = sum(values[tooth] for tooth in teeth)
                    edges = []
                    for q in range(len(cycles[i][j])):
                        #for s in range(q+1,len(cycles[i][j])):
                                
                        mi = min(cycles[i][j][q],cycles[i][j][q-1])
                        ma = max(cycles[i][j][q],cycles[i][j][q-1])
                        #if (mi,ma) in frac_edges:
                        lhs += values[mi,ma]
                        edges.append((mi,ma))
                    rhs = len(cycles[i][j]) + (len(teeth)-1)/2
                    if lhs > rhs:
                        expr = sum(self.mdl.x[z] for z in teeth) + sum(self.mdl.x[c] for c in edges) <= rhs
                        cons = self.linear_ct_to_cplex(expr)
                        self.add(cons[0],cons[1],cons[2])
                        break
        
                       
                         

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    