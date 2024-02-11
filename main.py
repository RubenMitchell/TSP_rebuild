#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 09:58:28 2024

@author: rubenmitchell
"""

#%% Import libraries
import matplotlib.pyplot as plt

#%% Import scripts
from data_classes import rnd_data, tsplib
from model_class import model
from solver_functions import solve

#%% Plotter
def plot(mdl):
    plt.figure()
    for i in mdl.data.loc:
        plt.scatter(mdl.data.loc[i][0], mdl.data.loc[i][1], c = 'black')
        plt.annotate(i, (mdl.data.loc[i][0], mdl.data.loc[i][1]))
    if mdl.model_instance.solve_status != None:
        for (i,j) in mdl.data.A:
            if mdl.x[i,j].solution_value > 0:
                plt.plot([mdl.data.loc[i][0],mdl.data.loc[j][0]], [mdl.data.loc[i][1],mdl.data.loc[j][1]], c = 'red')
    plt.axis([0, mdl.data.width, 0, mdl.data.height])
    plt.grid()
    plt.figtext(0.5,0.06,"Solution value = "+str(round(mdl.solution.objective_value,3))+" with MIP gap: "+str(round(mdl.solution.solve_details.gap*100,2))+"%", wrap=True, horizontalalignment='center', fontsize=12)
    fig = plt.gcf()
    fig.set_size_inches(8, 8)
#%% Main 

random_data = rnd_data(400, 100, 0)
TSPlib_data = tsplib('d493.tsp')

datatset = random_data
mdl = model("TSP", datatset)

warmstart = True
do_incumbent_cb = True
do_customcut_cb = True
do_mincut = True  #only functions if doing custom cut callback 
do_comb = True    #only functions if doing custom cut callback 
do_heu_cb = True 

solve(mdl, warmstart, do_incumbent_cb, do_customcut_cb, do_mincut, do_comb, do_heu_cb)
plot(mdl)
print(mdl.solution.objective_value)
print(mdl.solution.solve_details.gap)
