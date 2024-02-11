#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 11:13:34 2024

@author: rubenmitchell
"""
#%% import callbacks
from call_back_classes import IncumbentCallback, CustomCutCallback, HeuristicsCallback
#%%

def solve(mdl, warmstart = False, inc = False, c_cut = False, mincut = False, comb = False, heu = False):
    
    if warmstart:
        mdl.build_warmstart()
    
    if inc: 
        cb_incumbent = mdl.model_instance.register_callback(IncumbentCallback)
        cb_incumbent.mdl = mdl
       
    
    if c_cut:
        cb_custom_cut = mdl.model_instance.register_callback(CustomCutCallback)
        cb_custom_cut.mdl = mdl
        if mincut:
            cb_custom_cut.do_mincut = True
        if comb:
            cb_custom_cut.do_comb = True
    if heu:
        cb_heu = mdl.model_instance.register_callback(HeuristicsCallback)
        cb_heu.mdl = mdl
        
    
    mdl.solve()
    
 
        
     
       