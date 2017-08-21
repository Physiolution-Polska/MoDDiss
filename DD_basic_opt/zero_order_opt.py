#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 23:26:41 2017

@author: edward
"""

from scipy.optimize import differential_evolution

def Find_PAR_DEv_zero_order (Texp, Cexp, k_0min=0., k_0max=100.):
    
    """A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter for zero order model by the minimazation of residual sum of squares RSS of experimentally estimated and theoretically calculated values of drug concentration.
    
    Reference to zero order model:
    Gurny R, Doelker E, Peppas NA. Modelling of sustained release of water-soluble drugs from porous, hydrophobic polymers. Biomaterials. 1982;3:27–32.

    bibtexkey: gurny1982"""
    
    
    k_0_bounds= [(k_0min, k_0max)]
    
    def RSS_ZO (k_0):
        
        Theor_C_zero_order = k_0*Texp
        RSS= sum((Cexp-Theor_C_zero_order)**2)
        return RSS
    
    DEv_result_ZO = differeOCntial_evolution(RSS_ZO, bounds = k_0_bounds , maxiter=10000)
    PAR_zero_order = {'k_0': DEv_result_ZO.x[0]}
    return PAR_zero_order

def Find_PAR_DEv_zero_order_T_lag (Texp, Cexp, k_0min=0., k_0max=10., T_lag_min = -10. , T_lag_max = 10.):
    
    """A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter for zero order model by the minimazation of residual sum of squares RSS of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to zero order with time lag model: 
    Borodkin S, Tucker FE. Linear drug release from laminated hydroxypropyl cellulose-polyvinyl acetate films. J Pharm Sci.1975;64:1289–94.

    bibtexkey: borodkin1975"""
    
    k_0_T_lag_bounds = [(k_0min, k_0max) , (T_lag_min, T_lag_max)]
    
    
    def RSS_ZO_T_lag (X_ZOTlag):
        
        Theor_C_zero_order_T_lag = X_ZOTlag[0] * (Texp - X_ZOTlag[1])
        RSSzotlag= sum((Cexp-Theor_C_zero_order_T_lag)**2)
        return RSSzotlag
    
    
    DEv_result_ZO_T_lag = differential_evolution(RSS_ZO_T_lag, bounds = k_0_T_lag_bounds , maxiter=100000000, strategy='best1bin', popsize=50, mutation=(0.9, 1.8))
    PAR_zero_order_T_lag = {'k_0': DEv_result_ZO_T_lag.x[0], 'T_lag': DEv_result_ZO_T_lag.x[1]}
    return PAR_zero_order_T_lag
    