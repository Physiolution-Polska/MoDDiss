#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 23:26:41 2017

@author: edward
"""

from scipy.optimize import differential_evolution

def Find_PAR_DEv_zero_order (Texp, Cexp, k_0min=0., k_0max=100.):
    
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_0 - the optimal zero order model parameter
    k_0min - an estimated minimal value of the k_0 parameter, which defines a lower boundary for the DEv algorithm
    k_0max - an estimated maximal value of the k_0 parameter, which defines an upper boundary for the DEv algorithm
       
       A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares RSS of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to zero order model:
    Gurny R, Doelker E, Peppas NA. Modelling of sustained release
    of water-soluble drugs from porous, hydrophobic polymers. Biomaterials. 1982;3:27–32.

    bibtexkey: gurny1982"""
    
    
    k_0_bounds= [(k_0min, k_0max)]
    
    def RSS_ZO (k_0):
        
        Theor_C_zero_order = k_0*Texp
        RSS= sum((Cexp-Theor_C_zero_order)**2)
        return RSS
    
    DEv_result_ZO = differential_evolution(RSS_ZO, bounds = k_0_bounds , maxiter=10000)
    PAR_zero_order = {'k_0': DEv_result_ZO.x[0]}
    return PAR_zero_order

def Find_PAR_DEv_zero_order_T_lag (Texp, Cexp, k_0min=0., k_0max=10., T_lag_min = -10. , T_lag_max = 100.):
    
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_0, T_lag - the optimal parameters for the first order model with Tlag 
    k_0min - an estimated minimal value of the k_0 parameter, which defines a boundary for the DEv algorithm
    k_0max - an estimated maximal value of the k_0 parameter, which defines a boundary for the DEv algorithm
    T_lag_min - an estimated minimal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    T_lag_max - an estimated maximal value of the T_lag parameter, which defines a boundary for the DEv algorithm
       
       A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares RSS of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to zero order with time lag model:
    Borodkin S, Tucker FE. Linear drug release from laminated
    hydroxypropyl cellulose-polyvinyl acetate films. J Pharm Sci.1975;64:1289–94.

    bibtexkey: borodkin1975"""
    
    k_0_T_lag_bounds = [(k_0min, k_0max) , (T_lag_min, T_lag_max)]
    
    
    def RSS_ZO_T_lag (X_ZOTlag):
        
        Theor_C_zero_order_T_lag = X_ZOTlag[0] * (Texp - X_ZOTlag[1])
        RSSzotlag= sum((Cexp-Theor_C_zero_order_T_lag)**2)
        return RSSzotlag
    
    
    DEv_result_ZO_T_lag = differential_evolution(RSS_ZO_T_lag, bounds = k_0_T_lag_bounds , maxiter=100000000, strategy='best1bin', popsize=50, mutation=(0.9, 1.8))
    PAR_zero_order_T_lag = {'k_0': DEv_result_ZO_T_lag.x[0], 'T_lag': DEv_result_ZO_T_lag.x[1]}
    return PAR_zero_order_T_lag
    
def Find_PAR_DEv_zero_order_F0 (Texp, Cexp, k_0min=0., k_0max=10.,  F0_min=0., F0_max=100):
    
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_0, T_lag - the optimal parameters for the first order model with Tlag 
    k_0min - an estimated minimal value of the k_0 parameter, which defines a boundary for the DEv algorithm
    k_0max - an estimated maximal value of the k_0 parameter, which defines a boundary for the DEv algorithm
    F0_min - an estimated minimal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    F0_max - an estimated maximal value of the T_lag parameter, which defines a boundary for the DEv algorithm
       
       A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares RSS of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to zero order with time lag model:
    Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33.

    bibtexkey: costa2001"""
    
    k_0_F0_bounds = [(k_0min, k_0max) , (F0_min, F0_max)]
    
    
    def RSS_ZO_F0 (X_ZOF0):
        
        Theor_C_zero_order_F0 = X_ZOF0[0] * Texp + X_ZOF0[1]
        RSSzoF0= sum((Cexp-Theor_C_zero_order_F0)**2)
        return RSSzoF0
    
    
    DEv_result_ZO_F0 = differential_evolution(RSS_ZO_F0, bounds = k_0_F0_bounds , maxiter=100000000, strategy='best1bin', popsize=50, mutation=(0.9, 1.8))
    PAR_zero_order_F0 = {'k_0': DEv_result_ZO_F0.x[0], 'F0':DEv_result_ZO_F0.x[1]}
    return PAR_zero_order_F0