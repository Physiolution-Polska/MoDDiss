#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 01:27:07 2017

@author: edward
"""

from scipy.optimize import differential_evolution

def Find_PAR_DEv_Higuchi (Texp, Cexp, k_H_min=0., k_H_max=150.):
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_H - the optimal Higuchi model parameter
    k_Hmin - an estimated minimal value of the k_H parameter, which defines a boundary for the DEv algorithm
    k_Hmax - an estimated maximal value of the k_H parameter, which defines a boundary for the DEv algorithm
    
      A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares (RSS) of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to the Higuchi model:
    Higuchi T. Rate of release of medicaments from ointment bases
    containing drugs in suspension. J Pharm Sci. 1961;50:874–5.
    
    bibtexkey: higuchi1961
    """
    k_H_bounds= [(k_H_min, k_H_max)]
    
    def RSS_Higuchi (k_H):
        
        Theor_C_Higuchi= k_H*(Texp**0.5)
        RSSHig= sum((Cexp-Theor_C_Higuchi)**2)
        
        return RSSHig

    DEv_result_Hig= differential_evolution(RSS_Higuchi, bounds= k_H_bounds, maxiter=50000, popsize= 25)
    PAR_Higuchi={'k_H':DEv_result_Hig.x[0]}
    return PAR_Higuchi

def Find_PAR_DEv_Higuchi_T_lag (Texp, Cexp, k_H_min=0., k_H_max=150., T_lag_min= 0., T_lag_max= 24.):
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_H, T_lag - the optimal Higuchi model parameters
    k_Hmin - an estimated minimal value of the k_H parameter, which defines a boundary for the DEv algorithm
    k_Hmax - an estimated maximal value of the k_H parameter, which defines a boundary for the DEv algorithm
    T_lag_min - an estimated minimal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    T_lag_max - an estimated maximal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    
      A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares (RSS) of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to the Higuchi with T_lag model:
    
        Tarvainen M, Peltonen S, Mikkonen H, Elovaara M, Tuunainen
    M, Paronen P et al. Aqueous starch acetate dispersion as a novel
    coating material for controlled release products. J Control
    Release. 2004;96:179–91.
    
    bibtexkey: tarvainen2004
    """
    k_H_T_lag_bounds= [(k_H_min, k_H_max), (T_lag_min, T_lag_max)]
    
    def RSS_Higuchi_T_lag (k_H, T_lag):
        Theor_C_HiguchiTlag= k_H*((Texp-T_lag)**0.5)
        RSSHigTlag= sum((Cexp-Theor_C_HiguchiTlag)**2)
        
        return RSSHigTlag
    
    DEv_result_Hig_T_lag= differential_evolution(RSS_Higuchi_T_lag, bounds= k_H_T_lag_bounds, maxiter=50000, popsize= 25)
    PAR_Higuchi_T_lag= {'k_H':DEv_result_Hig_T_lag[0], 'T_lag': DEv_result_Hig_T_lag[1]}
    return PAR_Higuchi_T_lag

def Find_PAR_DEv_Higuchi_F0 (Texp, Cexp, k_H_min=0., k_H_max=150., F0_min= 0., F0_max= 24.):
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_H, T_lag - the optimal Higuchi model parameters
    k_Hmin - an estimated minimal value of the k_H parameter, which defines a boundary for the DEv algorithm
    k_Hmax - an estimated maximal value of the k_H parameter, which defines a boundary for the DEv algorithm
    F0_min - an estimated minimal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    F0_max - an estimated maximal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    
      A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares (RSS) of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to the Higuchi with F0 model:
        
        Ford JL, Mitchell K, Rowe P, Armstrong DJ, Elliott PNC,
    Rostron C et al. Mathematical modelling of drug release from
    hydroxypropylmethylcellulose matrices: effect of temperature.
    Int J Pharm. 1991;71:95–104
    
    bibtexkey: ford1991
    """
    k_H_T_lag_bounds= [(k_H_min, k_H_max), (F0_min, F0_max)]
    
    def RSS_Higuchi_F0 (k_H, F0):
        Theor_C_HiguchiF0= F0+ k_H*(Texp**0.5)
        RSSHigF0= sum((Cexp-Theor_C_HiguchiF0)**2)
        
        return RSSHigF0
    
    DEv_result_Hig_F0= differential_evolution(RSS_Higuchi_F0, bounds= k_H_T_lag_bounds, maxiter=50000, popsize= 25)
    PAR_Higuchi_F0= {'k_H':DEv_result_Hig_F0[0], 'F0': DEv_result_Hig_F0[1]}
    return PAR_Higuchi_F0
    