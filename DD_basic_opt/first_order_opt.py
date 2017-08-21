#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 17:56:23 2017

@author: edward
"""

from scipy.optimize import differential_evolution
import numpy as np

def Find_PAR_DEv_first_order (Texp, Cexp, k_1min=0., k_1max=100.):
    
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_1 - the optimal first order model parameter
    k_1min - an estimated minimal value of the k_1 parameter, which defines a boundary for the DEv algorithm
    k_1max - an estimated maximal value of the k_1 parameter, which defines a boundary for the DEv algorithm
       
       A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares RSS of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to the first order model:
    Polli JE, Rekhi GS, Augsburger LL, Shah VP. Methods to compare dissolution 
    profiles and a rationale for wide dissolution specifications for metoprolol
    tartrate tablets. J Pharm Sci.1997;86:690–700.

    bibtexkey: polly1997"""
    
    
    k_1_bounds= [(k_1min, k_1max)]
    
    def RSS_FO (k_1):
        
        Theor_C_first_order = 100*(1-np.exp(-k_1*Texp))
        RSSFO= sum((Cexp-Theor_C_first_order)**2)
        return RSSFO
    
    DEv_result_FO = differential_evolution(RSS_FO, bounds = k_1_bounds , maxiter=10000)
    PAR_First_order = {'k_1': DEv_result_FO.x[0]}
    return PAR_First_order

def Find_PAR_DEv_first_order_T_lag (Texp, Cexp, k_1min=0., k_1max=100., T_lag_min=0., T_lag_max=50):
    
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_1, T_lag - the optimal parameters for the first order model with Tlag 
    k_1min - an estimated minimal value of the k_1 parameter, which defines a boundary for the DEv algorithm
    k_1max - an estimated maximal value of the k_1 parameter, which defines a boundary for the DEv algorithm
    T_lag_min - an estimated minimal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    T_lag_max - an estimated maximal value of the T_lag parameter, which defines a boundary for the DEv algorithm
       
       A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares RSS of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to the first order model with T_lag:
    Phaechamud T, Pitaksantayothin K, Kositwattanakoon P, Seehapong P, Jungvivatanavong S.
    Sustainable release of propranolol hydrochloride tablet using chitin as 
    press-coating material.Silpakorn Univ Int J. 2002;2:147–59.

    bibtexkey: phaechamud2002"""
    
    
    k_1_T_lag_bounds= [(k_1min, k_1max), (T_lag_min, T_lag_max)]
    
    def RSS_FO_T_lag (XFOTlag):
        
        Theor_C_First_order_T_lag = 100*(1- np.exp(-XFOTlag[0] * (Texp- XFOTlag[1])))
        RSSFOTLAG = sum((Cexp -Theor_C_First_order_T_lag)**2)
        return RSSFOTLAG
    
    DEv_result_FOTlag = differential_evolution(RSS_FO_T_lag, bounds = k_1_T_lag_bounds , maxiter=10000)
    PAR_First_order_T_lag = {'k_1': DEv_result_FOTlag.x[0], 'T_lag': DEv_result_FOTlag.x[1]}
    return PAR_First_order_T_lag

def Find_PAR_DEv_first_order_F_max (Texp, Cexp, k_1min=0., k_1max=100., F_max_min=0., F_max_max=1000):
    
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_1, F_max - the optimal parameters for the first order model with F_max
    k_1min - an estimated minimal value of the k_1 parameter, which defines a boundary for the DEv algorithm
    k_1max - an estimated maximal value of the k_1 parameter, which defines a boundary for the DEv algorithm
    F_max_min - an estimated minimal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    F_max_max - an estimated maximal value of the T_lag parameter, which defines a boundary for the DEv algorithm
       
       A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares RSS of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to the first order model with F_max :
    Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and 
    acceptance sampling rule based on profile modeling and principal component 
    analysis. J Biopharm Stat. 1997;7:423–39.

    bibtexkey: tsong1997"""
    
    
    k_1_F_max_bounds= [(k_1min, k_1max), (F_max_min, F_max_max)]
    
    def RSS_FO_F_max (XFOFmax):
        
        Theor_C_First_order_F_max = XFOFmax[1]*(1- np.exp(-XFOFmax[0] * Texp))
        RSSFOFmax = sum((Cexp -Theor_C_First_order_F_max)**2)
        return RSSFOFmax
    
    DEv_result_FOFmax = differential_evolution(RSS_FO_F_max, bounds = k_1_F_max_bounds , maxiter=10000)
    PAR_First_order_F_max = {'k_1': DEv_result_FOFmax.x[0], 'F_max': DEv_result_FOFmax.x[1]}
    return PAR_First_order_F_max

def Find_PAR_DEv_first_order_F_max_T_lag (Texp, Cexp, k_1min=0., k_1max=100., F_max_min=0., F_max_max=1000, T_lag_min=0., T_lag_max=50):
    
    """
    Texp - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    Cexp - an 1-D np.array of experimental data corresponding to the drug concentration in the fixed moments defined by Texp.
    k_1, F_max, T_lag - the optimal parameters for the first order model with F_max and T_lag
    k_1min - an estimated minimal value of the k_1 parameter, which defines a boundary for the DEv algorithm
    k_1max - an estimated maximal value of the k_1 parameter, which defines a boundary for the DEv algorithm
    F_max_min - an estimated minimal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    F_max_max - an estimated maximal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    T_lag_min - an estimated minimal value of the T_lag parameter, which defines a boundary for the DEv algorithm
    T_lag_max - an estimated maximal value of the T_lag parameter, which defines a boundary for the DEv algorithm
       
       A differential evolution algorithm ( https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html )
    is used to find an optimal parameter
    for zero order model by the minimazation of residual sum of squares RSS of experimentally estimated 
    and theoretically calculated values of drug concentration.
    
    Reference to the first order model with Fmax and Tlag:
     Berry MR, Likar MD. Statistical assessment of dissolution and drug release 
    profile similarity using a model-dependent approach. J Pharm Biomed Anal. 
    2007;45:194–200.
    
    bibtexkey: berry2007"""
    
    
    k_1_F_max_T_lag_bounds= [(k_1min, k_1max), (F_max_min, F_max_max), (T_lag_min, T_lag_max)]
    
    def RSS_FO_F_max_T_lag (XFOFmaxTlag):
        
        Theor_C_First_order_F_max_T_lag = XFOFmaxTlag[1]*(1- np.exp(-XFOFmaxTlag[0] * (Texp- XFOFmaxTlag[2])))
        RSSFOFmaxTlag = sum((Cexp -Theor_C_First_order_F_max_T_lag)**2)
        return RSSFOFmaxTlag
    
    DEv_result_FOFmaxTlag = differential_evolution(RSS_FO_F_max_T_lag, bounds = k_1_F_max_T_lag_bounds , maxiter=50000, popsize= 25)
    PAR_First_order_F_max_T_lag = {'k_1': DEv_result_FOFmaxTlag.x[0], 'F_max': DEv_result_FOFmaxTlag.x[1], 'T_lag': DEv_result_FOFmaxTlag.x[2] }
    return PAR_First_order_F_max_T_lag