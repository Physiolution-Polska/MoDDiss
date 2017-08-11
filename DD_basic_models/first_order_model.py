#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 15:26:58 2017

@author: edward
"""
import numpy as np

def C_first_order(k_1, t):
    """
    k_1- parameter of the model
    t - an 1-D np.array of data corresponding to the time elapsed from the beginning of the experiment
       
    Polli JE, Rekhi GS, Augsburger LL, Shah VP. Methods to compare dissolution 
    profiles and a rationale for wide dissolution specifications for metoprolol
    tartrate tablets. J Pharm Sci.1997;86:690–700.

    bibtexkey: polly1997
    """    
    Theor_C_first_order = 100*(1- np.exp(-k_1*t))
    
    return Theor_C_first_order


def C_first_order_T_lag(k_1, T_lag, t):
    """
    k_1, T_lag - parameters of the model
    t - an 1-D np.array of data corresponding to the time elapsed from the beginning of the experiment
       
    Phaechamud T, Pitaksantayothin K, Kositwattanakoon P, Seehapong P, Jungvivatanavong S.
    Sustainable release of propranolol hydrochloride tablet using chitin as 
    press-coating material.Silpakorn Univ Int J. 2002;2:147–59.

    bibtexkey: phaechamud2002
    """   
    Theor_C_first_order_T_lag = 100*(1- np.exp( - k_1*(t - T_lag)))
    
    return Theor_C_first_order_T_lag


def C_first_order_F_max( k_1, F_max, t):
    """
    k_1, F_max- parameters of the model
    t - an 1-D np.array of data corresponding to the time elapsed from the beginning of the experiment
        
    Tsong Y, Hammerstrom T, Chen JJ. Multipoint dissolution specification and 
    acceptance sampling rule based on profile modeling and principal component 
    analysis. J Biopharm Stat. 1997;7:423–39.

    bibtexkey: tsong1997
    """    
    Theor_C_first_order_F_max = F_max *(1- np.exp( -k_1* t))
    
    return Theor_C_first_order_F_max


def C_first_order_F_max_T_lag( k_1, F_max, T_lag, t):
    """
    k_1, F_max, T_lag- parameters of the model
    t - an 1-D np.array of data corresponding to the time elapsed from the beginning of the experiment
        
    Berry MR, Likar MD. Statistical assessment of dissolution and drug release 
    profile similarity using a model-dependent approach. J Pharm Biomed Anal. 
    2007;45:194–200.
    
    bibtexkey: berry2007
    """    
    Theor_C_first_order_F_max_T_lag = F_max *(1- np.exp( -k_1*(t - T_lag)))
    
    return Theor_C_first_order_F_max_T_lag