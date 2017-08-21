#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 12 00:15:10 2017

@author: edward
"""

def C_Higuchi (k_H,t):
    """
    t - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    k_H - parameter of the model
    
    Higuchi T. Rate of release of medicaments from ointment bases
    containing drugs in suspension. J Pharm Sci. 1961;50:874–5.
    
    bibtexkey: higuchi1961
    """
    C_theor_Higuchi= k_H*(t**0.5)
    
    return C_theor_Higuchi

def C_Higuchi_T_lag (k_H, T_lag, t):
    """
    t - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    k_H, T_lag - parameters of the model
    
    Tarvainen M, Peltonen S, Mikkonen H, Elovaara M, Tuunainen
    M, Paronen P et al. Aqueous starch acetate dispersion as a novel
    coating material for controlled release products. J Control
    Release. 2004;96:179–91.
    
    bibtexkey: tarvainen2004
    """
    C_theor_Higuchi_T_lag= k_H*((t-T_lag)**0.5)
    
    return C_theor_Higuchi_T_lag

def C_Higuchi_F0 (k_H, F0, t):
    """
    t - an 1-D np.array of experimental data corresponding to the time elapsed from the beginning of the experiment.
    k_H, F0 - parameters of the model
    
    Ford JL, Mitchell K, Rowe P, Armstrong DJ, Elliott PNC,
    Rostron C et al. Mathematical modelling of drug release from
    hydroxypropylmethylcellulose matrices: effect of temperature.
    Int J Pharm. 1991;71:95–104
    
    bibtexkey: ford1991
    """
    C_theor_Higuchi_F0= F0+ k_H*(t**0.5)
    
    return C_theor_Higuchi_F0
