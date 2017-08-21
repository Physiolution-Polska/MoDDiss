# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 14:19:29 2017
@author: gbanach
"""

def C_zero_order(k_0, t):
    """
    k_0- parameter of the model
    t - an 1-D np.array of data corresponding to the time elapsed from the beginning of the experiment
    
    Gurny R, Doelker E, Peppas NA. Modelling of sustained release
    of water-soluble drugs from porous, hydrophobic polymers. Biomaterials. 1982;3:27–32.
    bibtexkey: gurny1982
    """    
    Theor_C_zero_order = k_0*t
    
    return Theor_C_zero_order


def C_zero_order_T_lag(k_0, T_lag, t):
    """
    k_0, T_lag- parameters of the model
    t - an 1-D np.array of data corresponding to the time elapsed from the beginning of the experiment
    
    Borodkin S, Tucker FE. Linear drug release from laminated
    hydroxypropyl cellulose-polyvinyl acetate films. J Pharm Sci.1975;64:1289–94.
    bibtexkey: borodkin1975
    """   
    Theor_C_zero_order_T_lag = k_0*(t - T_lag)
    
    return Theor_C_zero_order_T_lag


def C_zero_order_F0( k_0, F_0, t):
    """
    k_0, F_0- parameters of the model
    t - an 1-D np.array of data corresponding to the time elapsed from the beginning of the experiment
    Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33.
    bibtexkey: costa2001
    """    
    Theor_C_zero_order_F_0 = F_0 + k_0*t
    
    return Theor_C_zero_order_F_0
