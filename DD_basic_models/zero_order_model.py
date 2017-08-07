# -*- coding: utf-8 -*-
"""
Created on Sat Aug  5 14:19:29 2017

@author: gbanach
"""

def zero_order(t, k_0):
    """
    Gurny R, Doelker E, Peppas NA. Modelling of sustained release
    of water-soluble drugs from porous, hydrophobic polymers. Biomaterials. 1982;3:27–32.

    bibtexkey: gurny1982
    """    
    F = k_0*t
    
    return F


def zero_order_T(t, k_0, T_lag):
    """
    Borodkin S, Tucker FE. Linear drug release from laminated
    hydroxypropyl cellulose-polyvinyl acetate films. J Pharm Sci.1975;64:1289–94.

    bibtexkey: borodkin1975
    """   
    F = k_0*(t - T_lag)
    
    return F


def zero_order_F0(t, k_0, F_0):
    """
    Costa P, Sousa Lobo JM. Modeling and comparison of dissolution profiles. Eur J Pharm Sci. 2001;13:123–33.

    bibtexkey: costa2001
    """    
    F = F_0 + k_0*t
    
    return F
