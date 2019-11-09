# -*- coding: utf-8 -*-

import numpy as np

#
# The parameters used in the functions below.
# 

standard_parameters = {
    'Gy_max20':2.,
    'theta_Gy':1.07, 
    'irr_opt':300.,
    'Kp':0.03, # g/m3
    'Kn':0.1,  # g/m3 
    
    'Ry20':0.025, # phy respiration
    'theta_Ry':1.045,
    
    'Cg20':8.3,   # zoo grazing
    'theta_Cg':1.045,
    
    'as':0.6,
    'Ky':0.051,
    
    'Rz20':0.04,  # zoo respiration
    'theta_Rz':1.045,
    
    'Dz':0.075,   # zoo mortal
    'ry':0.2,
    'rz':0.2,
    
    'Kd20':0.02, # det remineralization
    'theta_Kd':1.08,
    
    'Wy':0.1,
    'Wd':0.1,
    
    'orgOP':142,
    'phyNP':10.,
    'zooNP':10.,
    }


def lightattenuation(irr0, phy, H):
    if irr0 > 0.0:
        k = 10.*phy + 0.6
        irr = [irr0]
        dH = np.append(H[0], H[1:]-H[:-1])
        for ki, dHi in zip(k, dH):
            irr.append(irr[-1]*np.exp(-ki*dHi))
        return irr[1:]
    else:
        return np.zeros_like(phy)

def _tempfunc(temp, theta):
    return theta**(temp-20.)

def phy_tempfunc(temp, p):
    fac1 = _tempfunc(temp, p['theta_Gy'])
    fac2 = np.min(np.array([temp/14., np.ones_like(temp), 2.-temp/20.]), axis=0)
    return np.max(np.array([fac1, fac2]), axis=0)

def lightresponse(irr, phy, H, p):
    k = 10.*phy + 0.6
    irr_norm = irr/p['irr_opt']
    return (np.exp(1.0-irr_norm*np.exp(-k*H))
            -np.exp(1.0-irr_norm))/(k*H)

def mm(C, K):
    # micachaelismenten
    return C / (K + C)

def nutrientuptake_PN(DIP, DIN, p):
    return np.min([mm(DIP, p['Kp']), mm(DIN, p['Kn'])], axis=0)

def grazing_linear(phy, temp, p):
    return p['Cg20'] * _tempfunc(temp, p['theta_Cg']) * phy

def r_assim(phy, p):
    return p['as'] * mm(p['Ky'], phy)

def phy_respiration(temp, p):
    return p['Ry20'] * _tempfunc(temp, p['theta_Ry'])

def zoo_respiration(temp, p):
    return p['Rz20'] * _tempfunc(temp, p['theta_Rz'])

def zoo_mortal(p):
    return p['Dz']

def det_remineralization(temp, p):
    return p['Kd20'] * _tempfunc(temp, p['theta_Kd'])