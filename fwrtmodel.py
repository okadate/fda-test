# -*- coding: utf-8 -*-

import modelcomponents
import numpy as np
from tqdm import tqdm

from diffusion import *
from semi_lagrangian import *

i_irr = 0
i_tem = 1
iphy = [i_irr, i_tem]
nphy = len(iphy)

i_DIN = nphy + 0
i_DIP = nphy + 1
i_phy = nphy + 2
i_zoo = nphy + 3
i_PON = nphy + 4
i_POP = nphy + 5
ibio = [i_DIN, i_DIP, i_phy, i_zoo, i_PON, i_POP]
nbio = len(ibio)


class FWRTmodel(object):
    
    """
    Fresh Water Red Tide model
    
    """

    def __init__(self, delta_t=0.01, **parameters):
        
        self.lightattenuation = modelcomponents.lightattenuation
        
        self.phy_tempfunc = modelcomponents.phy_tempfunc
        self.lightresponse = modelcomponents.lightresponse
        self.nutrientuptake = modelcomponents.nutrientuptake_PN
        
        self.grazing = modelcomponents.grazing_linear
        self.r_assim = modelcomponents.r_assim
        
        self.phy_respiration = modelcomponents.phy_respiration
        
        self.zoo_respiration = modelcomponents.zoo_respiration
        self.zoo_mortal = modelcomponents.zoo_mortal
        
        self.det_remineralization = modelcomponents.det_remineralization
        
        # load standard parameters and update them with user-supplied values
        self.p = modelcomponents.standard_parameters.copy()
        self.p.update(parameters)
        
        self.delta_t = delta_t
        
        # the stepsize h for the FDA-based tangent linear and adjoint models
        self.h_tl = 1e-6
        self.h_ad = 1e-6

    #
    # the nonlinear (nl) code
    #
        
    def p_growth(self, x_nl, diag=None):
        #
        # a phytoplankton response to light
        # and nutrient uptake formulation
        #
        growth  = self.p['Gy_max20']
        growth *= self.phy_tempfunc(x_nl[i_tem], self.p)
        growth *= self.lightresponse(x_nl[i_irr], x_nl[i_phy], H=np.arange(0.5, 10, 1), p=self.p)
        growth *= self.nutrientuptake(DIP=x_nl[i_DIP], DIN=x_nl[i_DIN], p=self.p)
        growth *= x_nl[i_phy]
        
        if diag is not None:
            diag[0] = self.phy_tempfunc(x_nl[i_tem], self.p)
            diag[1] = self.lightresponse(x_nl[i_irr], x_nl[i_phy], H=0.5, p=self.p)
            diag[2] = self.nutrientuptake(DIP=x_nl[i_DIP], DIN=x_nl[i_DIN], p=self.p)
        
        # apply growth to nutrients and phytoplankton
        x_nl[i_DIP] -= growth*self.delta_t
        x_nl[i_DIN] -= growth*self.delta_t * self.p['phyNP']
        x_nl[i_phy] += growth*self.delta_t
        
        if diag is not None:
            return x_nl, diag
        else:
            return x_nl

    def z_grazing(self, x_nl):
        #
        # a zooplankton grazing on
        # phytoplankton formulation
        #
        grazing  = self.grazing(phy=x_nl[i_phy], temp=x_nl[i_tem], p=self.p)
        grazing *= x_nl[i_zoo]
        
        assim = self.r_assim(x_nl[i_phy], self.p)
        
        # apply grazing to phytoplankton and zooplankton
        x_nl[i_phy] -= grazing*self.delta_t
        x_nl[i_zoo] += grazing*self.delta_t * assim
        x_nl[i_POP] += grazing*self.delta_t * (1.0-assim)
        x_nl[i_PON] += grazing*self.delta_t * (1.0-assim) * self.p['phyNP']
        
        return x_nl
    
    def p_loss(self, x_nl):
        #
        # a phytoplankton loss formulation
        #
        
        loss  = self.phy_respiration(x_nl[i_tem], self.p)
        loss *= x_nl[i_phy]
        
        # apply loss to phytoplankton and nutrients
        x_nl[i_phy] -= loss*self.delta_t
        x_nl[i_POP] += loss*self.delta_t * (1.0-self.p['ry'])
        x_nl[i_PON] += loss*self.delta_t * (1.0-self.p['ry']) * self.p['phyNP']
        x_nl[i_DIP] += loss*self.delta_t * self.p['ry']
        x_nl[i_DIN] += loss*self.delta_t * self.p['ry'] * self.p['phyNP']
        
        return x_nl
    
    def z_loss(self, x_nl):
        #
        # a zooplankton loss formulation
        #
        loss  = self.zoo_respiration(x_nl[i_tem], self.p)
        loss *= x_nl[i_zoo]
        
        # apply loss to zooplankton and nutrients
        x_nl[i_zoo] -= loss*self.delta_t
        x_nl[i_POP] += loss*self.delta_t * (1.0-self.p['rz'])
        x_nl[i_PON] += loss*self.delta_t * (1.0-self.p['rz']) * self.p['zooNP']
        x_nl[i_DIP] += loss*self.delta_t * self.p['rz']
        x_nl[i_DIN] += loss*self.delta_t * self.p['rz'] * self.p['zooNP']
        
        #
        # a zooplankton loss formulation
        #
        loss  = self.zoo_mortal(self.p)
        loss *= x_nl[i_zoo]
        
        # apply loss to zooplankton and nutrients
        x_nl[i_zoo] -= loss*self.delta_t
        x_nl[i_POP] += loss*self.delta_t
        x_nl[i_PON] += loss*self.delta_t * self.p['zooNP']
        
        return x_nl
    
    def d_loss(self, x_nl):
        
        loss  = self.det_remineralization(x_nl[i_tem], self.p)
        loss *= x_nl[i_POP]
        
        x_nl[i_POP] -= loss*self.delta_t
        x_nl[i_DIP] += loss*self.delta_t
        
        loss  = self.det_remineralization(x_nl[i_tem], self.p)
        loss *= x_nl[i_PON]
        
        x_nl[i_PON] -= loss*self.delta_t
        x_nl[i_DIN] += loss*self.delta_t
        
        return x_nl
    
    def settling(self, H, x_nl):
        
        w = np.zeros(shape=(nphy+nbio))
        w[i_phy] = self.p['Wy']
        w[i_POP] = self.p['Wd']
        w[i_PON] = self.p['Wd']
        
        for i in [i_phy, i_POP, i_PON]:
            u = np.zeros_like(x_nl[i])
            u[:-1] = w[i]
            
            sl = SemiLagrangian1D(H, x_nl[i], u, self.delta_t)
            x_nl[i] = sl.linear_evolve(nt=1)
            #x_nl[i] = sl.cubic_evolve(nt=1)
            
        return x_nl
        
    def run_nl(self, ini, num_t, irr, temp, dH, H, diag=False):
        
        _, nx = ini.shape
        
        # create history
        x_nl_history = np.full(shape=(num_t+1,nphy+nbio,nx), fill_value=np.nan)
        x_nl_history[0,i_irr] = self.lightattenuation(irr[0], x_nl_history[0,i_phy], H)
        x_nl_history[0,i_tem] = temp[0]
        x_nl_history[0,2:] = ini
        
        # criate diagnostics
        if diag:
            diag = np.zeros(shape=(num_t+1, 3, nx))
            
        # settling speeds
        
        # current conditions
        x_nl = x_nl_history[0,:].copy()
        
        for t in tqdm(range(num_t)):
            
            for i in ibio:
                x_nl[i] = diffusion(x_nl[i], 1e-5, self.delta_t, dH, nx)
            
            # call each model segment in sequence
            #x_nl, diag[t] = self.p_growth(x_nl, diag[t])
            #self.z_grazing(x_nl)
            #self.p_loss(x_nl)
            #self.z_loss(x_nl)
            #self.d_loss(x_nl)
            self.settling(H, x_nl)
   
            # obtain light
            x_nl[i_irr] = self.lightattenuation(irr[t+1], x_nl[i_phy], H)
            
            x_nl[i_tem] = temp[t+1]
            
            # record history
            x_nl_history[t+1,:] = x_nl
            
        if diag is not False:
            return x_nl_history, diag
        else:
            return x_nl_history

    #
    # the FDA-based tangent linear (tl) code
    #
    
    # the FDA-based tangent linear phytoplankton growth function is written out 
    # explicitly, all others use the _generic_tl function below
    def p_growth_tl(self, x_nl, x_tl):
        x_pos = np.zeros(4)
        x_neg = np.zeros(4)
        
        # step (1): split the tangent linear 
        # model state and compute 
        # x_pos and x_neg
        for i in range(4):
            if (x_nl[i] + self.h_tl*x_tl[i]) >= 0:
                x_pos[i] = x_nl[i] + self.h_tl*x_tl[i]
                x_neg[i] = x_nl[i] 
            else:
                x_pos[i] = x_nl[i]
                x_neg[i] = x_nl[i] - self.h_tl*x_tl[i]
        
        # step (2): 2 evaluations of the 
        # nonlinear model
        x_pos = self.p_growth(x_pos)
        x_neg = self.p_growth(x_neg)
        
        # step (3): compute difference between 
        # the nonlinear model outputs and 
        # divide it by h (self.h_tl)
        for i in range(4):
            x_tl[i] = (x_pos[i]-x_neg[i])/self.h_tl
        
        return x_tl
    
    # A function providing the FDA-based tangent linear code for any 
    # nonlinear model segment function (generic_nl).
    # For usage see definitions of z_grazing_tl, p_loss_tl, etc. below.
    def _generic_tl(self, generic_nl, x_nl, x_tl):
        x_pos = np.zeros(4)
        x_neg = np.zeros(4)
        
        # step (1): split the tangent linear 
        # model state and compute 
        # x_pos and x_neg
        for i in range(4):
            if (x_nl[i] + self.h_tl*x_tl[i]) >= 0:
                x_pos[i] = x_nl[i] + self.h_tl*x_tl[i]
                x_neg[i] = x_nl[i] 
            else:
                x_pos[i] = x_nl[i]
                x_neg[i] = x_nl[i] - self.h_tl*x_tl[i]
        
        # step (2): 2 evaluations of the 
        # nonlinear model
        x_pos = generic_nl(x_pos)
        x_neg = generic_nl(x_neg)
        
        # step (3): compute difference between 
        # the nonlinear model outputs and 
        # divide it by h (self.h_tl)
        for i in range(4):
            x_tl[i] = (x_pos[i]-x_neg[i])/self.h_tl
        
        return x_tl
    
    def z_grazing_tl(self, x_nl, x_tl):
        return self._generic_tl(self.z_grazing, x_nl, x_tl)
        
    def p_loss_tl(self, x_nl, x_tl):
        return self._generic_tl(self.p_loss, x_nl, x_tl)
        
    def z_loss_tl(self, x_nl, x_tl):
        return self._generic_tl(self.z_loss, x_nl, x_tl)
        
    def run_tl(self, npz_ini, x_tl_ini, num_t):        
        '''Run the FDA-based tangent linear model.

Parameters
----------
npz_ini: array-like
    The initial conditions for nutrients, phytoplankton, and zooplankton.
x_tl_ini: array-like
    The initial conditions for the tangent linear model (including light, 
    nutrients, phytoplankton, and zooplankton).
num_t: int
    The number of time steps to run the model for.
    
Returns
-------
x_nl_history: array
    An `num_t`-by-4 array containing the nonlinear model state at each time 
    step.
x_tl_history: array
    An `num_t`-by-4 array containing the tangent linear model state at each 
    time step.
        '''
        
        # generate light
        light = self.generate_light(num_t+1)
        
        # create histories
        x_nl_history = np.full(shape=(num_t+1,4), fill_value=np.nan)
        x_nl_history[0,1:] = npz_ini
        x_nl_history[0,i_irr] = light[0]
        
        x_tl_history = np.full(shape=(num_t+1,4), fill_value=np.nan)
        x_tl_history[0,:] = x_tl_ini
        
        # current conditions
        x_nl = x_nl_history[0,:].copy()
        x_tl = x_tl_history[0,:].copy()
        
        for t in range(num_t):
            # for each of the 4 segments: 
            # step (1) update the tangent linear 
            # state x_tl using the correct 
            # value of the nonlinear state x_nl
            # step (2) propagate the nonlinear 
            # state forward
            
            # segment: P growth
            x_tl = self.p_growth_tl(x_nl, x_tl)
            x_nl = self.p_growth(x_nl)
            
            # segment: Z grazing
            x_tl = self.z_grazing_tl(x_nl, x_tl)
            x_nl = self.z_grazing(x_nl)
            
            # segment: P loss
            x_tl = self.p_loss_tl(x_nl, x_tl)
            x_nl = self.p_loss(x_nl)
            
            # segment: Z loss
            x_tl = self.z_loss_tl(x_nl, x_tl)
            x_nl = self.z_loss(x_nl)
            
            # obtain light
            x_nl[i_irr] = light[t+1]
            
            # record history
            x_nl_history[t+1,:] = x_nl
            x_tl_history[t+1,:] = x_tl
            
        return x_nl_history, x_tl_history
    
    #
    # the FDA-based adjoint (ad) code
    #
    
    # A function providing the FDA-based adjoint code for any 
    # nonlinear model segment function (generic_nl).
    # For usage see definitions of p_growth_ad, z_grazing_ad etc. below.
    def _generic_ad(self, generic_nl, ivars_in, x_nl, x_ad, x_nl_ref):
        
        # Note that nonlinear reference 
        # solution for the unperturbed 
        # nonlinear state is provided as 
        # an input argument here but could 
        # be computed instead.
        
        # create copy of adjoint state
        x_ad_orig = x_ad.copy()
        
        for i in ivars_in:
            # create a copy of x_nl 
            # with variable i perturbed by h (self.h_ad)
            x_nl_h = x_nl.copy()
            x_nl_h[i] = x_nl_h[i] + self.h_ad
            
            # nonlinear model evaluation for 
            # perturbed state
            x_nl_h = generic_nl(x_nl_h)
            
            # compute matrix product and save 
            # result in delta_ad
            delta_ad = 0.0
            for j in range(4):
                delta_ad = delta_ad + (x_nl_h[j] - x_nl_ref[j])*x_ad_orig[j]
                
            x_ad[i] = delta_ad/self.h_ad
        return x_ad
    
    # for p_growth: light, N, and P act as input
    def p_growth_ad(self, x_nl, x_ad, x_nl_ref):
        return self._generic_ad(self.p_growth, (i_irr,i_nut,i_phy), x_nl, x_ad, x_nl_ref)
    
    # for z_grazing: P and Z act as input
    def z_grazing_ad(self, x_nl, x_ad, x_nl_ref):
        return self._generic_ad(self.z_grazing, (i_phy,i_zoo), x_nl, x_ad, x_nl_ref)
    
    # for p_loss: P acts as input
    def p_loss_ad(self, x_nl, x_ad, x_nl_ref):
        return self._generic_ad(self.p_loss, (i_phy,), x_nl, x_ad, x_nl_ref)
    
    # for z_loss: Z acts as input
    def z_loss_ad(self, x_nl, x_ad, x_nl_ref):
        return self._generic_ad(self.z_loss, (i_zoo,), x_nl, x_ad, x_nl_ref)
    
    def run_ad(self, npz_ini, x_ad_ini, num_t=None, x_nl_history=None):
        '''Run the FDA-based adjoint model.

Parameters
----------
npz_ini: array-like
    The initial conditions for nutrients, phytoplankton, and zooplankton.
x_ad_ini: array-like
    The initial conditions for the adjoint model (including light, 
    nutrients, phytoplankton, and zooplankton).
num_t: int, optional
    The number of time steps to run the model for (optional if `x_nl_history`
    is supplied).
x_nl_history: array
    An `num_t`-by-4 array containing the nonlinear model state at each time 
    step. If `x_nl_history` is provided, it is used as the nonlinear model 
    state in the computation of the adjoint. Otherwise `x_nl_history` is 
    computed using `npz_ini` and `num_t`.
    
Returns
-------
x_nl_history: array
    An `num_t`-by-4 array containing the nonlinear model state at each time 
    step.
x_ad_history: array
    An `num_t`-by-4 array containing the adjoint model state at each time 
    step.
        '''
        
        # run nonlinear 
        if x_nl_history is None:
            x_nl_history = self.run_nl(npz_ini=npz_ini, num_t=num_t)
        
        #print(x_nl_history)
        
        x_ad_history = np.full(shape=(num_t+1,4), fill_value=np.nan)
        x_ad_history[-1,:] = x_ad_ini
        
        x_ad = x_ad_history[-1,:].copy()
        
        for t in reversed(range(num_t)):
            # load appropriate nonlinear state
            x_nl = x_nl_history[t,:].copy()
            
            # First, propagate the nonlinear model 
            # solution and save it at every stage.
            
            # segment: P growth
            x_nl_p_growth = x_nl.copy()
            x_nl = self.p_growth(x_nl)
            
            # segment: Z grazing
            x_nl_z_grazing = x_nl.copy()
            x_nl = self.z_grazing(x_nl)
            
            # segment: P loss
            x_nl_p_loss = x_nl.copy()
            x_nl = self.p_loss(x_nl)
            
            # segment: Z loss
            x_nl_z_loss = x_nl.copy()
            x_nl = self.z_loss(x_nl)
            
            # Second, propagate adjoint model 
            # solution in reverse order using 
            # the appropriate saved nonlinear 
            # model solution.
            # Use saved nonlinear model solution 
            # for the next segment as the 
            # reference solution for the current 
            # segment (3rd input argument).
            
            # segment: Z loss
            x_ad = self.z_loss_ad(x_nl_z_loss, x_ad, x_nl)
            
            # segment: P loss
            x_ad = self.p_loss_ad(x_nl_p_loss, x_ad, x_nl_z_loss)
            
            # segment: Z grazing
            x_ad = self.z_grazing_ad(x_nl_z_grazing, x_ad, x_nl_p_loss)
            
            # segment: P growth
            x_ad = self.p_growth_ad(x_nl_p_growth, x_ad, x_nl_z_grazing)
            
            # record history
            x_ad_history[t,:] = x_ad
        return x_nl_history, x_ad_history
    
if __name__=='__main__':
    m = NPZModel()
    
    x_nl = m.run_nl(npz_ini=(10.,5.,5.), num_t=500)
    
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    ax.plot(x_nl)
    ax.legend(('I','N','P','Z'))
    ax.grid(True)
    ax.set(title='nonlinear model simulation', xlabel='time step', 
           ylabel='variable value')
    
    # plot tl
    
    x_nl, x_tl = m.run_tl(npz_ini=(10.,5.,5.), x_tl_ini=(1.0,1.0,0.2,0.2), num_t=500)
        
    colors = ('C0','C1','C2','C3')
    
    fig, ax = plt.subplots()
    
    for i,name in enumerate('INPZ'):
        ax.plot(x_nl[:,i], color=colors[i], label=name)
        ax.plot(x_nl[:,i]+x_tl[:,i], color=colors[i], alpha=0.5)
        
    ax.legend()
    ax.grid(True)
    ax.set(title='nonlinear and tangent linear model simulation', 
           xlabel='time step', ylabel='variable value')
    
    plt.show()
    
    
    
    
