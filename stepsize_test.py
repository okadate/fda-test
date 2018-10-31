'''
Check sensitivity of the FDA technique to the value of the stepsize for the
tangent linear and adjoint model (h_tl and h_ad, respectively). Sensitivity
is assessed by performing num_tests dot-product tests for random initial 
conditions and tangent linear and adjoint state vectors for every combination
of h_tl and h_ad.
'''

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from NPZModel import NPZModel

m = NPZModel()
num_t = 10

num_tests = 100

npz_inis = np.random.uniform(size=(num_tests,3), low=1.0, high=10.0)
x_tls = np.random.uniform(size=(num_tests,4), low=-1.0, high=1.0)
x_ads = np.random.uniform(size=(num_tests,4), low=-1.0, high=1.0)

# values for h_tl and h_ad 
h_tls = np.array([1e-1, 1e-3, 1e-5, 1e-7, 1e-9, 1e-11]) 
h_ads = np.array([1e-1, 1e-3, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-11]) 

# arrays to store error values
mean_errors = np.full(shape=(h_tls.size,h_ads.size), fill_value=np.nan)
errors = np.empty(shape=num_tests)

for itl, h_tl in enumerate(h_tls):
    for iad, h_ad in enumerate(h_ads):
        
        # set model stepsize
        m.h_lt = h_tl
        m.h_ad = h_ad
        
        # clean old results
        errors[:] = np.nan
        
        print('performing {} tests for h_tl={:.1e}, h_ad={:.1e}'.format(num_tests, h_tl, h_ad), end='')
        for itest in range(num_tests):
            # load initial conditions and tangent linear and adjoint state vectors
            npz_ini = npz_inis[itest,:]
            x_tl = x_tls[itest,:]
            x_ad = x_ads[itest,:]
            
            res_tl = m.run_tl(npz_ini=npz_ini, x_tl_ini=x_tl, num_t=num_t)[1][-1,:]
            res_ad = m.run_ad(npz_ini=npz_ini, x_ad_ini=x_ad, num_t=num_t)[1][0,:]
            
            errors[itest] = np.abs(x_ad.dot(res_tl) - x_tl.dot(res_ad))
            
        mean_errors[itl,iad] = np.mean(errors)
        print('; mean error: {:.3e}'.format(mean_errors[itl,iad]))
        
# print results

print('\nsummary (num_t={}):'.format(num_t))
print('{:7s} | {}'.format('', 'h_ad'))
print('{:7s} | {}'.format('h_tl', ' '.join(['{:7.1e}'.format(x) for x in h_ads])))
print((2+8*(len(h_ads)+1))*'-')
for irow,row in enumerate(mean_errors):
    print('{:7.1e} | {}'.format(h_tls[irow], ' '.join(['{:7.1e}'.format(x) for x in row])))


