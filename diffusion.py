import numpy as np

from TDMAsolver import TDMAsolver


def diffusion(x, coef, dt, dx, nx):
    
    # CFL condition
    k = coef * dt / (dx * dx)
    assert (k < 0.5).all(), k
    cff = k
    
    a = -0.5 * cff[1:]
    b = (1.0 + cff)
    c = -0.5 * cff[:-1]
    
    b[0]  = (1.0 + 0.5*cff[0])  # boundary conditions: ddC/dxdx=0 at boundaris
    b[-1] = (1.0 + 0.5*cff[-1])
    
    d = (1.0 - cff) * x
    d[0]  = (1.0 - 0.5 * cff[0]) * x[0]
    d[-1] = (1.0 - 0.5 * cff[-1]) * x[-1]
    
    d[1:]  += 0.5 * cff[:-1] * x[:-1]
    d[:-1] += 0.5 * cff[1:] * x[1:]
    
    return TDMAsolver(a, b, c, d)


if __name__ == '__main__':
    nt = 5000
    nx = 20

    diff = 1e-3
    dt = 0.1
    dx = 0.1

    print diff*dt/dx < 1., diff*dt/dx

    x_his = np.zeros(shape=(nt+1,nx))
    x_his[0] = np.arange(nx)

    x = x_his[0].copy()

    for t in range(nt):

        x = diffusion(x, diff, dt, dx, nx)

        x_his[t+1,:] = x

    import matplotlib.pyplot as plt
    
    plt.pcolor(x_his.T)