import numpy as np

def L1NRM(x, mu):
    newX = x - ((x - ((1 - mu) / ((x + mu)**2)) + (mu / ((x - 1 + mu)**2))) / (1 + ((2 * (1 - mu)) / (np.abs(x + mu)**3)) + ((2 * mu)/((np.abs(x - 1 + mu))**3))))
    return newX

def L2NRM(x, mu):
    newX = x - ((x - ((1 - mu) / ((x + mu)**2)) - (mu / ((x - 1 + mu)**2))) / (1 + ((2 * (1 - mu)) / (np.abs(x + mu)**3)) + ((2 * mu)/((np.abs(x - 1 + mu))**3))))
    return newX

def L3NRM(x, mu):
    newX = x - ((x + ((1 - mu) / ((x + mu)**2)) + (mu / ((x - 1 + mu)**2))) / (1 + ((2 * (1 - mu)) / (np.abs(x + mu)**3)) + ((2 * mu)/((np.abs(x - 1 + mu))**3))))
    return newX

def getLPoints(mu, tolerance=1e-8):
    converged = False
    XForL1 = (1/2) - mu
    XForL2 = 1 - mu + 0.3
    XForL3 = -mu - 0.3
    L4 = [(1/2) - mu, np.sqrt(3) / 2]
    L5 = [(1/2) - mu, -np.sqrt(3) / 2]

    while converged == False:   #LPoints are where the gradient of OmegaX = 0, zero movement.
        if np.abs(XForL1 - L1NRM(XForL1, mu)) < tolerance:   #Netwon-Raphson!
            converged = True
        XForL1 = L1NRM(XForL1, mu)
    converged = False

    while converged == False:
        if np.abs(XForL2 - L2NRM(XForL2, mu)) < 1e-8:
            converged = True
        XForL2 = L2NRM(XForL2, mu)
    converged = False

    while converged == False:
        if np.abs(XForL3 - L3NRM(XForL3, mu)) < 1e-8:
            converged = True
        XForL3 = L3NRM(XForL3, mu)

    return (XForL1, 0), (XForL2, 0), (XForL3, 0), L4, L5