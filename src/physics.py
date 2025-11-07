import numpy as np

def omegaValues(x, y, mu, collisionChecker):
    collided = False
    r1 = np.sqrt(((x + mu)**2) + y**2) # Distance to Pluto from particle
    r2 = np.sqrt(((x - 1 + mu)**2) + y**2) # Distance to Charon from particle
    if collisionChecker == True:
        if r1 < 1e-6:
            r1 = 1e-6
            collided = True
        if r2 < 1e-6:
            r2 = 1e-6
            collided = True
    # Scalar field strength at a point
    OmegaAtxy = ((1/2) * ((x**2) + (y**2))) + ((1 - mu) / r1) + (mu/r2)
    # Gradient of Omega at x
    OmegaX = x - ((1 - mu) * ((x + mu) / (r1**3))) - (mu * ((x - 1 + mu) / (r2**3)))
    # Gradient of Omega at y
    OmegaY = y - ((1 - mu) * (y / (r1**3))) - ((mu) * (y / (r2**3)))
    return [OmegaAtxy, OmegaX, OmegaY, collided] # [0] will be used for the Z Grid.

def equationsOfMotion(x, y, vx, vy, mu):
    dOmegaX, dOmegaY, collided = omegaValues(x, y, mu, True)[1], omegaValues(x, y, mu, True)[2], omegaValues(x, y, mu, True)[3]
    accelerationX = (2 * vy) + dOmegaX
    accelerationY = (-2 * vx) + dOmegaY
    return [vx, vy, accelerationX, accelerationY, collided]

def getJacobiC(Omega, vx, vy):  #Derivation available in paper
    C = (2 * Omega) - ((vx**2) + (vy ** 2))
    return C
