import numpy as np

from .physics import equationsOfMotion, omegaValues, getJacobiC


def positionChanger(startX, startY, startVX, startVY, time, mu): #Runge-Kutta Fourth Order Method
    #Returns a tuple with (k1vx, k1vy, k1ax, k1ay, collided)
    k1 = equationsOfMotion(startX, startY, startVX, startVY, mu)
    k2 = equationsOfMotion(startX + ((time/2) * k1[0]), startY + ((time/2) * k1[1]), startVX + ((time/2) * k1[2]), startVY + ((time/2) * k1[3]), mu)
    k3 = equationsOfMotion(startX + ((time/2) * k2[0]), startY + ((time/2) * k2[1]), startVX + ((time/2) * k2[2]), startVY + ((time/2) * k2[3]), mu)
    k4 = equationsOfMotion(startX + (time * k3[0]), startY + (time * k3[1]), startVX + (time * k3[2]), startVY + (time * k3[3]), mu)
    posX = startX + ((time/6) * (k1[0] + (2*k2[0]) + (2*k3[0]) + k4[0]))
    posY = startY + ((time/6) * (k1[1] + (2*k2[1]) + (2*k3[1]) + k4[1]))
    vx = startVX + ((time/6) * (k1[2] + (2*k2[2]) + (2*k3[2]) + k4[2]))
    vy = startVY + ((time/6) * (k1[3] + (2*k2[3]) + (2*k3[3]) + k4[3]))

    return posX, posY, vx, vy, k1[4]


def propagate(x, y, vx, vy, time, steps, startingC, mu):   #positionChanger returns x,y,vx,vy,collision
    trajXArray = np.zeros(steps)
    trajYArray = np.zeros(steps)
    jacobiDriftTest = []
    for i in range (steps):
        x, y, vx, vy, collided = positionChanger(x, y, vx, vy, time, mu)
        trajXArray[i], trajYArray[i] = x, y
        #DIFFERENCE IN C TEST
        C = getJacobiC(omegaValues(x, y, mu, True)[0], vx, vy)
        jacobiDriftTest.append(np.abs(C - startingC))
        
        if collided == True:
            print("BREAKING")
            break
    jacobiMax = max(jacobiDriftTest)
    return trajXArray, trajYArray, jacobiMax, collided

    
