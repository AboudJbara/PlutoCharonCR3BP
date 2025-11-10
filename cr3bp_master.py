#============================================================
# NOTE: This is a standalone, master version. No functions from other
# files are needed.  Current version takes ~10 seconds to run.
# Requires numpy and matplotlib
# 
# This version produces the textbook horseshoe shape, to get the 
# chaotic attempt at a horseshoe under the Pluto-Charon system, 
# set USE_TEXTBOOK_MU = False.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# === System Selection ===
USE_TEXTBOOK_MU = True  # True produces a clean horseshoe. False produces a chaotic horseshoe under Pluto-Charon.
# =========================

massPluto = 1.3025e22 #Uncertainty = ±0.0006
massCharon = 1.5897e21 #Uncertainty = ±0.0045

if USE_TEXTBOOK_MU:
    mu = 0.001
    print("Mode: textbook μ=0.001")
else:
    mu = massCharon/(massCharon+massPluto)
    print("Mode: Pluto-Charon μ")



#FUNCTIONS

#Everything will be simplified, separation= 1, angular rate n = 1, total mass = 1

def omegaValues(x,y,collisionChecker):
    collided = False
    r1 = np.sqrt(((x + mu)**2) + y**2) #distance to pluto from particle
    r2 = np.sqrt(((x - 1 + mu)**2) + y**2) #distance to charon from particle
    if collisionChecker == True:
        if r1 < 1e-6:
            r1 = 1e-6
            collided = True
        if r2 < 1e-6:
            r2 = 1e-6
            collided = True
    #Scalar field strength at a point
    OmegaAtxy = ((1/2) * ((x**2) + (y**2))) + ((1 - mu) / r1) + (mu/r2)
    #Gradient of Omega at x
    OmegaX = x - ((1 - mu) * ((x + mu) / (r1**3))) - (mu * ((x - 1 + mu) / (r2**3)))
    #Gradient of Omega at y
    OmegaY = y - ((1 - mu) * (y / (r1**3))) - ((mu) * (y / (r2**3)))
    return [OmegaAtxy, OmegaX, OmegaY, collided] #[0] is used for the Z Grid.

def equationsOfMotion(x, y, vx, vy):
    omegaX, omegaY, collision = omegaValues(x, y, True)[1], omegaValues(x, y, True)[2], omegaValues(x, y, True)[3]
    accelerationX = (2 * vy) + omegaX
    accelerationY = (-2 * vx) + omegaY
    return (vx, vy, accelerationX, accelerationY, collision)

def getJacobiC(Omega, xDot, yDot):  #Derived formula in paper.
    C = (2 * Omega) - ((xDot**2) + (yDot ** 2))
    return C

# Three different SIGN CHANGES because of the way absolute values work,
# whether the numerator is positive or negative depends entirely on
# the interval!

def L1NRM(x):
    newX = x - ((x - ((1 - mu) / ((x + mu)**2)) + (mu / ((x - 1 + mu)**2))) / (1 + ((2 * (1 - mu)) / (np.abs(x + mu)**3)) + ((2 * mu)/((np.abs(x - 1 + mu))**3))))
    return newX

def L2NRM(x):
    newX = x - ((x - ((1 - mu) / ((x + mu)**2)) - (mu / ((x - 1 + mu)**2))) / (1 + ((2 * (1 - mu)) / (np.abs(x + mu)**3)) + ((2 * mu)/((np.abs(x - 1 + mu))**3))))
    return newX

def L3NRM(x):
    newX = x - ((x + ((1 - mu) / ((x + mu)**2)) + (mu / ((x - 1 + mu)**2))) / (1 + ((2 * (1 - mu)) / (np.abs(x + mu)**3)) + ((2 * mu)/((np.abs(x - 1 + mu))**3))))
    return newX

def getLPoints():
    converged = False
    XForL1 = (1/2) - mu
    XForL2 = 1 - mu + 0.3
    XForL3 = -mu - 0.3
    L4 = [(1/2) - mu, np.sqrt(3) / 2]
    L5 = [(1/2) - mu, -np.sqrt(3) / 2]

    while converged == False:   # LPoints are where the gradient of OmegaX = 0, zero movement.
        if np.abs(XForL1 - L1NRM(XForL1)) < 1e-8:   # Newton-Raphson!
            converged = True
        XForL1 = L1NRM(XForL1)
    converged = False

    while converged == False:
        if np.abs(XForL2 - L2NRM(XForL2)) < 1e-8:
            converged = True
        XForL2 = L2NRM(XForL2)
    converged = False

    while converged == False:
        if np.abs(XForL3 - L3NRM(XForL3)) < 1e-8:
            converged = True
        XForL3 = L3NRM(XForL3)

    return (XForL1, 0), (XForL2, 0), (XForL3, 0), L4, L5

def positionChanger(startX, startY, startVX, startVY, time):
    # equationsOfMotion returns k1vx, k1vy, k1ax, k1ay, collided
    k1 = equationsOfMotion(startX, startY, startVX, startVY)
    k2 = equationsOfMotion(startX + ((time/2) * k1[0]), startY + ((time/2) * k1[1]), startVX + ((time/2) * k1[2]), startVY + ((time/2) * k1[3]))
    k3 = equationsOfMotion(startX + ((time/2) * k2[0]), startY + ((time/2) * k2[1]), startVX + ((time/2) * k2[2]), startVY + ((time/2) * k2[3]))
    k4 = equationsOfMotion(startX + (time * k3[0]), startY + (time * k3[1]), startVX + (time * k3[2]), startVY + (time * k3[3]))
    posX = startX + ((time/6) * (k1[0] + (2*k2[0]) + (2*k3[0]) + k4[0]))
    posY = startY + ((time/6) * (k1[1] + (2*k2[1]) + (2*k3[1]) + k4[1]))
    vx = startVX + ((time/6) * (k1[2] + (2*k2[2]) + (2*k3[2]) + k4[2]))
    vy = startVY + ((time/6) * (k1[3] + (2*k2[3]) + (2*k3[3]) + k4[3]))

    return posX, posY, vx, vy, k1[4]

# MAIN CODE
h = (5e-3)
T = 300     # Increase T for more librations/full horseshoe shape.
repetitions = int(T/h)
gridLength = 300 # Grid goes from bottom left to bottom right to top right

ZGridValues = np.zeros((gridLength, gridLength))

xPointsArray = np.zeros((gridLength, gridLength))
yPointsArray = np.zeros((gridLength, gridLength))

plotPointsX = []
plotPointsY = []

L1, L2, L3, L4, L5 = getLPoints()

#               #
currentPoint = L3
#               #

CAtLagrangePoint = getJacobiC(omegaValues(currentPoint[0], currentPoint[1], True)[0], 0, 0)

# Initial conditions (Vertical Velocity Case)
TargetC = CAtLagrangePoint - 5e-5

x = currentPoint[0] - 0.027
y = currentPoint[1] 

# #VELOCITY CALCULATIONS#
velocityAbs = np.sqrt((2 * omegaValues(x,y,True)[0]) - TargetC)

# TANGENTIAL
# vDirX, vDirY = -y, x #CLOCKWISE
# vDirX, vDirY = vDirX / (np.sqrt((vDirX**2) + (vDirY**2))), vDirY / (np.sqrt((vDirX**2) + (vDirY**2)))
# vx, vy = velocityAbs * vDirX, velocityAbs * vDirY

vx = 0 
vy = velocityAbs

Om = omegaValues(x,y, True)[0]
currentC = getJacobiC(Om,vx,vy)  #Gets C, using above functions

xPoints = np.linspace(-1.5, 1.5, gridLength)
yPoints = np.linspace(-1.5, 1.5, gridLength)

for i in range(gridLength): #ARRAYS
    for j in range(gridLength):
        ZGridValues[i][j] = 2 * omegaValues(xPoints[j], yPoints[i], True)[0] #i AND j WERE SWAPPED
        xPointsArray[i][j] = xPoints[j]
        yPointsArray[i][j] = yPoints[i]


trajXArray = np.zeros(repetitions)
trajYArray = np.zeros(repetitions)
count = 0
tempList = []
while count < repetitions: 
    # positionChanger returns x,y,vx,vy,collision
    x, y, vx, vy, collided = positionChanger(x, y, vx, vy, h)

    trajXArray[count], trajYArray[count] = x, y
    count += 1 
    # DIFFERENCE IN C TEST
    C = getJacobiC(omegaValues(x, y, True)[0], vx, vy)
    tempList.append(np.abs(C - currentC))

    if collided == True:
        print("BREAKING")
        break
print("Max drift in C: " + str(max(tempList)))


#PLOT
fig, ax = plt.subplots()

plt.plot(trajXArray[:count], trajYArray[:count], "-", linewidth = 0.6)
# plt.contour(xPoints, yPoints, ZGridValues, levels = [currentC]) # Add back for Zero-Velocity Curves
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
plt.xlabel("x-Axis")
plt.ylabel("y-Axis")
plt.title(r"Traditional Horseshoe Orbit with $\mu=0.001$")

SPACECRAFT = Circle((currentPoint[0] + 1e-2, currentPoint[1]), 0.02, facecolor = "red")


pluto = Circle((-mu, 0), 0.1, facecolor = "#000000", edgecolor = "#000000")
charon = Circle((1 - mu, 0), 0.05, facecolor = "#000000", edgecolor = "#000000")
# L1Patch = Circle((L1[0], L1[1]), 0.01, facecolor = "black")   # L1 and L2 are unnecessary to
# L2Patch = Circle((L2[0], L2[1]), 0.01, facecolor = "black")   # show, since they are on the
L3Patch = Circle((L3[0], L3[1]), 0.01, facecolor = "black")     # smaller primary.
L4Patch = Circle((L4[0], L4[1]), 0.01, facecolor = "black")
L5Patch = Circle((L5[0], L5[1]), 0.01, facecolor = "black")
# ax.add_patch(L1Patch)
# ax.add_patch(L2Patch)
ax.add_patch(L3Patch)
ax.add_patch(L4Patch)
ax.add_patch(L5Patch)
# plt.text(L1[0], L1[1], "L1", fontsize = 7)
# plt.text(L2[0], L2[1], "L2", fontsize = 7)
plt.text(L3[0], L3[1], "L3", fontsize = 7)
plt.text(L4[0], L4[1], "L4", fontsize = 7)
plt.text(L5[0], L5[1], "L5", fontsize = 7)

ax.add_patch(pluto)
ax.add_patch(charon)
ax.set_aspect(1)


plt.show()
