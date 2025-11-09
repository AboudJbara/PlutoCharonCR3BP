# ============= Chaotic Transit at Neck Figures ================
# Current configuration is on ReflectionAtL1, for TransitAtL1, change
# the following:
#
# - TransitAtL1 
#       - line 48: Change 1e-2 -> -(1e-2)


import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from src.physics import *
from src.lpoints import *
from src.integrator import *

# Parameters
massPluto = 1.3025e22 #Uncertainty = +-0.0006
massCharon = 1.5897e21 #Uncertainty = +-0.0045
mu = massCharon/(massCharon+massPluto)

time = (5e-4)
T = 50     # Increase this value for more librations.
repetitions = int(T/time)
gridLength = 300 # Grid goes from bottom left to bottom right to top right

#Arrays
ZGridValues = np.zeros((gridLength, gridLength))

xPointsArray = np.zeros((gridLength, gridLength))
yPointsArray = np.zeros((gridLength, gridLength))

plotPointsX = []
plotPointsY = []

#Lagrange
L1, L2, L3, L4, L5 = getLPoints(mu)

#               #
currentPoint = L1
#               #

CAtLagrangePoint = getJacobiC(omegaValues(currentPoint[0], currentPoint[1], mu, True)[0], 0, 0)

# Initial conditions (Vertical Velocity Case)
TargetC = CAtLagrangePoint + 1e-2   # NOTE: This is the only value that will need to be changed

x = currentPoint[0] - 0.04
y = currentPoint[1]

velocityAbs = np.sqrt((2 * omegaValues(x, y, mu, True)[0]) - TargetC)

vx = 0            
vy = velocityAbs      

# Building Z-Grid
Om = omegaValues(x,y, mu, True)[0]
currentC = getJacobiC(Om,vx,vy)

xPoints = np.linspace(-1.5, 1.5, gridLength)
yPoints = np.linspace(-1.5, 1.5, gridLength)

for i in range(gridLength): 
    for j in range(gridLength):
        ZGridValues[i][j] = 2 * omegaValues(xPoints[j], yPoints[i], mu, True)[0]  # Loop produces Z-Grid
        xPointsArray[i][j] = xPoints[j]
        yPointsArray[i][j] = yPoints[i]

# Creates an Array with the trajectories
trajectoryX, trajectoryY, maxC, collided = propagate(x, y, vx, vy, time, repetitions, currentC, mu)

#Checks collisions
print("Drift in C: " + str(maxC)) # Note: if value here is too high, then RK4 has failed; results are innacurate.
if collided == True:
    print("OBJECT HAS COLLIDED.")


# Plotting
fig, ax = plt.subplots()

plt.plot(trajectoryX[:repetitions], trajectoryY[:repetitions], "-", linewidth = 0.6)
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
plt.xlabel("x-Axis")
plt.ylabel("y-Axis")
plt.title(r"Transit/Reflection at L1 Neck")

SPACECRAFT = Circle((currentPoint[0] + 1e-2, currentPoint[1]), 0.02, facecolor = "red")


pluto = Circle((-mu, 0), 0.1, facecolor = "#778899", edgecolor = "#F08080")
charon = Circle((1 - mu, 0), 0.05, facecolor = "#C0C0C0", edgecolor = "#D2BB9D")
L1Patch = Circle((L1[0], L1[1]), 0.01, facecolor = "black")
L2Patch = Circle((L2[0], L2[1]), 0.01, facecolor = "black")
L3Patch = Circle((L3[0], L3[1]), 0.01, facecolor = "black")
L4Patch = Circle((L4[0], L4[1]), 0.01, facecolor = "black")
L5Patch = Circle((L5[0], L5[1]), 0.01, facecolor = "black")
ax.add_patch(L1Patch)
ax.add_patch(L2Patch)
ax.add_patch(L3Patch)
ax.add_patch(L4Patch)
ax.add_patch(L5Patch)
plt.text(L1[0], L1[1], "L1", fontsize = 7)
plt.text(L2[0], L2[1], "L2", fontsize = 7)
plt.text(L3[0], L3[1], "L3", fontsize = 7)
plt.text(L4[0], L4[1], "L4", fontsize = 7)
plt.text(L5[0], L5[1], "L5", fontsize = 7)

ax.add_patch(pluto)
ax.add_patch(charon)
ax.set_aspect(1)


plt.show()
