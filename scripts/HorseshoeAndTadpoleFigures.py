# ============= Tadpole and Horseshoe Figures ================
# Current configuration is on Figure 3, for Figure 4, change
# the following:
#
# - Figure 4   NOTE: To see full horseshoe, increase T (line 26) 
#       - line 43: Change L4 -> L3
#       - line 49: Change 5e-4 -> 5e-5
#       - line 51: Change -0.04 -> -0.027
#       - line 52: Change -0.03 -> 0
#       - REPLACE lines 57-59 with lines 62-63

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from src.physics import *
from src.lpoints import *
from src.integrator import *

# Parameters
mu = 0.001

time = (5e-4)
T = 150     # Increase this value for more librations.
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
currentPoint = L4
#               #

CAtLagrangePoint = getJacobiC(omegaValues(currentPoint[0], currentPoint[1], mu, True)[0], 0, 0)

# Initial conditions (Vertical Velocity Case)
TargetC = CAtLagrangePoint - 5e-4

x = currentPoint[0] - 0.04
y = currentPoint[1] - 0.03

velocityAbs = np.sqrt((2 * omegaValues(x, y, mu, True)[0]) - TargetC)

# Tangential Velocity   NOTE: Remove this for Figure 4, replace with lines 61 & 62
vDirX, vDirY = -y, x # Clockwise Configuration  
vDirX, vDirY = vDirX / (np.sqrt((vDirX**2) + (vDirY**2))), vDirY / (np.sqrt((vDirX**2) + (vDirY**2)))
vx, vy = velocityAbs * vDirX, velocityAbs * vDirY

# Manual Velocity NOTE: Replace tangential VELOCITY with this for Figure 4. 
# vx = 0            
# vy = velocityAbs      

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

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

plt.plot(trajectoryX[:repetitions], trajectoryY[:repetitions], "-", linewidth = 0.6)
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
plt.xlabel("x-Axis")
plt.ylabel("y-Axis")
plt.title(r"Horseshoe/Tadpole Orbits")

SPACECRAFT = Circle((currentPoint[0] + 1e-2, currentPoint[1]), 0.02, facecolor = "red")


pluto = Circle((-mu, 0), 0.1, facecolor = "#000000", edgecolor = "#000000")
charon = Circle((1 - mu, 0), 0.05, facecolor = "#000000", edgecolor = "#000000")
L3Patch = Circle((L3[0], L3[1]), 0.01, facecolor = "black")
L4Patch = Circle((L4[0], L4[1]), 0.01, facecolor = "black")
L5Patch = Circle((L5[0], L5[1]), 0.01, facecolor = "black")
ax.add_patch(L3Patch)
ax.add_patch(L4Patch)
ax.add_patch(L5Patch)
plt.text(L3[0], L3[1], "L3", fontsize = 7)
plt.text(L4[0], L4[1], "L4", fontsize = 7)
plt.text(L5[0], L5[1], "L5", fontsize = 7)

ax.add_patch(pluto)
ax.add_patch(charon)
ax.set_aspect(1)


plt.show()
