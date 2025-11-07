# ============= Zero-Velocity Curves ================
# Current configuration is on Figure 1, for other
# figures change the following:
# Figure 2a: 


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

gridLength = 300 #Grid goes from bottom left to bottom right to top right

#Arrays
ZGridValues = np.zeros((gridLength, gridLength))

xPointsArray = np.zeros((gridLength, gridLength))
yPointsArray = np.zeros((gridLength, gridLength))

plotPointsX = []
plotPointsY = []

#Lagrange
L1, L2, L3, L4, L5 = getLPoints(mu)

currentPoint = L4

# Building Z-Grid
Om = omegaValues(currentPoint[0], currentPoint[1], mu, True)[0]

currentC = getJacobiC(Om,0,0) + 1e-3

xPoints = np.linspace(-1.5, 1.5, gridLength)
yPoints = np.linspace(-1.5, 1.5, gridLength)

for i in range(gridLength): 
    for j in range(gridLength):
        ZGridValues[i][j] = 2 * omegaValues(xPoints[j], yPoints[i], mu, True)[0] #i AND j WERE SWAPPED
        xPointsArray[i][j] = xPoints[j]
        yPointsArray[i][j] = yPoints[i]

# Creates Array with trajectory

# Plotting
fig, ax = plt.subplots()

plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

plt.contour(xPoints, yPoints, ZGridValues, levels = [currentC]) #ALTER C
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
plt.xlabel("x-Axis")
plt.ylabel("y-Axis")
# plt.title(f"Zero-Velocity Curves with C = {currentC:.4f}")
plt.title(r"Chaotic trajectory near L3 with $\mu$ = " + f"{mu:.4f} (Pluto-Charon)")

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

