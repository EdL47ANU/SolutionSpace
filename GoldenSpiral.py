import numpy as np
import matplotlib.pyplot as plt
#This file shows the difference btween the old and new Theta and Phi generation methods. 
# It used to favour the Z axis and negate a perfect negatve z.
#There are as many blue points as red, but the red points are more evenly distributed.

NOrient = 10
n = NOrient**2
i = np.arange(0, n, dtype=float) + 0.5
phi = np.arccos(1 - 2*i/n)
goldenRatio = (1 + 5**0.5)/2
theta = 2 * np.pi * i / goldenRatio
x, y, z = np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)

#OLd Scatter plot
OldTheta = np.arange(0, 2*np.pi, 2*np.pi/NOrient)
OldPhi = np.arange(0, np.pi, np.pi/NOrient)

fig1 = plt.figure()
fig2 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')
ax2 = fig2.add_subplot(111, projection='3d')
ax.scatter(x, y, z, 'o', color='r', label = 'Golden Spiral')
LegendFixer = True
for i in range(len(OldTheta)):
    #ax2.plot(np.cos(OldTheta[i])*np.sin(OldPhi), np.sin(OldTheta[i])*np.sin(OldPhi), np.cos(OldPhi), color='b')
    for j in range(len(OldPhi)):
        if LegendFixer:
            ax2.scatter(np.cos(OldTheta[i])*np.sin(OldPhi[j]), np.sin(OldTheta[i])*np.sin(OldPhi[j]), np.cos(OldPhi[j]), 'o', color='b', label = 'Old Spiral')
            LegendFixer = False
        else:
            ax2.scatter(np.cos(OldTheta[i])*np.sin(OldPhi[j]), np.sin(OldTheta[i])*np.sin(OldPhi[j]), np.cos(OldPhi[j]), 'o', color='b')

plt.legend()
plt.show()

