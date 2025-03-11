from HyperfineHamilClass import HyperH0
from Consts import Planck, MHzmT, SpinOp4, SpinOp2, SpinOp3, OrthoVecs, rotate_vectors
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def RotatingB1Data(axis, H, NOrient = 100):
    Z = axis/np.linalg.norm(axis)
    X,Y = OrthoVecs(axis)
    angles = np.arange(0, 2*np.pi, 2*np.pi/NOrient)
    New_x, New_y = [],[]
    bools = []
    for angle in angles:
        new_x, new_y = rotate_vectors(Z, X, Y, angle)
        New_x.append(new_x)
        New_y.append(new_y)
    for i in range(len(New_x)):
        H.FibbonachiRotateHamil(Z, New_x[i], New_y[i])
        H.MkEigensystem()
        H.MkTransFrq()
        H.MkTransInt()
        #print(np.round(H.TransInt,3))
        if H.FrqCheckIndependantPulses() and H.TransCheck() == True:
            bools.append(True)
        else:
            bools.append(False)
    print(np.sum(bools)/len(bools))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(New_x)):
        if bools[i]:
            ax.scatter(New_x[i][0], New_x[i][1], New_x[i][2], c = 'r')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.set_title('Solution Space for different orientations of the perturbing microwave field at a known successful B value to show the neccesity of chi.')
    plt.show()

def HyperSearch(Ax, Ay, Az, B, Orientations, gi):
    AxCoords = Ax
    AyCoords = Ay
    AzCoords = Az
    AxboolsTot = []
    #Thetas, Phis, Chis = GoldenSpiralAnglesRoll(NOrient)
    for ax in Ax:
        AyboolsTot = []
        for ay in Ay:
            AzboolsTot = []
            for az in Az:
                Bools = []
                for b in B:
                    H = HyperH0(A = np.diag([ax, ay, az]), coil = b, SpinOp = SpinOp2, NucSpinOp=SpinOp2, gi=np.diag([gi,gi,gi]), gs=np.diag([2.00232,2.00232,2.00232]))
                    bools = HyperSolChecker1Bool(Orientations, H=H) #(thetas, phis, chis, H)
                    Bools.append(bools)
                AzboolsTot.append(np.sum(Bools)/len(B))
                #print(ax, ay, az, "             ", np.sum(BF)/len(B), np.sum(BT)/len(B))
            AyboolsTot.append(AzboolsTot)
            print("Ax, Ay, boolean = ", ax, ay, np.sum(AyboolsTot[-1])/len(Az))
        AxboolsTot.append(AyboolsTot)
    return AxCoords, AyCoords, AzCoords, AxboolsTot

def HyperSolChecker1Bool(Orientations, H=HyperH0): #(thetas, phis, chis, H):
    bools = []
    for OrientationVec, Ortho1, Ortho2 in Orientations:
        H.FibbonachiRotateHamil(OrientationVec, Ortho1, Ortho2)
        H.MkEigensystem()
        H.MkTransFrq()
        H.MkTransInt()
        if H.FrqCheckIndependantPulses() and H.TransCheck() == True:
            bools.append(True)
    boolsPercent = np.sum(bools)/len(Orientations)
        #OPTIMISE - check for which statements are least often true. put them earlier to avoid unneccesary checking. 
    #print(np.sum(boolsFrq)/n, "  rp Hits out of ", n)
    #print(np.sum(boolsTrans)/n, "Trans Hits out of ", n)
    return boolsPercent

def HyperPlotter1Bool(Ax, Ay, Az, AxboolsTot, fidelity):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #x Extract the x, y, z coordinates of the solution space,
    success = []
    x_coords = []
    y_coords = []
    z_coords = []
    for i in range(len(Ax)):
        for j in range(len(Ay)):
            for k in range(len(Az)):
                if AxboolsTot[i][j][k]>fidelity:
                    x_coords.append(Ax[i])
                    y_coords.append(Ay[j])
                    z_coords.append(Az[k])
                    success.append(AxboolsTot[i][j][k]*100)
    M = np.argmax(success)
    print("Ax = ", x_coords[M], "Ay = ", y_coords[M], "Az = ", z_coords[M], "Success = ", success[M])
    scatter = ax.scatter(x_coords, y_coords, z_coords, 'D',c=success, cmap=plt.cm.rainbow)
    # Add color bar to show the color map legend
    colorbar = fig.colorbar(scatter, ax=ax, label='% Success Rate for Hyperfine Coupling Paramater')
    colorbar.set_label('% Success Rate for Hyperfine Coupling Paramater', fontsize=24)
    ax.set_xlabel(r'$A_x$ [MHz]', fontsize=24)
    ax.set_ylabel(r'$A_y$ [MHz]', fontsize=24)
    ax.set_zlabel(r'$A_z$ [MHz]', fontsize=24)
    plt.title('Solution Space for a Ytrium Hyperfine System', fontsize=32)
    plt.show()

def BvsOrientationPlotterMultiple(file_paths, labels = None):
    plt.figure()
    for i, file_path in enumerate(file_paths):
        if file_path.endswith('.pkl'):
            with open(file_path, 'rb') as f:
                B, Orientations, Bools = pickle.load(f)
        else:
            data = np.loadtxt(file_path)
            B = data[:, 0]
            Bools = data[:, 1:]
        success_rate = np.sum(Bools, axis=1)
        # Use custom label if provided, otherwise use file path
        print(Orientations[np.argmax(success_rate)], "Max Success Rate = ", np.max(success_rate))
        label = labels[i] if labels and i < len(labels) else file_path
        plt.plot(B, success_rate, label=label)
        print(f"Max Success Rate for {file_path} = {np.max(success_rate)}")
    plt.xlabel('B Field [mT]', fontsize = 22)
    plt.ylabel('Success Rate [%]', fontsize = 22)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.title('Success Rate vs B Field for Multiple Files', fontsize = 26)
    plt.legend(fontsize = 20)
    plt.show()
# BvsOrientationPlotterMultiple(['file1.txt', 'file2.txt', 'file3.txt'])

def HyperBRotChecker1Bool(B, Orientations, H): 
    Bools = []
    OrientationVecs = []
    a = 0
    for b in B:
        H.coil = b
        bools = []
        for OrientationVec, Ortho1Vec, Ortho2Vec in Orientations:
            if a == 0:
                OrientationVecs.append(OrientationVec)
            H.FibbonachiRotateHamil(OrientationVec, Ortho1Vec, Ortho2Vec)
            H.MkEigensystem()
            H.MkTransFrq()
            H.MkTransInt()
            if H.FrqCheckIndependantPulses() and H.TransCheck() == True:
                bools.append(True)
            else:
                bools.append(False)
        a += 1
        Bools.append(bools)
        print("B = ", b, ", Bool = ", np.sum(Bools[-1])/len(Bools[-1]) )
    return Bools, OrientationVecs

def BRotPlotter1Bool(BVals, OrientationVecs, Bools, fidelity):
    Bools = np.array(Bools)
    NewBools = np.zeros((len(OrientationVecs)))
    for i in range(len(Bools[0])):
        NewBools[i] = np.sum(Bools[:, i])*100*0.4 / len(Bools)
    max = np.argmax(NewBools)
    print(OrientationVecs[max], "Max Success Rate = ", np.max(NewBools))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(OrientationVecs)):
        if NewBools[i] > fidelity:
            #print(OrientationVecs[i], NewBools[i])
            #scatter = ax.scatter(OrientationVecs[i][0], OrientationVecs[i][1], OrientationVecs[i][2], marker = 'D', c = NewBools[i], cmap=plt.cm.rainbow)
            ax.scatter(OrientationVecs[i][0], OrientationVecs[i][1], OrientationVecs[i][2], marker = 'D', color = 'black', s=40, zorder =10)
            ax.scatter(OrientationVecs[i][0], OrientationVecs[i][1], -1, marker = 'o', alpha = 0.2, color = 'grey', s = 40, zorder=1)
            ax.scatter(OrientationVecs[i][0], 1, OrientationVecs[i][2], marker = 'o', alpha = 0.2, color = 'grey', s = 40, zorder=1)
            ax.scatter(-1, OrientationVecs[i][1], OrientationVecs[i][2], marker = 'o', alpha = 0.2, color = 'grey', s = 40, zorder=1)
    #colorbar = fig.colorbar(scatter, ax=ax, label='% Success Rate for across Magnetic Field Range')
    #colorbar.set_label('% Success Across Magnetic Field Range', fontsize=18)
    #Max = np.argmax(success)
    ax.set_xlabel(r'$\vec{B_X}/|B|$', fontsize=18)
    ax.set_ylabel(r'$\vec{B_Y}/|B|$', fontsize=18)
    ax.set_zlabel(r'$\vec{B_Z}/|B|$', fontsize=18)
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    #plt.title('Orientations of Magnetic Field with Respect to Hydrogen System', fontsize=20)
    plt.show()

def PullCorrectHyper(Ax, Ay, Az, AxboolsTot):
    for i in range(len(Ax)):
        for j in range(len(Ay)):
            for k in range(len(Az)):
                if i == 0 and j == 0 and k == 0:
                    max = [i,j,k,AxboolsTot[i][j][k][0]]
                else:
                    max = [i,j,k,AxboolsTot[i][j][k][0]] if AxboolsTot[i][j][k][0] > max[3] else max
    print(max)
    print("Ax = ", Ax[max[0]], "Ay = ", Ay[max[1]], "Az = ", Az[max[2]], "Frq = ", max[3], "Trans = ", AxboolsTot[max[0]][max[1]][max[2]][1])

def PullCorrect(BVals, Thetas, Phis, BBolsFrqTot, BboolsTransTot):
    for i in range(len(BVals)):
        for j in range(len(Thetas)):
            if BBolsFrqTot[i][j] and BboolsTransTot[i][j] > 0:
                print("BValue = ", BVals[i], "Theta = ", (Thetas[j]%(2*np.pi))*360/(2*np.pi), "Phi = ", 180*Phis[j]/np.pi)