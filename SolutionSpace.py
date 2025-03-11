from HamiltonianClass import H0
from Consts import Planck, MHzmT, SpinOp4, SpinOp2, SpinOp3, DEcm_ZfsMHZ, rotate_vectors, Orientations, OrthoVecs
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
import sys
import select
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec

# Check if Hamiltonian at each orientation obeys a set of rules, 
# Then append a true false boolean to a list of solution space, 
# That being the Mag Field's orientation,  
# And then plot the solution space.
#Don't go off what wubits exist in lit review. go off what would make sense.
def DESearch1Bool(D, B, gs, Orientations):
    DCoords = []
    ECoords = []
    DboolsTot = []
    for d in D:
        E = np.arange(-abs(d)/3, abs(d)/3, abs(d)/len(D))
        EBoolsTot = []
        eCoords = []
        for e in E:
            eCoords.append(e)
            H = H0(D = DEcm_ZfsMHZ(d, e), gs = np.diag([gs,gs,gs]), coil = 0, SpinOp = SpinOp4, TransIntVal = 0.075, FrqMidVal = 9000, FrqRangeVal = 1000)
            BBoolsTot = []
            for b in B:
                H.coil = b
                Bools = DESolChecker1Bool(Orientations, H)
                BBoolsTot.append(Bools)
            EBoolsTot.append(np.sum(BBoolsTot)/len(BBoolsTot))
            print("D,E,Success = ", round(d,4),round(e,4),EBoolsTot[-1])
        ECoords.append(eCoords)
        DCoords.append(d)
        DboolsTot.append(EBoolsTot)
        print("D, Success = ", d, np.sum(DboolsTot[-1])/len(EBoolsTot))
        # Check for user input without blocking
        print("Press 'stop' to stop the script or wait to continue...")
        if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
            user_input = sys.stdin.readline().strip()
            if user_input.lower() == 'stop':
                print("Stopping the script...")
                break
    return DCoords, ECoords, DboolsTot

def DESolChecker1Bool(Orientations, H=H0):
    bools = []
    for OrientationVec, Ortho1, Ortho2 in Orientations:
        boolt = []
        for angle in np.linspace(0, 2*np.pi, 8):
            Ortho1t, Ortho2t = rotate_vectors(OrientationVec, Ortho1, Ortho2, angle)
            H.FibbonachiRotateHamil(OrientationVec, Ortho1t, Ortho2t)
            H.MkEigensystem()
            H.MkTransFrq()
            H.MkTransInt()
            if H.FreqCheck() and H.TransCheck() == True:
                boolt.append(True)
            else:
                boolt.append(False)
        bools.append(np.sum(boolt)/len(boolt))
    boolsPercent = np.sum(bools)/len(Orientations)
    return boolsPercent

def BOrientSearch1Bool(B = np.arange(0,500,1), NOrient=200, Nangles= 20, D = 0.13, E = 0.0346, gs = 2.0023):
    OrientationVecs = Orientations(NOrient)
    print(D, E, gs)
    H = H0(D = DEcm_ZfsMHZ(D,E), gs = np.diag([gs, gs, gs]), coil = 0, SpinOp = SpinOp4, TransIntVal = 0.27, FrqMidVal = 42000, FrqRangeVal = 8000)
    print("H.D =", H.D)
    print("H.gs =", H.gs)
    BboolsTot = []
    a = 0
    Orientations2 = []
    angles_array = np.linspace(0, 2*np.pi, Nangles)
    for b in B:
        H.coil = b
        bools = []
        for OrientationVec, Ortho1, Ortho2 in OrientationVecs:
            """
            H.FibbonachiRotateHamil(OrientationVec, Ortho1, Ortho2)
            H.MkEigensystem()
            H.MkTransFrq()
            H.MkTransInt()
            if H.FreqCheck() and H.TransCheck() == True:
                bools.append(True)
            else:
                bools.append(False)
            """
            boolt = []
            for angle in angles_array:
                Ortho1t, Ortho2t = rotate_vectors(OrientationVec, Ortho1, Ortho2, angle)
                H.FibbonachiRotateHamil(OrientationVec, Ortho1t, Ortho2t)
                H.MkEigensystem()
                H.MkTransFrq()
                H.MkTransInt()
                if H.FreqCheck() and H.TransCheck() == True:
                    boolt.append(True)
                else:
                    boolt.append(False)
            bools.append(np.sum(boolt)/len(boolt))
            if a == 0:
                Orientations2.append(OrientationVec)
        a +=1
        BboolsTot.append(bools)
        print("B, Success = ", b, np.sum(BboolsTot[-1]))
    return BboolsTot, Orientations2

def BOrientationPlotter1Bool(BVals, BVecs, Bools):
    BVecs = []
    success = []
    for i in range(len(BVals)):
        for j in range(len(BVecs)):
            if Bools[i][j] == True:
                success.append(Bools[i][j])
                BVecs.append(BVals[i]*BVecs[j])
    best = np.argmax(success)
    print("Best Solution: B = ", BVecs[best])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(BVecs, 'o')
    colorbar = plt.colorbar(scatter)
    ax.set_xlabel(r'$B_X')
    ax.set_ylabel(r'B_Y')
    ax.set_zlabel(r'B_Z')
    #plt.title('Solution Space for Orientation and Magnetic  Field')
    plt.show()

def TiFrqvsMagField(H, BVals, OrientationVec):
    Ortho1, Ortho2 = OrthoVecs(OrientationVec)
    Frqs = []
    TIs = []
    BasisMix = []
    for i in BVals:
        H.coil = i
        H.FibbonachiRotateHamil(OrientationVec, Ortho1, Ortho2)
        H.MkEigensystem()
        H.MkTransFrq()
        H.MkTransInt()
        #print(H.TransInt)
        Frqs.append([H.TransFrq[0,1], H.TransFrq[0,2], H.TransFrq[2,3]])
        TIs.append([H.TransInt[0,1], H.TransInt[0,2], H.TransInt[2,3]])
        BasisMix.append([H.Evecs[0,1],H.Evecs[2,1]])
    Frqs = np.array(Frqs)
    TIs = np.array(TIs)
    BasisMix = np.array(BasisMix)
    return Frqs, TIs, BasisMix

def TiFrqvsMagFieldPlot(Ti, Frqs, BasisMix, BVals, Frq_Lower = 8000, Frq_Upper = 10000):
    #Two plots, one of Ti vs MagField, one of Frq vs MagField
    fig, ax = plt.subplots(2,1, sharex=True, figsize=(10, 8))
    #plotting Ti 0,1 0,2 2,3 where Ti is a 4x4 mstrix
    ax[0].plot(BVals, Ti[:,0], label = r'$ | \langle 00 | \hat{S}_x + \hat{S}_y |10 \rangle|^2 $')
    ax[0].plot(BVals, Ti[:,1], label = r'$ | \langle 10 | \hat{S}_x + \hat{S}_y |11 \rangle |^2 $')
    ax[0].plot(BVals, Ti[:,2], label = r'$ | \langle 01 | \hat{S}_x + \hat{S}_y |11 \rangle |^2 $')
    Threshold = 0.15*np.ones(len(BVals))
    ax[0].fill_between(BVals, Threshold, color = 'red', alpha = 0.3)
    ax[0].plot(BVals, Threshold, color = 'black')
    B_above_Lower = BVals[Frqs[:,0] > Frq_Lower]
    B_below_Upper = BVals[Frqs[:,1] < Frq_Upper]
    #print(B_above_Lower)
    #print(B_below_Upper)

    # Draw green vertical shaded region between B_above_8k and B_below_10k
    """
    if B_above_Lower.size > 0 and B_below_Upper.size > 0:
        ax[0].axvspan(B_above_Lower[0], B_below_Upper[-1], color='green', alpha=0.3)
        ax[1].axvspan(B_above_Lower[0], B_below_Upper[-1], color='green', alpha=0.3)
        
    
    #Vert line that boarders the green region. 
    ax[0].axvline(B_above_Lower[0], color='black', linestyle='--')
    ax[0].axvline(B_below_Upper[-1], color='black', linestyle='--')
    ax[1].axvline(B_above_Lower[0], color='black', linestyle='--')
    ax[1].axvline(B_below_Upper[-1], color='black', linestyle='--')
    """
    TenkHz = Frq_Upper*np.ones(len(BVals))
    EightkHz = Frq_Lower*np.ones(len(BVals))
    HzMax = max(Frqs[:,1])*np.ones(len(BVals))
    ax[1].tick_params(axis='x', which='major', labelsize=16)
    ax[1].plot(BVals, Frqs[:,0], label = r'$\omega_{10-00}$')
    ax[1].plot(BVals, Frqs[:,1], label = r'$\omega_{11-10}$')
    ax[1].plot(BVals, Frqs[:,2], label = r'$\omega_{11-01}$')
    ax[1].fill_between(BVals, EightkHz, color = 'red', alpha = 0.3)
    ax[1].fill_between(BVals, HzMax, TenkHz, color = 'red', alpha = 0.3)
    ax[1].plot(BVals, TenkHz, color = 'black')
    ax[1].plot(BVals, EightkHz, color = 'black')
    # Remove minor ticks
    ax[0].minorticks_off()
    
    # Set x-axis limits to remove extra space
    ax[0].set_xlim(0, 300)
    ax[1].set_xlim(0, 300)
    
    ax[0].set_ylabel('Transition Intensity [a.u.]', size = 16)
    ax[1].set_xlabel('Magnetic Field [mT]', size = 18)
    ax[1].set_ylabel('Transition Frequency [MHz]', size = 16)
    # Add legends to each subplot
    ax[0].legend(fontsize = 14)
    ax[1].legend(fontsize = 14)
    
    plt.tight_layout()
    plt.show()

def TiFrqvsOrientation(H, Orientations):
    Frqs = []
    TIs = []
    vecs = []
    index = []
    for Orientation, Ortho1, Ortho2 in Orientations:
        H.FibbonachiRotateHamil(Orientation, Ortho1, Ortho2)
        H.MkEigensystem()
        H.MkTransFrq()
        H.MkTransInt()
        if H.FreqCheck() and H.TransCheck() == True:
            index.append(True)
        else:
            index.append(False)
        Frqs.append([H.TransFrq[0,1], H.TransFrq[0,2], H.TransFrq[2,3]])
        TIs.append([H.TransInt[0,1],  H.TransInt[0,2], H.TransInt[2,3]])
        vecs.append(Orientation)
    Frqs = np.array(Frqs)
    TIs = np.array(TIs)
    vecs = np.array(vecs)
    index = np.array(index)
    print(np.sum(index))
    return Frqs, TIs, vecs, index

def TiFrqvsOrientationPlot(Ti, Frqs, vecs, index):
    # Create a figure with a GridSpec layout
    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    # Create the 3D subplot for the Fibonacci sphere on the LHS
    ax_sphere = fig.add_subplot(gs[:, 0], projection='3d')
    cmap = plt.get_cmap('viridis')
    norm = Normalize(vmin=0, vmax=len(vecs))
    for i in range(len(vecs)):
        ax_sphere.scatter(vecs[i, 0], vecs[i, 1], vecs[i, 2], color=cmap(norm(i)), alpha=0.5)
        if index[i] == True:
            ax_sphere.scatter(vecs[i, 0], vecs[i, 1], vecs[i, 2], color='green', s=100)
            ax_sphere.scatter(-1, vecs[i, 1], vecs[i, 2], color='grey', alpha=0.4, s=80)
            ax_sphere.scatter(vecs[i, 0], 1, vecs[i, 2], color='grey', alpha=0.4, s=80)
            ax_sphere.scatter(vecs[i, 0], vecs[i, 1], -1, color='grey', alpha=0.4, s=80)
    ax_sphere.set_title(r'Mag Vector Relative to Molecular Z-axis $\vec{B}/|\vec{B}|$', size=20)
    ax_sphere.set_xlabel(r'$B_{X}$')
    ax_sphere.set_ylabel(r'$B_{Y}$')
    ax_sphere.set_zlabel(r'$B_{Z}$')
    # Create the first subplot for Ti on the top RHS
    TiMax = np.max(Ti)
    ax0 = fig.add_subplot(gs[0, 1])
    ax0.plot(range(len(vecs)), Ti[:,0], label=r'$| \langle 00 | S_x + S_y |01 \rangle |^{2} $', alpha=0.5, color='blue')
    ax0.plot(range(len(vecs)), Ti[:,1], label=r'$| \langle 00 | S_x + S_y |10 \rangle |^{2} $', alpha=0.5, color='yellow')
    ax0.plot(range(len(vecs)), Ti[:,2], label=r'$| \langle 10 | S_x + S_y |11 \rangle |^{2} $', alpha=0.5, color = 'purple')
    Threshold = 0.15 * np.ones(len(vecs))
    ax0.fill_between(range(len(vecs)), 0, Threshold, color='red', alpha=0.3)
    ax0.set_ylabel('Transition Intensity')
    ax0.set_xticks([])
    ax0.legend()
    # Create the second subplot for Frqs on the bottom RHS
    ax1 = fig.add_subplot(gs[1, 1])
    Frq_Higher = 10000 * np.ones(len(vecs))
    Frq_Lower = 8000 * np.ones(len(vecs))
    HzMax = max(Frqs[:,1]) * np.ones(len(vecs))
    ax1.plot(range(len(vecs)), Frqs[:,0], label=r'$\omega_{01-00}$', alpha=0.5, color='blue')
    ax1.plot(range(len(vecs)), Frqs[:,1], label=r'$\omega_{10-00}$', alpha=0.5, color='yellow')
    ax1.plot(range(len(vecs)), Frqs[:,2], label=r'$\omega_{11-10}$', alpha=0.5, color='purple')
    ax1.fill_between(range(len(vecs)), 0, Frq_Lower, color='red', alpha=0.3)
    ax1.fill_between(range(len(vecs)), Frq_Higher, HzMax, color='red', alpha=0.3)
    # Add vertical bands where index is True
    for i in range(len(index)):
        if index[i] == True:
            ax1.fill_betweenx([0,np.max(HzMax)], i-0.5, i+0.5, color='green', alpha=0.3)
            ax0.fill_betweenx([0, TiMax], i-0.5, i+0.5, color='green', alpha=0.3)
    # Finals for plots.
    ax1.set_xlabel('Orientation Index')
    ax1.set_ylabel('Transition Frequency')
    ax1.legend()
    
    # Add a color bar below the plots
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax_sphere, orientation='horizontal', fraction=0.02, pad=0.1)
    cbar.set_label('Orientation Index')
    
    plt.tight_layout()
    plt.show()

def DEPlotter1Bool(DCoords,ECoords,DboolsTot, Verbose = 0):
    List = []
    D = []
    E = []
    bools = []
    for i in range(len(DCoords)):
        for j in range(len(ECoords[i])):
            if DboolsTot[i][j] > 0:
                D.append(DCoords[i])
                E.append(ECoords[i][j])
                bools.append(DboolsTot[i][j]*100)
                if Verbose == 1:
                    List.append([DCoords[i], ECoords[i][j], DboolsTot[i][j]])
    scatter = plt.scatter(D, E, c=bools, cmap='rainbow', linewidths=1)
    D_line = np.linspace(min(DCoords), max(DCoords), 10)
    E_line = D_line / 3
    plt.plot(D_line, E_line, color='black', label = '|E| = |D/3|', linewidth = 3)
    plt.plot(D_line, -E_line, color='black', linewidth = 3)
    plt.scatter(0.0977, 0.025, color = 'black', marker = 'x', linewidths=3, s=100)
    if Verbose == 1:
        if List == []:
            print("No successful solutions")
        else:
            Success = np.array(List)[:,2]
            MaxIndex = np.argmax(Success)
            print(List[MaxIndex])
            print("Total Possibilities = ",len(List))
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    colorbar = plt.colorbar(scatter, label = '% Success Metric')
    colorbar.set_label('% Success Metric', fontsize = 20)
    plt.legend(by_label.values(),by_label.keys(), fontsize = 18, loc = 'upper center')
    plt.xlabel(r'D [$cm^{-1}$]', fontsize = 20)
    plt.ylabel(r'E [$cm^{-1}$]', fontsize = 20)
    #plt.title('Solution Space for D and E', fontsize = 32)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.show()

#H = H0(D = DEcm_ZfsMHZ(0.13, -0.0346), gs = np.diag([2.0023, 2.0023, 2.0023]), coil = 162, SpinOp = SpinOp4, TransIntVal = 0.075, FrqMidVal = 9000, FrqRangeVal = 1000)
#Orientations = Orientations(1600)
#Frqs, TIs, vecs, index = TiFrqvsOrientation(H, Orientations)
#TiFrqvsOrientationPlot(TIs, Frqs, vecs, index)