import numpy as np
from SolutionSpace import DESearch1Bool, BOrientSearch1Bool, DEPlotter1Bool, BOrientSearch1Bool, TiFrqvsMagField, TiFrqvsMagFieldPlot, TiFrqvsOrientation, TiFrqvsOrientationPlot, BOrientationPlotter1Bool
from HyperfineHamilClass import HyperH0
from HamiltonianClass import H0
from SolutionSpaceHyperfine import HyperSearch, PullCorrectHyper, HyperPlotter1Bool, HyperBRotChecker1Bool, BRotPlotter1Bool, RotatingB1Data, BvsOrientationPlotterMultiple
from Consts import SpinOp2, Orientations, FibbonaciSphere, DEcm_ZfsMHZ, SpinOp4, SpinOp3
import pickle
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

H = HyperH0(A = np.diag([-5598, -6008, -3180]), coil = 340, SpinOp = SpinOp2, NucSpinOp=SpinOp2, gi=np.diag([1.335, 1.410, 1.693]), gs=np.diag([1.335, 1.410, 1.693]))
#H = HyperH0(A = np.diag([-4929, -4929, -4929]), coil = 340, SpinOp = SpinOp2, NucSpinOp=SpinOp2, gi=np.diag([0,0,0]), gs=np.diag([1.479, 1.479, 1.479]))
#H = HyperH0(A=np.diag([-10,10,-10]), coil = 340, SpinOp = SpinOp2, NucSpinOp=SpinOp2, gi=np.diag([2,2,2]), gs=np.diag([2,2,2]))

#H.FreePlotNegative(0, 10, 1)
H.FreePlot(0, 3000, 1)
exit()
 
 # hmmmm

Hydrogen = 5.58569468
Ytrium = -0.2748308
gs = 1.99


D = 0.13
E = -0.0346
#D=0.0977
#E=0.025

B = np.arange(0, 300, 1)
H = H0(D = DEcm_ZfsMHZ(D,E), gs = np.diag([gs, gs, gs]), coil = 0, SpinOp = SpinOp4)
Frqs, TIs, BasisMix = TiFrqvsMagField(H, B, OrientationVec = np.array([0.74888958, 0.05988857, 0.65998315]))
TiFrqvsMagFieldPlot(TIs, Frqs, BasisMix, B)
exit()


with open('D-0.150.1450.005_B1002455_N200QuditFull.pkl', 'rb') as f:
    DCoords, ECoords, DboolsFrqTot = pickle.load(f)
DEPlotter1Bool(DCoords, ECoords, DboolsFrqTot, Verbose = 1)
exit("Plotted?")



with open('IdealDE0.13_-0.0346_B1002991_N1600.pkl', 'rb') as f:
    B, OrientationVecs, Bools = pickle.load(f)
BRotPlotter1Bool(B, OrientationVecs, Bools, 0.0)
exit("Plotted?")


D=0.13
E=-0.0346
B = np.arange(50, 400, 1)
NOrient = 3200
Bools, OrientationVecs = BOrientSearch1Bool(B, NOrient, Nangles = 8, D=D, E=E, gs=gs)
try:
    with open('RealDE'+str(D)+'_'+str(E)+'_B'+str(B[0])+str(B[-1])+str(B[2]-B[1])+'_N'+str(NOrient)+'.pkl', 'wb') as f:
        pickle.dump((B, OrientationVecs, Bools), f)
    print("File saved successfully")
except Exception as e:
    print("Error saving file:", e)
BRotPlotter1Bool(B, OrientationVecs, Bools, 0.0)
exit("Well Done!")

with open('IdealDE0.611_0.1572_B1009991_N1600.pkl', 'rb') as f:
    B, OrientationVecs, Bools = pickle.load(f)
#BRotPlotter1Bool(B, OrientationVecs, Bools, 0)
BOrientationPlotter1Bool(B, OrientationVecs, Bools)
exit("Plotted?")

NOrient = 200 #Now with beetter distribution of Points that dont favour the poles.
D = np.arange(-.15, -.05, .005)
D = np.append(D, np.arange(.05, .15, .005))
#print(len(D))
#exit()
B = np.arange(100, 250, 5)
OrientationsVecs = Orientations(NOrient)
DCoords, ECoords, Bools = DESearch1Bool(D, B, gs, OrientationsVecs)
try:
    #with open('DEWithRolland1Bool2.pkl', 'wb') as f:
    with open('D'+str(round(D[0],3))+str(round(D[-1],3))+str(round(D[2]-D[1],3))+'_B'+str(B[0])+str(round(B[-1],3))+str(round(B[2]-B[1],3))+'_N'+str(NOrient)+'QuditFull.pkl', 'wb') as f:
        pickle.dump((DCoords, ECoords, Bools), f)
    print("File saved successfully")
except Exception as e:
    print("Error saving file:", e)
DEPlotter1Bool(DCoords, ECoords, Bools)
exit("OMG thank god, wtf is wrong with you.")


B = np.arange(300, 1500, 2)
NOrient = 800
Bools, OrientationVecs = BOrientSearch1Bool(B, NOrient, Nangles = 8, D=D, E=E, gs=gs)
try:
    with open('RealDE'+str(D)+'_'+str(E)+'_B'+str(B[0])+str(B[-1])+str(B[2]-B[1])+'_N'+str(NOrient)+'.pkl', 'wb') as f:
        pickle.dump((B, OrientationVecs, Bools), f)
    print("File saved successfully")
except Exception as e:
    print("Error saving file:", e)
BRotPlotter1Bool(B, OrientationVecs, Bools, 0.0)
exit("Is it over? cann i go now?")

#H = H0(D = DEcm_ZfsMHZ(D,E), gs = np.diag([1.99,1.99,1.99]), coil = 0, SpinOp = SpinOp4)

Nucleus = Ytrium
B = np.arange(100, 500, 1)
NOrient = 800
AxAyAz = np.diag([300, 180, 30])
OrientationsZipList = Orientations(NOrient)
H = HyperH0(A = AxAyAz, coil = 0, SpinOp = SpinOp2, NucSpinOp=SpinOp2, gi=np.diag([Nucleus, Nucleus, Nucleus]), gs=np.diag([gs,gs,gs]))
Bools, OrientationsVec = HyperBRotChecker1Bool(B, OrientationsZipList, H)
try:
    with open(str(Nucleus)+'_AxAyAz'+str(AxAyAz[0,0])+str(AxAyAz[1,1])+str(AxAyAz[2,2])+'_B'+str(B[0])+str(B[-1])+str(B[2]-B[1])+'_N'+str(NOrient)+'_2Qubit.pkl', 'wb') as f:
        pickle.dump((B, OrientationsVec, Bools), f)
    print("File saved successfully")
except Exception as e:
    print("Error saving file:", e)
BRotPlotter1Bool(B, OrientationsVec, Bools, 0)
exit()

#-0.2748308_AxAy039010_Az039010_B20039010_N100_2Qubit.pkl - 5.58569468_AxAy039010_Az039010_B20039010_N100_2Qubit.pkl
with open('5.58569468_AxAy039010_Az039010_B20039010_N100_2Qubit.pkl', 'rb') as f:
    Ax, Ay, Az, AxboolsTot = pickle.load(f)
print("File loaded successfully")
HyperPlotter1Bool(Ax, Ay, Az, AxboolsTot, fidelity=0.0)
#print(AxboolsTot)
exit()

Ax = np.arange(0, 400, 10)
Ay = np.arange(0, 400, 10)
Az = np.arange(0, 400, 10)
B = np.arange(200, 400, 10)
NOrient = 100
OrientationVecs = Orientations(NOrient) #This is a list of tuples of the form (OrientationVector, Ortho1Vector, Ortho2Vector)
Nucleus = Ytrium
Ax, Ay, Az, AxboolsTot = HyperSearch(Ax, Ay, Az, B, Orientations=OrientationVecs, gi = Nucleus)
try:
    with open(str(Nucleus)+'_AxAy'+str(Ax[0])+str(Ax[-1])+str(Ax[2]-Ax[1])+'_Az'+str(Az[0])+str(Az[-1])+str(Az[2]-Az[1])+'_B'+str(B[0])+str(B[-1])+str(B[2]-B[1])+'_N'+str(NOrient)+'_2Qubit.pkl', 'wb') as f:
        pickle.dump((Ax, Ay, Az, AxboolsTot), f)
    print("File saved successfully")
except Exception as e:
    print("Error saving file:", e)
HyperPlotter1Bool(Ax, Ay, Az, AxboolsTot, fidelity=0)
exit("Search Complete")

# BvsOrientationPlotterMultiple(['file1.txt', 'file2.txt', 'file3.txt'])
filepaths = ['-0.2748308_AxAyAz3107010_B2004982_N800_NewRotation.pkl',
            '5.58569468_AxAyAz2706010_B2003982_N800_NewRotation.pkl',
            'D013_E00346_B1005002_N256_Roll_1Bool.pkl']
labels = ['Ytrium', 'Hydrogen', 'Chromium']
BvsOrientationPlotterMultiple(filepaths, labels)
exit()


H = HyperH0(A = np.diag([-10, 150, -160]), coil = 320, SpinOp = SpinOp2, NucSpinOp=SpinOp2, gi=np.diag([Hydrogen, Hydrogen, Hydrogen]), gs=np.diag([gs,gs,gs]))
RotatingB1Data(np.array([ 0.95228413, -0.02568605, -0.30413016]), H, NOrient = 100)
exit()


with open('D-0.20.1950.005_B107955_N200NewRotation.pkl','rb') as f:
    DCoords, ECoords, Bools = pickle.load(f)
DEPlotter1Bool(DCoords, ECoords, Bools)
exit()


with open('DESearchWithRoll.pkl', 'rb') as f:
    DCoords, ECoords, DboolsFrqTot, DboolsTransTot = pickle.load(f)
DEPlotter(DCoords, ECoords, DboolsFrqTot, DboolsTransTot, Verbose = 1)
exit()


OrientationVectors = Orientations(10)
H1 = HyperH0(A = np.diag([12, 16, 27]), coil = 340, SpinOp = SpinOp2, NucSpinOp=SpinOp2, gi=np.diag([Hydrogen, Hydrogen, Hydrogen]), gs=np.diag([gs,gs,gs]))
H2 = HyperH0(A = np.diag([16, 12, 27]), coil = 340, SpinOp = SpinOp2, NucSpinOp=SpinOp2, gi=np.diag([Hydrogen, Hydrogen, Hydrogen]), gs=np.diag([gs,gs,gs]))
for Orientation, Ortho1, Ortho2 in OrientationVectors:
    print("Orientation = ", Orientation)
    H1.FibbonachiRotateHamil(Orientation, Ortho1, Ortho2)
    H1.MkEigensystem()
    H1.MkTransFrq()
    H1.MkTransInt()
    H2.FibbonachiRotateHamil(Orientation, Ortho1, Ortho2)
    H2.MkEigensystem()
    H2.MkTransFrq()
    H2.MkTransInt()
    print("Evals 1 = \n",np.round(H1.Evals,3))
    print("Evals 2 = \n", np.round(H2.Evals,3))
    input()
    print("Evecs 1 = \n", np.round(H1.Evecs,3))
    print("Evecs 2 = \n", np.round(H2.Evecs,3))
    input()
    print("Trans 1 = \n", np.round(H1.TransInt,3))
    print("Trans 2 = \n, ", np.round(H2.TransInt,3))
    input()
exit()