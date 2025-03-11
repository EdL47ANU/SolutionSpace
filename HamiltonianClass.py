import numpy as np
import matplotlib.pyplot as plt
from Consts import MHzmT, SpinOp4, SpinOp2, SpinOp3, DEcm_ZfsMHZ, PolarToCartesian, unit_vector_to_spherical_coords
print("Starting ZFS....     ----------------------------------------------------------")

# B field with sinesoidal patterns for different oriantations, return the collection of diagonalised matrices' energy differences, 
# and the probability of the transition.
class H0:
    def __init__(self, 
                 D = np.diag([600, 800, 1200]), #in MHz, use DEcm_ZfsMHZ to convert from D,E in cm^-1 to a matrix in MHz,
                 gs = np.diag([2.0, 2.0, 2.0]),
                 coil = 160,        # In mT
                 SpinOp = SpinOp4, 
                 TransIntVal = 0.075,
                 FrqMidVal = 9000,
                 FrqRangeVal = 1000,
                 OrientationVec = np.array([0,0,1])):
        self.D = D
        SpinOpLen = len(SpinOp[0])
        self.gs = gs
        self.SpinOp = 0.5*SpinOp
        self.ZFS = np.einsum('iab,ij,jbc -> ac', self.SpinOp, self.D, self.SpinOp) # D matrix is expressed in MHz already.
        #self.ZFS = np.einsum('iab,ij,jbc -> ac', self.SpinOp, D, self.SpinOp) # D matrix is expressed in MHz already.
        # print("ZFS = \n", self.ZFS)
        self.coil = coil
        self.Rotation = np.zeros((3,3))
        self.Evals = np.zeros(SpinOpLen)
        self.Evecs = np.zeros((SpinOpLen, SpinOpLen), dtype=complex)
        self.TransFrq = np.zeros((SpinOpLen, SpinOpLen))
        self.Zeeman=np.zeros((SpinOpLen,SpinOpLen))
        self.TransInt = np.zeros((SpinOpLen, SpinOpLen))
        self.TransFrq = np.zeros((SpinOpLen, SpinOpLen))
        self.Pulse = 0.5*SpinOp
        self.TransIntVal = TransIntVal
        self.FrqMidVal = FrqMidVal
        self.FrqRangeVal = FrqRangeVal
        self.OrientationVec = OrientationVec

    def FibbonachiRotateHamil(self, OrientationVec, Ortho1, Ortho2):
        B = OrientationVec*self.coil
        # Find a vector not parallel to OrientationVec
        self.Pulse_Check = np.einsum('iab,i -> ab', self.Pulse, Ortho1) + np.einsum('iab,i -> ab', self.Pulse, Ortho2)
        self.Zeeman = np.einsum('iab,ji,j -> ab', self.SpinOp, self.gs, B) * MHzmT
        return self.Zeeman, self.ZFS
    
    # Diagonalises the Hamiltonain and rearranges EVals & Evecs in DESCENDING order of energy.
    def MkEigensystem(self):
        self.Evals, self.Evecs = np.linalg.eig(self.ZFS + self.Zeeman)
        #self.Evals = np.real(self.Evals)
        indices = np.argsort(self.Evals)[::-1] #sorted biggest to smalles to immitate a Z axis Zeeman from +3/2 first (top left element) to -3/2 last (Bottom right element)
        self.Evals = self.Evals[indices]
        self.Evecs = self.Evecs[indices,:]
        return self.Evals, self.Evecs
    
    # Calculates the transition frequencies between the energy levels.
    def MkTransFrq(self):
        for i in range(len(self.Evals)):
            for j in range(len(self.Evals)):
                if i != j:
                    self.TransFrq[i,j] = abs(self.Evals[i] - self.Evals[j])
        return self.TransFrq 
    
    # Calculates the transition intensities between the energy levels through the Eigenvalues over a Sx + Sy MW Field.
    def MkTransInt(self):
        for i in range(len(self.Evecs)):
            for j in range(len(self.Evecs)):
                if i != j:
                    #a = abs(np.einsum('i,ij,j ->',np.array([1,0,0,0]),self.Pulse[0]+self.Pulse[1],np.array([0,1,0,0])))**2
                    #print(a)
                    self.TransInt[i,j] = abs(np.einsum('i,ij,j ->',np.conj(self.Evecs[i]),self.Pulse_Check,self.Evecs[j]))**2
        return self.TransInt
    
    def FreqCheck(self):
        if abs(self.TransFrq[2,3]-self.TransFrq[0,1]) > 100:
            if np.isclose(self.TransFrq[0,2], self.FrqMidVal, atol = self.FrqRangeVal):
                if np.isclose(self.TransFrq[0,1], self.FrqMidVal, atol = self.FrqRangeVal):
                    if np.isclose(self.TransFrq[2,3], self.FrqMidVal, atol = self.FrqRangeVal):
                        return True
        else:
            return False

    def TransCheck(self):
        if self.TransInt[0,1]> self.TransIntVal:
            if self.TransInt[0,2] > self.TransIntVal:
                if self.TransInt[2,3] > self.TransIntVal:
                    return True
        else:
            return False

    def PlotEvecs(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        colourmap = {0:'r', 1:'g', 2:'b', 3:'y'}
        for i in range(len(self.Evecs)):
            for j in range(len(self.Evecs[i])):
                colour = colourmap[i]
                ax.bar3d(i,j,0,1,1,self.Evecs[i,j], color=colour)
        plt.show()
        
    def SplittingPlot(self, Bmin,Bmax,Bstep,theta=0, phi=0):
        Evals = []
        zerofieldBVals = range(-int(round(Bmax/10,0)),Bmin,Bstep)
        for i in zerofieldBVals:
            self.coil = np.array([0,0,0])
            self.RotHamil(theta,phi)
            self.MkEigensystem()
            Evals.append(self.Evals)
        BSample = range(Bmin,Bmax,Bstep)
        PlotB=range(min(zerofieldBVals), Bmax, Bstep)
        for i in BSample:
            self.coil = i
            self.MkEigensystem()
            Evals.append(self.Evals)
        plt.ylabel("Energy (MHz)")
        plt.xlabel("B Field (mT)")
        plt.title("Energy Splitting of Ideal Molecule for Super Dense Coding")
        plt.plot(PlotB, Evals)
        plt.show()
        
    def TransmissionFreqPlot(self, Bmin, BMax, BStep):
        Frqs = []
        BSample = range(Bmin,BMax,BStep)
        for i in BSample:
            self.coil = np.array([0,0,i])
            self.RotHamil(np.pi/4,np.pi/4)
            self.MkEigensystem()
            self.MkTransFrq()
            Temp = []
            for k in range(len(self.Evals)):
                for j in range(len(self.Evals)):
                    if k > j:
                        Temp.append(self.TransFrq[k][j])
            Frqs.append(Temp)
        plt.plot(BSample, Frqs)
        plt.show()
        
    def FreePlot(self, Bmin,Bmax,Bstep, MagVec = np.array([0,0,1])):
        Evals = []
        zerofieldBVals = range(-int(round(Bmax/4,0)),Bmin,Bstep)
        HighFieldBVals = np.arange(Bmax, Bmax*3, 1)
        for i in zerofieldBVals:
            self.coil = 0
            self.FibbonachiRotateHamil(MagVec, np.array([1,0,0]), np.array([0,1,0]))
            self.MkEigensystem()
            Evals.append(self.Evals)
        BSample = range(Bmin,Bmax,Bstep)
        PlotB=range(min(zerofieldBVals), max(HighFieldBVals)+1, Bstep)
        for i in BSample:
            self.coil = i
            self.FibbonachiRotateHamil(MagVec, np.array([1,0,0]), np.array([0,1,0]))
            self.MkEigensystem()
            Evals.append(self.Evals)
        for i in HighFieldBVals:
            Evals.append(self.Evals)
        plt.ylabel("Energy [MHz]", size = 26)
        plt.yticks([-10000, -8000, -6000,-4000,-2000, 0, 2000, 4000, 6000, 8000, 10000], size = 18)
        plt.xlabel("Magnetic Field [mT], Time", size = 26)
        
        # Define custom tick positions and labels
        tick_positions = [0,50,100,150,200, 250, 300, 350, 400]
        tick_labels = [str(pos) for pos in tick_positions[:4]] + ['Time', 'evolution', 'at', 'static', 'field.']
        plt.xticks(tick_positions, tick_labels, size=18)
        
        #plt.xticks([0,50,100,162], size = 18) 
        #plt.title("Energy Splitting of Ideal Molecule for Super Dense Coding", size = 32)
        plt.plot(PlotB, Evals, color = 'black')
        plt.show()


#H = H0(D = DEcm_ZfsMHZ(0.13,-0.0326), gs = np.diag([2,2,2]), coil = 162, SpinOp = SpinOp4, )

#H.FreePlot(0, 162, 1, MagVec = np.array([0.75, 0.0622, 0.658]))


# [0.629, 0.0199, 0.777] 0.13 0.0326 768