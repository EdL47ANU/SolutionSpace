import numpy as np
import matplotlib.pyplot as plt
from Consts import MHzmT, SpinOp4, SpinOp2, SpinOp3, DEcm_ZfsMHZ, Orientations
print("Starting Hyperfine....     ----------------------------------------------------------")

# B field with sinesoidal patterns for different oriantations, return the collection of diagonalised matrices' energy differences, 
# and the probability of the transition.
class HyperH0:
    def __init__(self, #in MHz, use DEcm_ZfsMHZ to convert from D,E in cm^-1 to a matrix in MHz,
                 gs = np.diag([2.0, 2.0, 2.0]),
                 A = np.diag([0, 0, 0]), 
                 gi = np.diag([2.0, 2.0, 2.0])/1836,
                 coil = 100,        # In mT
                 SpinOp = SpinOp2,
                 NucSpinOp = SpinOp2):
        self.A = A
        self.gs = gs
        self.gi = gi
        self.coil = coil
        self.SpinOp = 0.5*SpinOp
        self.NucSpinOp = 0.5*NucSpinOp
        MatSpinLen = len(SpinOp[0])*len(NucSpinOp[0])
        self.TransInt = np.zeros((MatSpinLen, MatSpinLen))
        self.TransFrq = np.zeros((MatSpinLen, MatSpinLen))
        self.Pulse = np.array([np.kron(np.eye(len(self.SpinOp[0])),self.NucSpinOp[0]) + np.kron(self.NucSpinOp[0], np.eye(len(self.SpinOp[0]))),
                               np.kron(np.eye(len(self.SpinOp[1])), self.NucSpinOp[1]) + np.kron(self.NucSpinOp[1], np.eye(len(self.SpinOp[1]))),
                               np.kron(np.eye(len(self.SpinOp[2])),self.NucSpinOp[2]) + np.kron(self.NucSpinOp[2], np.eye(len(self.SpinOp[2])))])
        def SpinKron(self):
            SpinKron = np.zeros((3,3, len(self.SpinOp[0])*len(self.NucSpinOp[0]), len(self.SpinOp[0])*len(self.NucSpinOp[0]) ), dtype=complex)
            for i in range(len(self.SpinOp)):
                for j in range(len(self.NucSpinOp)):
                    SpinKron[i,j] = np.kron(np.conj(self.SpinOp[i]), self.NucSpinOp[j])
            return SpinKron
        self.SpinKron = SpinKron(self)
        self.Hyperfine = np.einsum('ij,jiab->ab', self.A, self.SpinKron) *MHzmT # Convert from mT to MHz
        #print(self.SpinKron)

    def ZaxisHamil(self):
        self.Zeeman = np.kron(np.einsum('iab,ji,j -> ab', np.conj(self.SpinOp), self.gs, self.coil), np.eye(len(self.NucSpinOp[0]))) * MHzmT + np.kron(np.eye(len(self.SpinOp[0])), np.einsum('iab,ji,j -> ab', self.NucSpinOp, self.gi, self.coil))*MHzmT/1832 # Convert from mT to MHz
    
    def FibbonachiRotateHamil(self, OrientationVec, Ortho1, Ortho2):
        B = OrientationVec * self.coil
        # Find a vector not parallel to OrientationVec
        self.Pulse_Check = np.einsum('iab,i -> ab', self.Pulse, Ortho1) + np.einsum('iab,i -> ab', self.Pulse, -Ortho2)
        self.Zeeman = np.kron(np.einsum('iab,ji,j -> ab', self.SpinOp, self.gs, B),np.eye(2))*MHzmT + np.kron(np.eye(2),np.einsum('iab,ji,j -> ab', self.NucSpinOp, self.gi, B))*MHzmT/1836 # Convert from mT to MHz
  
    # Diagonalises the Hamiltonain and rearranges EVals & Evecs in DESCENDING order of energy.
    def MkEigensystem(self):
        self.Evals, self.Evecs = np.linalg.eig(self.Hyperfine + self.Zeeman)
        indices = np.argsort(self.Evals) #sorted biggest to smalles to immitate a Z axis Zeeman from +3/2 first (top left element) to -3/2 last (Bottom right element)
        self.Evals = self.Evals[indices]
        self.Evecs = self.Evecs.T[indices,:] # Adding this transpose was the FIX?!

    # Calculates the transition frequencies between the energy levels.
    def MkTransFrq(self):
        for i in range(len(self.Evals)):
            for j in range(len(self.Evals)):
                if i != j:
                    self.TransFrq[i,j] = abs(self.Evals[i] - self.Evals[j])
        #print(self.TransFrq)

    def MkTransFrqCheat(self):
        self.TransFrq[0,1] = abs(self.Evals[1] - self.Evals[0])
        self.TransFrq[0,2] = abs(self.Evals[2] - self.Evals[0])
        self.TransFrq[2,3] = abs(self.Evals[3] - self.Evals[2])

    # Calculates the transition intensities between the energy levels through the Eigenvalues over a Sx + Sy MW Field.
    def MkTransInt(self):
        for i in range(len(self.Evals)):
            for j in range(len(self.Evals)):
                if i != j:
                    self.TransInt[i,j] = abs(np.einsum('i,ji,j ->',np.conj(self.Evecs[i]),self.Pulse_Check, self.Evecs[j]))**2
        #print(self.TransInt)
    
    def MkTransIntCheat(self):
        self.TransInt[0,1]  = abs(np.einsum('i,ij,j ->',np.conj(self.Evecs[0]),self.Pulse_Check,self.Evecs[1])**2)
        self.TransInt[0,2] = abs(np.einsum('i,ij,j ->',np.conj(self.Evecs[0]),self.Pulse_Check,self.Evecs[2])**2)
        self.TransInt[2,3] = abs(np.einsum('i,ij,j ->',np.conj(self.Evecs[2]),self.Pulse_Check,self.Evecs[3])**2)

    def PlotEvecs(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        colourmap = {0:'r', 1:'g', 2:'b', 3:'y'}
        for i in range(len(self.Evecs)):
            for j in range(len(self.Evecs[i])):
                colour = colourmap[i]
                ax.bar3d(i,j,0,1,1,self.Evecs[i,j], color=colour)
        plt.show()
        
    def FrqCheckIndependantPulses(self):
        if np.isclose(self.TransFrq[1,3], 9000, atol = 1000):
            if np.isclose(self.TransFrq[0,2], 9000, atol = 1000):
                if np.isclose(self.TransFrq[0,1], 210, atol = 190):
                        if ~np.isclose(self.TransFrq[0,2], self.TransFrq[1,3], atol = 200):
                            return True
        else:
           return False

    def TransCheck(self):
        if self.TransInt[0,1]> 0.075:
            if self.TransInt[1,3] > 0.075:
                if self.TransInt[0,2] > 0.075:
                    #print("Pass")
                    return True
        else:
            return False

    def SplittingPlot(self, Bmin,Bmax,Bstep,theta=0, phi=0, chi=0):
        Evals = []
        zerofieldBVals = range(-int(round(Bmax/10,0)),Bmin,Bstep)
        for i in zerofieldBVals:
            self.coil = np.array([0,0,0])
            self.RotHamil(theta,phi,chi)
            self.MkEigensystem()
            Evals.append(self.Evals)
        BSample = range(Bmin,Bmax,Bstep)
        PlotB=range(min(zerofieldBVals), Bmax, Bstep)
        for i in BSample:
            self.coil = np.array([0,0,i])
            self.RotHamil(theta=theta,phi=phi, chi=chi)
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
            self.RotHamil(np.pi/4,np.pi/4,np.pi/4)
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

    def FreePlot(self, Bmin,Bmax,Bstep,theta=0, phi=0, chi=0):
        Evals = []
        Evecs = []
        zerofieldBVals = range(-int(round(Bmax/10,0)),Bmin,Bstep)
        #HighFieldBVals = np.arange(Bmax, Bmax+1, 0.002)
        for i in zerofieldBVals:
            self.coil = np.array([0,0,0])
            #self.RotHamil(theta,phi)
            self.ZaxisHamil()
            self.MkEigensystem()
            Evals.append(self.Evals)
            Evecs.append(self.Evecs)
        BSample = range(Bmin,Bmax,Bstep)
        PlotB=range(-Bstep*len(zerofieldBVals), Bmax, Bstep)
        for i in BSample:
            self.coil = np.array([0,0,i])
            #self.RotHamil(theta=theta,phi=phi, chi=chi)
            self.ZaxisHamil()
            self.MkEigensystem()
            Evals.append(self.Evals)
            Evecs.append(self.Evecs)
        fig  = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)
        ax1.set_ylabel("Energy (MHz)")
        ax1.set_xlabel("B Field (mT) up to Then Pulse Sequence")
        ax1.set_title("Energy Splitting of Ideal Molecule for Super Dense Coding")
        ax1.plot(PlotB, Evals)
        ax2.set_ylabel("Mixing")
        ax2.set_xlabel("B Field (mT) up to Then Pulse Sequence")
        #Colour list = to dimension of Evecs
        colors = ['r', 'g', 'b', 'y']
        Evecs = np.array(Evecs)
        #Extract
        A,B,C,D = Evecs[:,3,0], Evecs[:,3,1], Evecs[:,3,2], Evecs[:,3,3]
        A,B,C,D = np.abs(A), np.abs(B), np.abs(C), np.abs(D)
        
        #Plot
        ax2.plot(PlotB, A, color = 'r', label = '|+1/2, +1/2>')
        ax2.plot(PlotB, B, color = 'g', label = '|+1/2, -1/2>')
        ax2.plot(PlotB, C, color = 'b', label = '|-1/2, +1/2>')
        ax2.plot(PlotB, D, color = 'y', label = '|-1/2, -1/2>') 
        ax2.legend()
        plt.show()
        """
        plt.subplot(2,1,1)
        plt.ylabel("Energy (MHz)")
        plt.xlabel("B Field (mT) up to Then Pulse Sequence")
        plt.title("Energy Splitting of Ideal Molecule for Super Dense Coding")
        plt.plot(PlotB, Evals)
        #plt.subplot(2,1,2)
        plt.ylabel("Mixing")
        plt.xlabel("B Field (mT) up to Then Pulse Sequence")
        plt.plot(PlotB, Evecs)
        plt.show()
        """
        
    def FreePlotNegative(self, Bmin,Bmax,Bstep,theta=0, phi=0, chi=0):
        Evals = []
        Evecs = []
        #HighFieldBVals = np.arange(Bmax, Bmax+1, 0.002)
        BSample = range(Bmin,Bmax,Bstep)
        for i in BSample:
            self.coil = np.array([0,0,-i])
            #self.RotHamil(theta,phi)
            self.ZaxisHamil()
            self.MkEigensystem()
            Evals.append(self.Evals)
            Evecs.append(self.Evecs)
        PlotB=range(-Bmax, Bmax, Bstep)
        for i in BSample:
            self.coil = np.array([0,0,i])
            #self.RotHamil(theta=theta,phi=phi, chi=chi)
            self.ZaxisHamil()
            self.MkEigensystem()
            Evals.append(self.Evals)
            Evecs.append(self.Evecs)
        fig  = plt.figure()
        ax1 = fig.add_subplot(2,1,1)
        ax2 = fig.add_subplot(2,1,2)
        ax1.set_ylabel("Energy (MHz)")
        ax1.set_xlabel("B Field (mT) up to Then Pulse Sequence")
        ax1.set_title("Energy Splitting of Ideal Molecule for Super Dense Coding")
        ax1.plot(PlotB, Evals)
        ax2.set_ylabel("Mixing")
        ax2.set_xlabel("B Field (mT) up to Then Pulse Sequence")
        #Colour list = to dimension of Evecs
        colors = ['r', 'g', 'b', 'y']
        Evecs = np.array(Evecs)
        #Extract
        A,B,C,D = Evecs[:,3,0], Evecs[:,3,1], Evecs[:,3,2], Evecs[:,3,3]
        E,F,G,H = Evecs[:,2,0], Evecs[:,2,1], Evecs[:,2,2], Evecs[:,2,3]
        I,J,K,L = Evecs[:,1,0], Evecs[:,1,1], Evecs[:,1,2], Evecs[:,1,3]
        M,N,O,P = Evecs[:,0,0], Evecs[:,0,1], Evecs[:,0,2], Evecs[:,0,3]
        #Plot
        
        ax2.plot(PlotB, A, color = 'r', label = '|+1/2, +1/2>', linestyle = ':')
        ax2.plot(PlotB, B, color = 'g', label = '|+1/2, -1/2>', linestyle = ':')
        ax2.plot(PlotB, C, color = 'b', label = '|-1/2, +1/2>', linestyle = ':')
        ax2.plot(PlotB, D, color = 'y', label = '|-1/2, -1/2>', linestyle = ':')
        """
        ax2.plot(PlotB, E, color = 'r', label = '|+1/2, +1/2>', linestyle = '--')
        ax2.plot(PlotB, F, color = 'g', label = '|+1/2, -1/2>', linestyle = '--')
        ax2.plot(PlotB, G, color = 'b', label = '|-1/2, +1/2>', linestyle = '--')
        ax2.plot(PlotB, H, color = 'y', label = '|-1/2, -1/2>', linestyle = '--')
        
        ax2.plot(PlotB, I, color = 'r', label = '|+1/2, +1/2>', linestyle = '--')
        ax2.plot(PlotB, J, color = 'g', label = '|+1/2, -1/2>', linestyle = '--')
        ax2.plot(PlotB, K, color = 'b', label = '|-1/2, +1/2>', linestyle = '--')
        ax2.plot(PlotB, L, color = 'y', label = '|-1/2, -1/2>', linestyle = '--')
        
        """
        ax2.plot(PlotB, M, color = 'r', label = '|+1/2, +1/2>', linestyle = '--')
        ax2.plot(PlotB, N, color = 'g', label = '|+1/2, -1/2>', linestyle = '--')
        ax2.plot(PlotB, O, color = 'b', label = '|-1/2, +1/2>', linestyle = '--')
        ax2.plot(PlotB, P, color = 'y', label = '|-1/2, -1/2>', linestyle = '--')
        
        ax2.legend()
        
        
        plt.show()

    #def RotHamil(self, theta, phi, chi):
        #Rotate = self.yaw_pitch_roll_matrix(theta, phi, chi)
        #Gs = (Rotate.dot(self.gs)).dot(np.transpose(Rotate))
        #Gi = (Rotate.dot(self.gi)).dot(np.transpose(Rotate))
        #A = (Rotate.dot(self.A)).dot(np.transpose(Rotate))
        #self.Zeeman = np.kron(np.einsum('iab,ji,j -> ab', self.SpinOp, Gs, self.coil),np.eye(2))*MHzmT + np.kron(np.eye(2),np.einsum('iab,ji,j -> ab', self.NucSpinOp, Gi, self.coil))*MHzmT/1836 # Convert from mT to MHz
        #self.Hyperfine = np.einsum('ij,jiab->ab', A, self.SpinKron) *MHzmT # Convert from mT to MHz
        #B = np.dot(Rotate, self.coil)
        #self.BVec = Rotate.dot(self.coil)
        #self.BVec = Rotate.dot(self.BVec)
        #return self.Hyperfine, self.Zeeman
    #Depreciated. Use FibbonachiRotateHamil instead.
#H = HyperH0(A = np.diag([0.05, 0.05, 0.1]), coil = 100, SpinOp = SpinOp2, NucSpinOp=SpinOp2, gi=-0.2748308, gs=2.00232)
#H.SplittingPlot(0, 1000, 5)