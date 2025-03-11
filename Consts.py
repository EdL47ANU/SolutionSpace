from scipy.constants import hbar, Planck, c
import numpy as np
muB = 9.274009994*10**-24 # J/T
MHzmT = (muB)/(Planck*(10**9))
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#…Protons, and neutrons is the nuclear magneton (equivalent to 5.05078 × 10−27 joule per tesla).
#The nuclear magneton, calculated by using the mass of the proton (rather than that of the electron, used to calculate the Bohr magneton) equals 1/1,836 Bohr magneton. See magnetic dipole.

Sx2 = np.array([ [0,    1],
                 [1,    0]])
Sy2 = np.array([ [0,    -1j],
                 [1j,   0]])
Sz2 = np.array([ [1,    0],
                 [0,    -1]])
SpinOp2 = np.array([Sx2, Sy2, Sz2])
Sx3 = np.array([[0,    1,   0],
                [1,    0,   1],
                [0,    1,   0]])
Sy3 = np.array([[0,    -1j,   0],
                [1j,    0,   -1j],
                [0,    1j,   0]])
Sz3 = np.array([[1,    0,   0],
                [0,    0,   0],
                [0,    0,   -1]])
SpinOp3 = np.array([Sx3, Sy3, Sz3])
Sx4 = np.array([[0,np.sqrt(3),  0,0],
                [np.sqrt(3),0,  2,0],
                [0,         2,0,np.sqrt(3)],
                [0,         0,np.sqrt(3),0]])
Sy4 = np.array([[0,-np.sqrt(3)*1j, 0,   0],
                [np.sqrt(3)*1j,0,   -2j,0],
                [0,         2j,0, -np.sqrt(3)*1j],
                [0,         0,np.sqrt(3)*1j, 0]])
Sz4 = np.array([[3,0,0,0],
                [0,1,0,0],
                [0,0,-1,0],
                [0,0,0,-3]])
SpinOp4 = np.array([Sx4, Sy4, Sz4])

def DEcm_ZfsMHZ(D, E):
    D,E = D * 29979, E * 29979# Convert from cm^-1 to MHz
    DMat = np.diag([(E-(D/3)), -E -(D/3), 2*D/3])  # Diagonal elements of D matrix
    # print("D Matrix = \n" ,DMatrix)
    return DMat
def PolarToCartesian(theta, phi, r):
    return np.sin(theta)*np.cos(phi)*r, np.sin(theta)*np.sin(phi)*r, np.cos(theta)*r
#thetas, phis, chis = GoldenSpiralAnglesRoll(100)

def FibbonaciSphere(samples):
    points = []
    phi = (1 + np.sqrt(5)) / 2  # golden ratio
    for i in range(samples):
        theta = 2 * np.pi * (i / phi - np.floor(i / phi))
        z = 1 - (2 * i) / (samples - 1)
        r = np.sqrt(1 - z * z)
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        points.append((x, y, z))
    points = np.array(points)
    # Check if all points are unitary
#    norms = np.linalg.norm(points, axis=1)
#    if not np.allclose(norms, 1):
#        print("Warning: Not all points are unitary")
#    else:
#        print("All points are unitary")
    return points

def OrthoVecs(OrientationVecs):
    Ortho1Vecs = []
    Ortho2Vecs = []
    if np.shape(OrientationVecs) == (3,):
        if np.allclose(OrientationVecs, [1, 0, 0]):
            temp_vec = np.array([0, 1, 0], dtype = float)
        else:
            temp_vec = np.array([1, 0, 0], dtype = float)
        # First orthogonal vector   
        orthogonal_vec1 = np.cross(OrientationVecs, temp_vec)
        orthogonal_vec1 /= np.linalg.norm(orthogonal_vec1)
        # Second orthogonal vector
        orthogonal_vec2 = np.cross(OrientationVecs, orthogonal_vec1)
        orthogonal_vec2 /= np.linalg.norm(orthogonal_vec2)
        return orthogonal_vec1, orthogonal_vec2
    else:
        for OrientationVec in OrientationVecs:
            if np.allclose(OrientationVec, [1, 0, 0]):
                temp_vec = np.array([0, 1, 0], dtype = float)
            else:
                temp_vec = np.array([1, 0, 0], dtype = float)
            # First orthogonal vector
            orthogonal_vec1 = np.cross(OrientationVec, temp_vec)
            orthogonal_vec1 /= np.linalg.norm(orthogonal_vec1)
            # Second orthogonal vector
            orthogonal_vec2 = np.cross(OrientationVec, orthogonal_vec1)
            orthogonal_vec2 /= np.linalg.norm(orthogonal_vec2)
            Ortho1Vecs.append(orthogonal_vec1)
            Ortho2Vecs.append(orthogonal_vec2)
        return np.array(Ortho1Vecs), np.array(Ortho2Vecs)

def Orientations(NOrient):
    OrientationVecs = FibbonaciSphere(NOrient)
    Ortho1Vecs, Ortho2Vecs = OrthoVecs(OrientationVecs)
    return list(zip(OrientationVecs, Ortho1Vecs, Ortho2Vecs))

def plot_fibonacci_sphere(NOrients):
    # Generate the vectors using FibonacciSphere
    vectors = FibbonaciSphere(NOrients)
    
    # Create a figure and a 3D axis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Create a colormap
    cmap = plt.get_cmap('viridis')
    
    # Normalize the range of colors
    norm = plt.Normalize(0, NOrients)
    
    # Plot the vectors as a continuous line with changing color
    for i in range(NOrients - 1):
        ax.plot(vectors[i:i+2, 0], vectors[i:i+2, 1], vectors[i:i+2, 2], color=cmap(norm(i)))
    
    # Show the plot
    plt.show()


def unit_vector_to_spherical_coords(unit_vector):
    """
    Convert a unit vector to spherical coordinates (theta, phi).
    Parameters:
    unit_vector (array-like): A unit vector [x, y, z].
    Returns:
    tuple: (theta, phi) where theta is the angle from the z-axis and phi is the angle from the x-axis in the xy-plane.
    """
    x, y, z = unit_vector
    # Calculate theta
    theta = np.arccos(z)
    # Calculate phi
    phi = np.arctan2(y, x)
    return theta, phi

def spherical_coords_to_unit_vector(theta, phi):
    """
    Convert spherical coordinates (theta, phi) to a unit vector.
    Parameters:
    theta (float): The angle from the z-axis.
    phi (float): The angle from the x-axis in the xy-plane.
    Returns:
    array-like: A unit vector [x, y, z].
    """
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    return np.array([x, y, z])

theta, phi = unit_vector_to_spherical_coords([0.74888958, 0.05988857, 0.65998315])
print(theta, phi)
# Example usage
#plot_fibonacci_sphere(50)


#xyz = FibbonaciSphere(100)
#abc = FibbonaciSphere(50)
#ax = plt.figure().add_subplot(projection='3d')
#xyz = xyz[:len(xyz)//2]
#ax.scatter(xyz[:,0], xyz[:,1], xyz[:,2], color = 'pink')
#ax.scatter(abc[:,0], abc[:,1], abc[:,2], color = 'blue')
#plt.show()

# Faith in Einsum restored.
#print(np.einsum('iab,ij,jbc -> ac', 0.5*SpinOp4,np.diag([1,1,1]), 0.5*SpinOp4))
#print(DEcm_ZfsMHZ(-0.23, 0.053))

def rotate_vectors(axis, vec1, vec2, angle):
    # Normalize the axis
    axis = axis / np.linalg.norm(axis)
    # Compute the rotation matrix using Rodrigues' rotation formula
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)
    # Rotate the vectors
    rotated_vec1 = np.dot(R, vec1)
    rotated_vec2 = np.dot(R, vec2)
    return rotated_vec1, rotated_vec2

def OrthogCircle(axis, NOrient):
    Z = axis/np.linalg.norm(axis)
    X,Y = OrthoVecs(axis)
    angles = np.linspace(0, 2*np.pi, NOrient)
    New_x, New_y = [],[]
    for angle in angles:
        new_x, new_y = rotate_vectors(Z, X, Y, angle)
        New_x.append(new_x)
        New_y.append(new_y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in New_x:
        ax.scatter(i[0], i[1], i[2], c = 'r')
    for i in New_y:
        ax.scatter(i[0], i[1], i[2], c = 'b')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    plt.show()




#OrthogCircle(np.array([1,0,0]), 100)
#OrthogCircle(np.array([np.sqrt(2)/2,np.sqrt(2)/2,0]), 100)