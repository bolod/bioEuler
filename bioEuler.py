from numpy import *



# Function to calculate the eigen vectors
# of an euler tensor in input, represented by an
# upper triangular 3x3 matrix. Returns a vector array.
def eigenvectors(E_a):
    w, v = linalg.eig(E_a)
    return v               



# Function to transform an array of vectors
# by the tensor R, it receives in input an
# array of vectors. Returns a new array
# composed of transformed vectors.
# type(atoms_xyz) = <type 'numpy.ndarray'>
def rotate(R, atoms_xyz):
    [dot(R, atom_xyz) for atom_xyz in atoms_xyz]


# Restituisce la versione matriciale
# del Tensore di Eulero per un singolo
# atomo.
# input:
#       - xyz: lista di coordinate dell'atomo.
#         type(xyz) = <type 'list'>
def matrice_E_a(xyz):
    xy = xyz[0] * xyz[1]
    xz = xyz[0] * xyz[2]
    yz = xyz[1] * xyz[2]

    E_a = matrix([xyz[0]**2, xy, xz],[xy, xyz[1]**2, yz], [xz, yz, xyz[2]**2])


# Restituisce la somma dei Tensori di Eulero di tutti gli
# atomi della proteina.
# input:
#       - atoms_xyz: array di terne di coordinate
#         type(atoms_xyz) = <type 'numpy.ndarray'>
def matrice_E(atoms_xyz):
#    E = sum([matrice_E_a(a[:-1]) for a in atoms])
    E = sum([matrice_E_a(a) for a in atoms_xyz])


# Normalizza un vettore
# input:
#       - v vettore n-dimensionale
#         type(v) = array_like
def normalizza(v):
    return v / linalg.norm(v)


