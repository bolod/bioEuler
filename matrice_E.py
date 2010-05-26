from numpy import *


def matrice_E_a(xyz):
    xy = xyz[0] * xyz[1]
    xz = xyz[0] * xyz[2]
    yz = xyz[1] * xyz[2]

    E_a = matrix([xyz[0]**2, xy, xz],[xy, xyz**2, yz], [xz, yz, xyz[2]**2])


def matrice_E(atoms):
    E = sum([matrice_E_a(a[:-1]) for a in atoms])
