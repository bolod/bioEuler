from numpy import array
from numpy import linalg as LA
import numpy as np

def eigenvectors(eulerTensor):
    a, b, c, d, e, f = eulerTensor
    matrix = array([[a, d, f],[0, b, e],[0, 0, c]])
    w, v = LA.eig(matrix)
    return v               


                   
                   
    
