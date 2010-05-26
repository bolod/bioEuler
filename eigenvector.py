#from numpy import array
import numpy as np
from numpy import linalg as LA

# Function to calculate the eigen vectors
# of an euler tensor in input, represented by an
# upper triangular 3x3 matrix. Returns a vector array.
def eigens(tensorMatrix):
    w, v = LA.eig(tensorMatrix)
    return v               

# Function to transform an array of vectors
# by the tensor R, it receives in input an
# array of vectors. Returns a new array
# composed of transformed vectors.
def rotate(R, vectors):
    result = [np.dot(R, vector) for vector in vectors]
    return result

# Not tested, use "rotate".
def rotate2(R, vectors):
    result = np.dot(R, vectors)
    return result
                   
                   
    
