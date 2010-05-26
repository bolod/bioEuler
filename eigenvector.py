#from numpy import array
from numpy import linalg as LA
import numpy as np

# Function to calculate the eigen vectors
# of an euler tensor in input, represented by an
# upper triangular 3x3 matrix, specified by
# 6 elements.
# Example:
# Given the A matrix as:
#   a d f
#   0 b e
#   0 0 c
# it's represented like a vector of 6 values:
# [a, b, c, d, e, f]
def eigenvectors(tensorMatrix):
    #a, b, c, d, e, f = eulerTensor
    #matrix = array([[a, d, f],[0, b, e],[0, 0, c]])
    w, v = LA.eig(tensorMatrix)
    return v               

# Function to transform an array of vectors
# by the tensor R, it receives in input an
# array of vectors and produces a new array
# composed of transformed vectors.
def rotate(R, vectors):
    result = [np.dot(R, vector) for vector in vectors]
    return result

# Not tested, use "rotate"
def rotate2(R, vectors):
    result = np.dot(R, vectors)
    return result
                   
                   
    
