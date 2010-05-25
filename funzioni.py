from numpy import *

def normalize(v):
    length = sqrt(power(v[0],2)+power(v[1],2)+power(v[2],2))
    v = [v[0]/length, v[1]/length,v[2]/length]
    return v

print normalize([1,2,3])