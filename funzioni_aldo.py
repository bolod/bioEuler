from numpy import *

## versione iniziale
def normalize1(v):
    length = sqrt(power(v[0],2)+power(v[1],2)+power(v[2],2))
    v = [v[0]/length, v[1]/length,v[2]/length]
    return v

## versione compatta usando un array
def normalize2(v):
    a = array(v)
    return a / sqrt(sum(a**2))

## versione con utilizzo di un modulo di numpy
def normalize3(v):
    return v / linalg.norm(v)

## test di correttezza delle tre implementazioni
print normalize1([1,2,3])
print normalize2([1,2,3])
print normalize3([1,2,3])
        

