from numpy import *
from operator import *

class BioNode():
    def __init__(self, euler, mass):
        ro_v = sqrt(diag(euler) / mass)
        self.value = [mass] + ro_v.tolist()
        self.weight = mass * sqrt(reduce(mul, ro_v) / sum(ro_v))

    def __repr__(self):
        return str(self.value)

    def get_value(self):
        return self.value

    def get_weight(self):
        return self.weight

    def encode(self, decoder):
        return str(decoder.encode_value(self.value[0], 1)) \
            + str(decoder.encode_value(self.value[1])) \
            + str(decoder.encode_value(self.value[2])) \
            + str(decoder.encode_value(self.value[3]))

if __name__ == "__main__":
    print "Hello World"
