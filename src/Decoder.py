from BioTree import *

class Decoder():
    def __init__(self, n_digit=5, precision=float(10**5)):
        self.n_digit = n_digit
        self.base = 36
        self.precision = precision
    
    def encode_num_base10_to_baseN(self, num, n):
        """
        Change a to a base-n number.
        Up to base-36 is supported without special notation.
        """
        new_num_string = ''
        current = num
        while current != 0:
            remainder = current % n
            if 36 > remainder > 9:
                remainder_string = chr(ord('a') + remainder - 10)
            elif remainder >= 36:
                remainder_string = '(' + str(remainder) + ')'
            else:
                remainder_string = str(remainder)
            new_num_string = remainder_string + new_num_string
            current = current / n
        return new_num_string

    def encode_vector(self, vector):
        return "".join([value.encode(self) for value in vector])

    def encode_value(self, value, prec=-1):
        if prec == -1:
            prec = self.precision
        if self.base**self.n_digit < value:
            print "Warning: Superato il valore massimo codificabile."
            enc = self.base**self.n_digit
        else:
            enc = str(self.encode_num_base10_to_baseN(int(value*prec), self.base))
        return enc.zfill(self.n_digit)

    def decode(self, value_string):
        n_digit = self.n_digit
        
        imp_sliced_l = [value_string[i*n_digit:i*n_digit+n_digit] for i in range(len(value_string)/n_digit)]

        result = []
        for i in range(0,len(imp_sliced_l),4):
            result.append(self.decode_value(imp_sliced_l[i],1))
            result.append(self.decode_value(imp_sliced_l[i+1]))
            result.append(self.decode_value(imp_sliced_l[i+2]))
            result.append(self.decode_value(imp_sliced_l[i+3]))

        return result

    def decode_value(self, value_string, precision=-1):
        if precision == -1:
            precision = self.precision
        base36 = {'0':0, '1':1, '2':2, '3':3, '4':4, '5':5,   \
              '6':6, '7':7, '8':8, '9':9, 'a':10, 'b':11,     \
              'c':12, 'd':13, 'e':14, 'f':15, 'g':16, 'h':17, \
              'i':18, 'j':19, 'k':20, 'l':21, 'm':22, 'n':23, \
              'o':24, 'p':25, 'q':26, 'r':27, 's':28, 't':29, \
              'u':30, 'v':31, 'w':32, 'x':33, 'y':34, 'z':35}
              
        return sum([ long(base36[value_string[i]])*self.base**(self.n_digit-(i+1)) for i in range(self.n_digit)])/precision
