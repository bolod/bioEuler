from BioTree import *
from Decoder import *
import os

class Imprint:
    def __init__(self, bio_tree):
        self.decoder = Decoder()
        self.bio_tree = bio_tree
        self.bio_vector = bio_tree.vectorialize()
        self.encoded_bio_vector = self.decoder.encode_vector(self.bio_vector)

    def get_imprint_vector(self):
        return self.bio_vector

    def get_encoded_imprint_vector(self):
        return self.encoded_bio_vector

    def write_on_file(self, filename):
        """
        Print imprint on file.

        Parameters
        ----------
        self: Imprint

        """
        path_name = "./impronte/"
        dir = os.path.dirname(path_name)
        try:
            os.stat(dir)
        except:
            os.mkdir(dir)
        FILE = open(path_name + filename,"w")
        FILE.write(self.encoded_bio_vector)
        FILE.close()






def read_from_file(file_name, decoder=Decoder()):
    path = "./impronte/"
    file = open(path+file_name,"r")
    impr = file.read()
    file.close()
    return decoder.decode(impr)

def create_vectors_from_file(filename):
    temp = []
    file = open(filename, "r")
    iter(file)

    for line in file:
        l = line.strip('\n\t')
        vector = read_from_file(l)
        couple = (l, vector)
        temp.append(couple)
    file.close()
    return temp