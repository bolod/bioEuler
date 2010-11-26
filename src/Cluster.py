from Clusterizer import *



class QTCluster():

    def __init__(self, vector_list, tuple_list):
        self.tuple_list = tuple_list
        self.vector_list = vector_list
        self.dictionary = {}
        for v in vector_list:
            self.dictionary.update({self.recognize(tuple_list, v):v})


    def __repr__(self):
        s = ""
        s += "------------------\n"
        s += "Cluster:\n"
        s += str(self.dictionary.keys())
        s += "\n------------------"
        return s

    def recognize(self, tuple_list, val):
        for i in tuple_list:
            if i[1]== val:
                return i[0]