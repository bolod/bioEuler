from numpy import *


def qt_clustering(G,d):

    while(len(G) > 0):
        if (len(G) == 1):  #caso base
            yield G
            return
        elif (len(G) == 0):
            return
        else:
            A_list = []
            for i in G:
                flag = True
                Ai = [i]
                while flag and not([x for x in G if x in Ai] == G):
                    app_j = [x for x in G if x not in Ai]
                    diameters = []
                    app_A = []
                    for x in app_j:
                        app_A = Ai[:]
                        app_A.append(x)
                        diameters.append(get_diameter(app_A))

                    min_indx = argmin(diameters)
                    j = app_j[min_indx]
                    diameter = diameters[min_indx]
                    if(diameter > d):
                        flag = False
                    else:
                        Ai.append(j)
                A_list.append(Ai)

            C = A_list[argmax([len(x) for x in A_list])]
            yield C
            G = [x for x in G if x not in C]
    return


def get_centroid(l):
    return sum([x for x in l], axis=0) / len(l)

def get_diameter(l):
    centroid = get_centroid(l)
    return max([get_weighted_distance(centroid, x) for x in l]) * 2

def get_distance(v1, v2):
    return linalg.norm(array(v2) - array(v1))

def get_weighted_distance(v1, v2):
    X = len(v1)/4
    return sqrt(sum([(v1[i] - v2[i])**2 * (X/int(log2(floor(i/4)+1)+1)) for i in range(len(v1))]))



def dictionary_clustering(input, R):
        dic = {}
        for i in range(len(input)):
            for j in range(i+1, len(input)):
                d = int(get_weighted_distance(input[i][1],input[j][1]))
                if (d <= R):
                    if input[i][0] in dic:
                        dic[input[i][0]].append((input[j][0],d))
                    else:
                        dic.update({input[i][0]:[(input[j][0],d)]})

                    if input[j][0] in dic:
                        dic[input[j][0]].append((input[i][0],d))
                    else:
                        dic.update({input[j][0]:[(input[i][0],d)]})
        return dic

