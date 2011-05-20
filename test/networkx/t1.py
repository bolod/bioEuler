import networkx as nx
import matplotlib.pyplot as plt
from numpy import sqrt, log2


def getPathFrequencies(g, label_name='label'):
	apsp = nx.all_pairs_shortest_path(g)

	l = [v0 for v in apsp.itervalues() for v0 in v.itervalues()]
	l_label = [[g.node[i][label_name] for i in item]  for item in l]


	prepared_l_label0 = [zip(path, ['/' for i in range(len(path))]) for path in l_label]
	prepared_l_label1 = [[el for subsublist in sublist for el in subsublist][:-1] for sublist in prepared_l_label0]

	all_paths = ['//' + ''.join(str_list) for str_list in prepared_l_label1]
	all_paths.sort()


	result = [[all_paths[0], 1]]

	for i in range(1, len(all_paths)):
		if all_paths[i] == all_paths[i-1]:
			result[-1][1] += 1
		else:
			result.append([all_paths[i], 1])

	freq = {}
	for el in result:
		freq[el[0]] = el[1]

	return freq


def mean(val1, val2):
	return sqrt(val1 * val2)

def mergeCargs(a,b):
	out = {}
	for k,v in a.iteritems():
		out[k] = v
	for k,v in b.iteritems():
		if out.has_key(k):
			out[k] = out[k] + v
		else:
			out[k] = v
	return out

def C(x, y=None):
	b = 2.
	log = log2

	if y:
		return C(mergeCargs(x, y))
	else:
		fps = [fp for fp in x.itervalues()]
		l = sum(fps)
		return b**(-1 * sum(map(lambda fp: fp*log(fp), [fp/float(l) for fp in fps])))

def D(a,b):
	"""
	a, b: strutture simili alla struttura freq attuale

	"""
	return (C(a,b)/mean(C(a), C(b))) - 1


if __name__ == '__main__':
	g1 = nx.DiGraph()
	g1.add_node('0', label='family')
	g1.add_node('01', label='surname')
	g1.add_node('02', label='person')
	g1.add_node('03', label='person')
	g1.add_node('020', label='name')
	g1.add_node('021', label='age')
	g1.add_node('030', label='name')
	g1.add_node('031', label='age')
	g1.add_node('032', label='shoeSize')

	# g1.add_node('0320', label='length')


	g2 = nx.DiGraph()
	g2.add_node('0', label='family')
	g2.add_node('01', label='surname')
	g2.add_node('02', label='person')
	g2.add_node('03', label='person')
	g2.add_node('020', label='name')
	g2.add_node('021', label='age')
	g2.add_node('030', label='name')
	g2.add_node('031', label='age')
	g2.add_node('032', label='shoeSize')


	g1.add_edges_from([('0','01'),('0','02'),('0','03'),('02','020'),('02','021'),('03','030'),('03','031'),('03','032')])
	g2.add_edges_from([('0','01'),('0','02'),('0','03'),('02','020'),('02','021'),('03','030'),('03','031'),('03','032')])	


	# g1.nodes(data=True)

	# nx.draw(g)
	# plt.show()


	pf1 = getPathFrequencies(g1)
	pf2 = getPathFrequencies(g2)

	print D(pf1, pf2)

	# print D(pf1, pf2)





