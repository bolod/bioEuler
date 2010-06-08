def z_split(atom_list):
	""" 
	splits a list of atoms along theirs z-axis
	"""
	filter(lambda atom: atom[2] <= 0, atom_list), 
	filter(lambda atom: atom[2] > 0, atom_list)