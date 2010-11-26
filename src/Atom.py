from numpy import *
from pyplasm import *
from AtomInfos import *

class Atom():

    def __init__(self, symbol, label, coords=array([0.,0.,0.]), res_number=0):
        
        def symbol_to_mass(atom): return atomic_radius[atom][5]
        def symbol_to_radius(atom): return atomic_radius[atom][1]

        self.symbol = symbol
        self.label = label
        self.coords = array(coords)
        self.res_number = res_number
        self.m = symbol_to_mass(symbol)
        self.radius = symbol_to_radius(symbol)

    def __repr__(self):
        """
        Gets the info of this atom.

        Returns
        -------
        info : String
            info of this atom

        """
        s = '{\n'
        s += 'symbol: \'' + self.symbol + '\',\n'
        s += 'label: \'' + self.label + '\',\n'
        s += 'coords: ' + str(self.coords) + ',\n'
        s += 'mass: ' + str(self.m) + ',\n'
        s += 'radius: ' + str(self.radius) + '\n'
        s += '}\n'
        return s

    def clone(self):
        return Atom(self.symbol, self.label, self.coords[:])

    def get_coords(self):
        """
        Gets the coordinates of this atom.

        Returns
        -------
        coords : ndarray, shape (3, )
            coords of this atom

        """
        return self.coords

    def get_x(self):
        """
        Gets the x coordinate of this atom.

        Returns
        -------
        x : float
            x coord of this atom

        """
        return self.coords[0]

    def get_y(self):
        """
        Gets the y coordinate of this atom.

        Returns
        -------
        y : float
            y coord of this atom

        """
        return self.coords[1]

    def get_z(self):
        """
        Gets the z coordinate of this atom.

        Returns
        -------
        z : float
            z coord of this atom

        """
        return self.coords[2]

    def get_mass(self):
        """
        Gets the mass of this atom.

        Returns
        -------
        mass : float
            mass of this atom

        """
        return self.m

    def get_radius(self):
        """
        Gets the x coordinate of this atom.

        Returns
        -------
        x : float
            x coord of this atom

        """
        return self.radius

    def get_res_number(self):
        return self.res_number

    def get_euler(self):
        """
        Gets the euler tensor of this atom.

        Returns
        -------
        euler : ndarray, shape (3, 3)
            euler tensor of this atom

            [ x*x, x*y, x*z,
              y*x, y*y, y*z,
              z*x, z*y, z*z ]

        """
        return array([ coord * self.coords for coord in self.coords ])

    def rotate(self, rotation):
        """
        Rotates this atom by the given rotation matrix.

        Parameters
        ----------
        rotation: ndarray, shape(3, 3)
            the rotation matrix

        Returns
        -------
        atom : Atom
            this atom rotated

        """
        self.coords = dot(rotation, self.coords)
        return self

    def translate(self, translation):
        """
        Translates this atom by the given translation vector.

        Parameters
        ----------
        translation: ndarray, shape(3, )
            the translation vector

        Returns
        -------
        atom : Atom
            this atom translated

        """
        self.coords = self.coords - translation
        return self

    def scale(self, scale):
        """
        Scales this atom by the given scale value vector.

        Parameters
        ----------
        scale: ndarray, shape(3, )
            the scale value vector

        Returns
        -------
        atom : Atom
            this atom scaled

        """
        self.coords = self.coords * scale
        return self

    def get_symbol(self):
        """
        Gets the symbol of this atom.

        Returns
        -------
        symbol : String
            symbol of this atom

        """
        return self.symbol

    def get_label(self):
        """
        Gets the label of this atom.

        Returns
        -------
        label : String
            symbol of this atom

        """
        return self.label

    def get_ammino_chain_seq(self):
        """
        Gets the chain sequence number of the amminoacid this atom belongs to.

        Returns
        -------
        seq : int
            the chain sequence number of the amminoacid this atom belongs to

        """
        return self.ammino_chain_seq

    def get_euclid_distance_to(self, atom):
        """
        Gets the euclid distance from this atom to the given atom.

        Returns
        -------
        distance : float
            the distance from this atom to the diven atom

        """
        return linalg.norm(self.get_coords() - atom.get_coords())

    def plasm_cube(self, size=0.1, color=WHITE):
        """
        Create a PLaSM cuboid with a color an put it on this atom coords.

        Parameters
        ----------
        size: float
            size of the cube
        color : integer
            number of the color in the dictionary "colors"

        Returns
        -------
        cuboid : colored PLaSM cuboid traslated to this atom coords

        """
        return COLOR(color)(T([1,2,3])(self.coords)(CUBOID([size, size, size])))