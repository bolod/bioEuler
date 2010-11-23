from pyplasm import *
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *

class Plasm():

    def create_cube(self, coords, size=0.1, color=WHITE):
        """
        Create a PLaSM cuboid with given size and coords.

        Parameters
        ----------
        size: float
            size of the cube

        Returns
        -------
        cuboid : PLaSM cuboid
            PLaSM cuboid with given size and coords.

        """
        return (T([1,2,3])(coords)(CUBOID([size, size, size])))

    def create_cube_list(self, coords_list, size=0.1, color=WHITE):
        """
        Create a PLaSM cuboid's list

        Parameters
        ----------
        size : float
            size of each cube
        color : integer
            number of the color in the dictionary "colors"

        Returns
        -------
        cuboid_list : colored PLaSM cuboid's list

        """
        return [ self.create_cube(coords, size, color) for coords in coords_list ]

    def create_ellipsoid(self, coords, size=0.1, color=WHITE):
        """
        Create a PLaSM ellipsoid of this atom list.

        Parameters
        ----------
        color : integer
            number of the color in the dictionary "colors"

        Returns
        -------
        ellipsoid : colored PLaSM ellipsoid

        """
        return COLOR(color)(T([1,2,3])(coords)((S([1,2,3])(size)(SPHERE(1)([16,16])))))

    def create_convex_hull(self, coords_list, size=0.1, color=WHITE):
        """
        Create a PLaSM convex hull of this atom list.

        Parameters
        ----------
        color : integer
            number of the color in the dictionary "colors"

        Returns
        -------
        convex_hull : colored PLaSM ellipsoid

        """
        return COLOR(color)(JOIN(self.create_cube_list(coords_list, size)))

    def create_polyline(self, coords_list):
        """
        Create a PLaSM polyline of this atom list.

        Parameters
        ----------
        coords_list : array_like, shape(dim, 3)
            coords list
        Returns
        -------
            polyline
        """
        return POLYLINE(coords_list)

    def color(self, plasm_object, color=WHITE):
        """
        Add color to a PLaSM objects.

        Parameters
        ----------
        plasm_object : PLaSM object
            PLaSM object to color
        color : PLaSM color
            PLaSM color
        """
        colors = {
            '1': WHITE,
            '2': RED,
            '3': GREEN,
            '4': BLUE }

        return COLOR(color)(plasm_object)

    def view(self, plasm_object_list):
        """
        View a list of PLaSM objects.

        Parameters
        ----------
        plasm_object_list : PLaSM object list
            list of PLaSM objects
        """
        VIEW(STRUCT(plasm_object_list))