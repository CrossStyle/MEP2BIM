import math
import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util
from ifcopenshell.util.placement import get_local_placement
import numpy as np


class Patcher:
    def __init__(self, file, args=None):
        self.file = file
        self.args = args

    def patch(self, product, transformation):
        absolute_placement = self.get_absolute_placement(product.ObjectPlacement)
        a = transformation @ ifcopenshell.util.placement.get_local_placement(absolute_placement)
        absolute_placement.RelativePlacement = self.get_relative_placement(a)

    def get_absolute_placement(self, object_placement):
        if object_placement.PlacementRelTo:
            return self.get_absolute_placement(object_placement.PlacementRelTo)
        return object_placement

    def z_rotation_matrix(self, angle):
        return [
            [math.cos(angle), -math.sin(angle), 0., 0.],
            [math.sin(angle), math.cos(angle), 0., 0.],
            [0., 0., 1., 0.],
            [0., 0., 0., 1.],
        ]

    def get_relative_placement(self, m):
        x = np.array((m[0][0], m[1][0], m[2][0]))
        z = np.array((m[0][2], m[1][2], m[2][2]))
        o = np.array((m[0][3], m[1][3], m[2][3]))
        object_matrix = ifcopenshell.util.placement.a2p(o, z, x)
        return self.create_ifc_axis_2_placement_3d(
            object_matrix[:, 3][0:3],
            object_matrix[:, 2][0:3],
            object_matrix[:, 0][0:3],
        )

    def create_ifc_axis_2_placement_3d(self, point, up, forward):
        return self.file.createIfcAxis2Placement3D(
            self.file.createIfcCartesianPoint(point.tolist()),
            self.file.createIfcDirection(up.tolist()),
            self.file.createIfcDirection(forward.tolist()),
        )