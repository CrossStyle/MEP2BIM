import math
import uuid
import time
import tempfile
import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeWedge, BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.gp import gp_Pnt, gp_Ax2, gp_Dir, gp_Circ
from ifcopenshell.util.element import copy_deep
import numpy as np
import pickle
from patcher import Patcher


O = 0., 0., 0.
X = 1., 0., 0.
Y = 0., 1., 0.
Z = 0., 0., 1.


def create_ifcaxis2placement(ifcfile, point=O, dir1=Z, dir2=X):
    point = ifcfile.createIfcCartesianPoint(point)
    dir1 = ifcfile.createIfcDirection(dir1)
    dir2 = ifcfile.createIfcDirection(dir2)
    axis2placement = ifcfile.createIfcAxis2Placement3D(point, dir1, dir2)
    return axis2placement


def create_ifcaxis2placement_2d(ifcfile, point=O, dir1=Z):
    point = ifcfile.createIfcCartesianPoint(point)
    dir1 = ifcfile.createIfcDirection(dir1)
    axis2placement = ifcfile.createIfcAxis2Placement2D(point, dir1)
    return axis2placement


# Creates an IfcLocalPlacement from Location, Axis and RefDirection, specified as Python tuples, and relative placement
def create_ifclocalplacement(ifcfile, point=O, dir1=Z, dir2=X, relative_to=None):
    axis2placement = create_ifcaxis2placement(ifcfile, point, dir1, dir2)
    ifclocalplacement2 = ifcfile.createIfcLocalPlacement(relative_to, axis2placement)
    return ifclocalplacement2


# Creates an IfcPolyLine from a list of points, specified as Python tuples
def create_ifcpolyline(ifcfile, point_list):
    ifcpts = []
    for point in point_list:
        point = ifcfile.createIfcCartesianPoint(point)
        ifcpts.append(point)
    polyline = ifcfile.createIfcPolyLine(ifcpts)
    return polyline


# Creates an IfcExtrudedAreaSolid from a list of points, specified as Python tuples
def create_ifcextrudedareasolid(ifcfile, point_list, ifcaxis2placement, extrude_dir, extrusion):
    polyline = create_ifcpolyline(ifcfile, point_list)
    ifcclosedprofile = ifcfile.createIfcArbitraryClosedProfileDef("AREA", None, polyline)
    ifcdir = ifcfile.createIfcDirection(extrude_dir)
    ifcextrudedareasolid = ifcfile.createIfcExtrudedAreaSolid(ifcclosedprofile, ifcaxis2placement, ifcdir, extrusion)
    return ifcextrudedareasolid


def create_ifcpolyloop(ifcfile, point_list):
    ifcpts = []
    for point in point_list:
        point = ifcfile.createIfcCartesianPoint(point)
        ifcpts.append(point)
    polyline = ifcfile.createIfcPolyLoop(ifcpts)
    return polyline


create_guid = lambda: ifcopenshell.guid.compress(uuid.uuid1().hex)
# IFC template creation
filename = "demo1.ifc"
timestamp = round(time.time())
timestring = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(timestamp))
creator = "Zaolin Pan"
organization = "UST"
application, application_version = "IfcOpenShell", "0.7"
project_globalid, project_name = create_guid(), "the MEP system 1"

template = """ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('ViewDefinition [CoordinationView]'),'2;1');
FILE_NAME('%(filename)s','%(timestring)s',('%(creator)s'),('%(organization)s'),'%(application)s','%(application)s','');
FILE_SCHEMA(('IFC4'));
ENDSEC;
DATA;
#1=IFCPERSON($,$,'%(creator)s',$,$,$,$,$);
#2=IFCORGANIZATION($,'%(organization)s',$,$,$);
#3=IFCPERSONANDORGANIZATION(#1,#2,$);
#4=IFCAPPLICATION(#2,'%(application_version)s','%(application)s','');
#5=IFCOWNERHISTORY(#3,#4,$,.ADDED.,$,#3,#4,%(timestamp)s);
#6=IFCDIRECTION((1.,0.,0.));
#7=IFCDIRECTION((0.,0.,1.));
#8=IFCCARTESIANPOINT((0.,0.,0.));
#9=IFCAXIS2PLACEMENT3D(#8,#7,#6);
#10=IFCDIRECTION((0.,1.,0.));
#11=IFCGEOMETRICREPRESENTATIONCONTEXT($,'Model',3,1.E-05,#9,#10);
#12=IFCDIMENSIONALEXPONENTS(0,0,0,0,0,0,0);
#13=IFCSIUNIT(*,.LENGTHUNIT.,.MILLI.,.METRE.);
#14=IFCSIUNIT(*,.AREAUNIT.,.MILLI.,.SQUARE_METRE.);
#15=IFCSIUNIT(*,.VOLUMEUNIT.,.MILLI.,.CUBIC_METRE.);
#16=IFCSIUNIT(*,.PLANEANGLEUNIT.,$,.RADIAN.);
#17=IFCMEASUREWITHUNIT(IFCPLANEANGLEMEASURE(0.017453292519943295),#16);
#18=IFCCONVERSIONBASEDUNIT(#12,.PLANEANGLEUNIT.,'DEGREE',#17);
#19=IFCUNITASSIGNMENT((#13,#14,#15,#18));
#20=IFCPROJECT('%(project_globalid)s',#5,'%(project_name)s',$,$,$,$,(#11),#19);
ENDSEC;
END-ISO-10303-21;
""" % locals()

# Write the template to a temporary file
temp_handle, temp_filename = tempfile.mkstemp(suffix=".ifc")
with open(temp_filename, "wb") as f:
    f.write(template.encode())


def create_pipe_material(_pipes):
    material = ifcfile.createIfcMaterial("pipe material")
    material_layer = ifcfile.createIfcMaterialLayer(material, 0.2, None)
    material_layer_set = ifcfile.createIfcMaterialLayerSet([material_layer], None)
    material_layer_set_usage = ifcfile.createIfcMaterialLayerSetUsage(material_layer_set, "AXIS2", "POSITIVE", -0.1)
    ifcfile.createIfcRelAssociatesMaterial(create_guid(), owner_history, RelatedObjects=_pipes,
                                           RelatingMaterial=material_layer_set_usage)


ifcfile = ifcopenshell.open(temp_filename)
owner_history = ifcfile.by_type("IfcOwnerHistory")[0]
project = ifcfile.by_type("IfcProject")[0]
context = ifcfile.by_type("IfcGeometricRepresentationContext")[0]
site_placement = create_ifclocalplacement(ifcfile)
site = ifcfile.createIfcSite(create_guid(), owner_history, "Site", None, None, site_placement, None, None, "ELEMENT",
                             None, None, None, None, None)
building_placement = create_ifclocalplacement(ifcfile, relative_to=site_placement)
building = ifcfile.createIfcBuilding(create_guid(), owner_history, 'Building', None, None, building_placement, None,
                                     None, "ELEMENT", None, None, None)
storey_placement = create_ifclocalplacement(ifcfile, relative_to=building_placement)
elevation = 0.0
building_storey = ifcfile.createIfcBuildingStorey(create_guid(), owner_history, 'Storey', None, None,
                                                  storey_placement,
                                                  None, None, "ELEMENT", elevation)
container_storey = ifcfile.createIfcRelAggregates(create_guid(), owner_history, "Building Container", None,
                                                  building,
                                                  [building_storey])
container_site = ifcfile.createIfcRelAggregates(create_guid(), owner_history, "Site Container", None, site,
                                                [building])
container_project = ifcfile.createIfcRelAggregates(create_guid(), owner_history, "Project Container", None, project,
                                                   [site])
flow_placement = create_ifclocalplacement(ifcfile, relative_to=storey_placement)


def colour(entities, color):
    if color is not None:
        r, g, b = color[0], color[1], color[2]
        owner_history = ifcfile.by_type("IfcOwnerHistory")[0]
        context = ifcfile.by_type("IfcGeometricRepresentationContext")[0]
        material = ifcfile.createIfcMaterial("my material")
        material_layer = ifcfile.createIfcMaterialLayer(material, 0.2, None)
        material_layer_set = ifcfile.createIfcMaterialLayerSet([material_layer], None)
        material_layer_set_usage = ifcfile.createIfcMaterialLayerSetUsage(material_layer_set, "AXIS2", "POSITIVE", -0.1)
        ifcfile.createIfcRelAssociatesMaterial(ifcopenshell.guid.new(), owner_history, RelatedObjects=entities,
                                               RelatingMaterial=material_layer_set_usage)
        # 0.93, 0.75, 0.38
        rgb = ifcfile.createIfcColourRgb("my color", r, g, b)
        factor = ifcfile.createIfcNormalisedRatioMeasure(0.65)
        factor1 = ifcfile.createIfcNormalisedRatioMeasure(0.67)
        rendering = ifcfile.createIfcSurfaceStyleRendering(rgb, 0.0, factor, None, None, None, factor1, None, "NOTDEFINED")
        style = ifcfile.createIfcSurfaceStyle("my style", "BOTH", [rendering])
        assignment = ifcfile.createIfcPresentationStyleAssignment([style])
        styled_item = ifcfile.createIfcStyledItem(None, [assignment], None)
        sub_context = ifcfile.createIFCGEOMETRICREPRESENTATIONSUBCONTEXT("Body", "Model",
                                                                         None, None, None, None,
                                                                         context, None, "MODEL_VIEW", None)
        styled_representation = ifcfile.createIFCSTYLEDREPRESENTATION(sub_context, None, None, [styled_item])
        ifcfile.createIfcMaterialDefinitionRepresentation(None, None, [styled_representation], material)


def generate_ifc_polygon(type, name, _contours, _bottom_height, _wall_height, _scale, color):
    _contours = _scale * _contours[:, [1, 0]]
    storey = ifcfile.by_type('IfcBuildingStorey')[0]
    OwnerHistory = ifcfile.by_type('IfcOwnerHistory')[0]
    context = ifcfile.by_type("IfcGeometricRepresentationContext")[0]
    space_placement = create_ifclocalplacement(ifcfile)
    extrusion_placement = create_ifcaxis2placement(ifcfile, point=(0.0, 0.0, float(_bottom_height)))
    point_list_extrusion_area = _contours.tolist()
    point_list_extrusion_area.append(_contours[0, :].tolist())
    solid = create_ifcextrudedareasolid(ifcfile,
                                        point_list_extrusion_area,
                                        extrusion_placement,
                                        (0.0, 0.0, 1.0),
                                        float(_wall_height))
    body_representation = ifcfile.createIfcShapeRepresentation(context, "Body", "SweptSolid", [solid])
    entity_shape = ifcfile.createIfcProductDefinitionShape(None, None, [body_representation])
    entity = ifcfile.create_entity(type, ifcopenshell.guid.new())
    entity.OwnerHistory = OwnerHistory
    entity.Name = name
    entity.Representation = entity_shape
    entity.ObjectPlacement = space_placement
    relatingObject = storey
    related_objects = []
    related_objects.append(entity)
    ifcfile.createIfcRelAggregates(ifcopenshell.guid.new(),
                                   OwnerHistory,
                                   None,
                                   None,
                                   relatingObject,
                                   related_objects)
    colour(related_objects, color)


def create_ifc_flow_terminal(length, width, cx, cy, cz, color):
    # length, width = 41, 41
    height = 0.3*length
    my_wedge = BRepPrimAPI_MakeWedge(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), length, height, width, height, height, 0.7*length,
                                     0.7*width).Shape()
    p1 = gp_Pnt(height, height, height)
    p2 = gp_Pnt(0.7*length, height+300, 0.7*width)
    my_box = BRepPrimAPI_MakeBox(p1, p2).Shape()
    myBody = BRepAlgoAPI_Fuse(my_wedge, my_box).Shape()

    rep = ifcopenshell.geom.tesselate(ifcfile.schema, myBody, 1.)
    ifcfile.add(rep)

    context = ifcfile.by_type("IfcGeometricRepresentationContext")[0]
    rep.Representations[0].ContextOfItems = context

    OwnerHistory = ifcfile.by_type('IfcOwnerHistory')[0]
    terminal = ifcfile.create_entity('IfcFlowTerminal', ifcopenshell.guid.new())
    terminal.OwnerHistory = OwnerHistory
    terminal.Representation = rep
    terminal.ObjectPlacement = create_ifclocalplacement(ifcfile)
    terminal.Name = 'air terminal'
    new_terminal = copy_deep(ifcfile, terminal)

    angle = 90.0
    movement = np.array([
        [1.0, 0.0, 0.0, 0.5*length+cx],
        [0.0, math.cos(math.radians(angle)), -math.sin(math.radians(angle)), 0.5*width +cy],
        [0.0, math.sin(math.radians(angle)), math.cos(math.radians(angle)), -height-300.0+cz],
        [0.0, 0.0, 0.0, 1.0]])

    patcher = Patcher(ifcfile)
    patcher.patch(new_terminal, movement)
    new_terminal.GlobalId = ifcopenshell.guid.new()

    ifcfile.createIfcRelAggregates(ifcopenshell.guid.new(),
                                   OwnerHistory,
                                   None,
                                   None,
                                   ifcfile.by_type('IfcBuildingStorey')[0],
                                   [new_terminal])
    colour([new_terminal], color)


def create_ifc_rect(x1, y1, x2, y2, name, color=None):
    p1 = gp_Pnt(x1, y1, 0)
    p2 = gp_Pnt(x2, y2, 500)
    myBody = BRepPrimAPI_MakeBox(p1, p2).Shape()
    rep = ifcopenshell.geom.tesselate(ifcfile.schema, myBody, 1.)
    ifcfile.add(rep)

    context = ifcfile.by_type("IfcGeometricRepresentationContext")[0]
    rep.Representations[0].ContextOfItems = context

    OwnerHistory = ifcfile.by_type('IfcOwnerHistory')[0]
    entity = ifcfile.create_entity('IfcAirTerminalBox', ifcopenshell.guid.new())
    entity.Name = name
    entity.OwnerHistory = OwnerHistory
    entity.Representation = rep
    entity.ObjectPlacement = create_ifclocalplacement(ifcfile)
    ifcfile.createIfcRelAggregates(ifcopenshell.guid.new(),
                                   OwnerHistory,
                                   None,
                                   None,
                                   ifcfile.by_type('IfcBuildingStorey')[0],
                                   [entity])
    colour([entity], color)


def new_pipe3d(input_skeleton):
    skeleton_length = input_skeleton.shape[0]
    array2 = TColgp_Array1OfPnt(1, skeleton_length)
    for idx in range(skeleton_length):
        ske = input_skeleton[idx]
        array2.SetValue(idx+1, gp_Pnt(int(ske[0]), int(ske[1]), int(ske[2])))
    ske = input_skeleton[0, :]
    ske1 = input_skeleton[1, :]
    ske2 = input_skeleton[-1, :]
    ske_dir = ske1 - ske
    bspline2 = GeomAPI_PointsToBSpline(array2).Curve()
    path_edge = BRepBuilderAPI_MakeEdge(bspline2).Edge()
    path_wire = BRepBuilderAPI_MakeWire(path_edge).Wire()

    dir = gp_Dir(ske_dir[0], ske_dir[1], ske_dir[2])
    circle = gp_Circ(gp_Ax2(gp_Pnt(int(ske[0]), int(ske[1]), int(ske[2])), dir), 100)
    profile_edge = BRepBuilderAPI_MakeEdge(circle).Edge()
    profile_wire = BRepBuilderAPI_MakeWire(profile_edge).Wire()
    profile_face = BRepBuilderAPI_MakeFace(profile_wire).Face()
    # pipe
    pipe = BRepOffsetAPI_MakePipe(path_wire, profile_face).Shape()
    sphere1 = BRepPrimAPI_MakeSphere(gp_Pnt(int(ske[0]), int(ske[1]), int(ske[2])), 100).Shape()
    sphere2 = BRepPrimAPI_MakeSphere(gp_Pnt(int(ske2[0]), int(ske2[1]), int(ske2[2])), 100).Shape()
    myBody1 = BRepAlgoAPI_Fuse(pipe, sphere1).Shape()
    myBody = BRepAlgoAPI_Fuse(myBody1, sphere2).Shape()
    return myBody


def create_ifc_entity(input_shpae, name, type, color):
    rep = ifcopenshell.geom.tesselate(ifcfile.schema, input_shpae, 1.)
    ifcfile.add(rep)

    context = ifcfile.by_type("IfcGeometricRepresentationContext")[0]
    rep.Representations[0].ContextOfItems = context

    OwnerHistory = ifcfile.by_type('IfcOwnerHistory')[0]
    pipe = ifcfile.create_entity(type, ifcopenshell.guid.new())
    pipe.Name = name
    pipe.OwnerHistory = OwnerHistory
    pipe.Representation = rep
    pipe.ObjectPlacement = create_ifclocalplacement(ifcfile)
    ifcfile.createIfcRelAggregates(ifcopenshell.guid.new(),
                                   OwnerHistory,
                                   None,
                                   None,
                                   ifcfile.by_type('IfcBuildingStorey')[0],
                                   [pipe])
    colour([pipe], color)


def generate_ifc_3d_pipe(input_data):
    sampled_skeleton = np.zeros((3, 2))
    sampled_skeleton[1, :] = input_data[input_data.shape[0] // 2, :]
    sampled_skeleton[0, :] = input_data[0, :]
    sampled_skeleton[-1, :] = input_data[int(input_data.shape[0]) - 1, :]
    z_coord = np.linspace(-70, 380, sampled_skeleton.shape[0])
    skeleton = np.insert(sampled_skeleton, 2, values=z_coord, axis=1)
    skeleton = skeleton[:, [1, 0, 2]]
    occ_pipe = new_pipe3d(skeleton)
    create_ifc_entity(occ_pipe, 'air pipe', 'IfcFlowSegment', [0.0, 1.0, 0.0])


def main():
    start = time.time()
    scale = 500/41

    with open("pkl/duct.pkl", "rb") as tf:
        duck = pickle.load(tf)
    contour = duck['contour']
    for key in contour.keys():
        generate_ifc_polygon('IfcFlowSegment', 'duct', contour[key], 0.0, 500.0, scale, [0.93, 0.75, 0.38])

    with open("pkl/equipment.pkl", "rb") as tf:
        equipment = pickle.load(tf)
    terminal_coords = equipment['air_terminal']
    terminal_width = terminal_height = 41
    for coord in terminal_coords:
        tx = coord[1] - int(0.5*terminal_width)
        ty = coord[0] + int(0.5*terminal_height)
        create_ifc_flow_terminal(500, 500, scale*tx, scale*ty, 0.0, [0.93, 0.75, 0.38])

    grille_coords = equipment['grill']
    for coord in grille_coords:
        coord = coord * scale
        x1, y1, x2, y2 = coord[0].item(), coord[1].item(), coord[2].item(), coord[3].item()
        create_ifc_rect(y1, x1, y2, x2, 'grille')

    fcu_10_coords = equipment['fcu_10']
    for coord in fcu_10_coords:
        coord = coord * scale
        x1, y1, x2, y2 = coord[0].item(), coord[1].item(), coord[2].item(), coord[3].item()
        create_ifc_rect(y1, x1, y2, x2, 'fcu_10', [1.0, 0.0, 0.0])

    fcu_4_coords = equipment['fcu_4']
    for coord in fcu_4_coords:
        coord = coord * scale
        x1, y1, x2, y2 = coord[0].item(), coord[1].item(), coord[2].item(), coord[3].item()
        create_ifc_rect(y1, x1, y2, x2, 'fcu_4', [1.0, 0.0, 0.0])

    fcu_8_coords = equipment['fcu_8']
    for coord in fcu_8_coords:
        coord = coord * scale
        x1, y1, x2, y2 = coord[0].item(), coord[1].item(), coord[2].item(), coord[3].item()
        create_ifc_rect(y1, x1, y2, x2, 'fcu_8', [1.0, 0.0, 0.0])

    ahu_coords = equipment['ahu']
    for coord in ahu_coords:
        coord = coord * scale
        x1, y1, x2, y2 = coord[0].item(), coord[1].item(), coord[2].item(), coord[3].item()
        create_ifc_rect(y1, x1, y2, x2, 'ahu', [0.0, 0.87, 0.0])

    contour = equipment['vav_4']
    for key in contour.keys():
        generate_ifc_polygon('IfcAirTerminalBox', 'vav_4', contour[key], 0.0, 500.0, scale, [0.0, 0.98, 0.99])
    contour = equipment['duct']
    for key in contour.keys():
        generate_ifc_polygon('IfcFlowSegment', 'duct', contour[key], 0.0, 500.0, scale, [0.93, 0.75, 0.38])

    air_pipe = equipment['air_pipe']
    for key in air_pipe.keys():
        generate_ifc_3d_pipe(air_pipe[key] * scale)
    ifcfile.write('demo1.ifc')
    end = time.time()
    print('used time: ', end-start)


if __name__ == '__main__':
    main()
