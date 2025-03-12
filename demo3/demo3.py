import math
import uuid
import time
import tempfile
import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_Cut
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeFace
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeWedge, BRepPrimAPI_MakeBox, BRepPrimAPI_MakeSphere, \
    BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeTorus
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.gp import gp_Pnt, gp_Ax2, gp_Dir, gp_Circ
from OCC.Display.SimpleGui import init_display
from ifcopenshell.util.element import copy_deep
from ifcopenshell.util.placement import get_local_placement
import numpy as np
import pickle
import random
from pipe_creation import create_sphere, Patcher


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
filename = "hello_duct.ifc"
timestamp = round(time.time())
timestring = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(timestamp))
creator = "Zaolin Pan"
organization = "UST"
application, application_version = "IfcOpenShell", "0.7"
project_globalid, project_name = create_guid(), "Hello Duct"

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
    else:
        print('keep original color')


def generate_ifc_polygon(type, name, _contours, _bottom_height, _wall_height, _scale, color):
    _contours = _scale * _contours[:, [1, 0]]
    # _contours = _scale * _contours
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


def generate_ifc_polygon_void(type, name, _contours, _contours1, _bottom_height, _wall_height, _scale, color):
    _contours = _scale * _contours[:, [1, 0]]
    _contours1 = _scale * _contours1[:, [1, 0]]
    storey = ifcfile.by_type('IfcBuildingStorey')[0]
    OwnerHistory = ifcfile.by_type('IfcOwnerHistory')[0]
    context = ifcfile.by_type("IfcGeometricRepresentationContext")[0]
    space_placement = create_ifclocalplacement(ifcfile)
    extrusion_placement = create_ifcaxis2placement(ifcfile, point=(0.0, 0.0, float(_bottom_height)))
    point_list_extrusion_area = _contours.tolist()
    point_list_extrusion_area.append(_contours[0, :].tolist())

    void_extrusion_area = _contours1.tolist()
    void_extrusion_area.append(_contours1[0, :].tolist())

    solid = create_ifcextrudedareasolid(ifcfile,
                                        point_list_extrusion_area,
                                        extrusion_placement,
                                        (0.0, 0.0, 1.0),
                                        float(_wall_height))

    void = create_ifcextrudedareasolid(ifcfile,
                                        void_extrusion_area,
                                        extrusion_placement,
                                        (0.0, 0.0, 1.0),
                                        float(_wall_height))

    opening_representation = ifcfile.createIfcShapeRepresentation(context, "Body", "SweptSolid", [void])
    opening_shape = ifcfile.createIfcProductDefinitionShape(None, None, [opening_representation])
    opening_element = ifcfile.createIfcOpeningElement(create_guid(), owner_history, "Opening", "An awesome opening",
                                                      None, space_placement, opening_shape, None)

    body_representation = ifcfile.createIfcShapeRepresentation(context, "Body", "SweptSolid", [solid])
    entity_shape = ifcfile.createIfcProductDefinitionShape(None, None, [body_representation])
    entity = ifcfile.create_entity(type, ifcopenshell.guid.new())
    entity.OwnerHistory = OwnerHistory
    entity.Name = name
    entity.Representation = entity_shape
    entity.ObjectPlacement = space_placement

    ifcfile.createIfcRelVoidsElement(create_guid(), owner_history, None, None, entity, opening_element)

    ifcfile.createIfcRelAggregates(ifcopenshell.guid.new(),
                                   OwnerHistory,
                                   None,
                                   None,
                                   storey,
                                   [entity])
    colour([entity], color)


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

    new_terminal = copy_deep(ifcfile, terminal)

    # 绕 x 轴旋转并平移到原点
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


def create_ifc_rect(x1, y1, x2, y2, name, h1, h2, color=None):
    p1 = gp_Pnt(x1, y1, float(h1))
    p2 = gp_Pnt(x2, y2, float(h2))
    myBody = BRepPrimAPI_MakeBox(p1, p2).Shape()
    rep = ifcopenshell.geom.tesselate(ifcfile.schema, myBody, 1.)
    ifcfile.add(rep)

    context = ifcfile.by_type("IfcGeometricRepresentationContext")[0]
    rep.Representations[0].ContextOfItems = context

    OwnerHistory = ifcfile.by_type('IfcOwnerHistory')[0]
    entity = ifcfile.create_entity('IfcBuildingElementProxy', ifcopenshell.guid.new())
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
    # 采样
    # idx = random.sample(range(1, (my_data.shape[0]//2)-1), 1)
    # idx.sort()
    # sampled_skeleton = np.zeros((len(idx)+2, 2))
    # for i, index in enumerate(idx):
    #     sampled_skeleton[i+1, :] = my_data[index, :]
    sampled_skeleton = np.zeros((3, 2))
    sampled_skeleton[1, :] = input_data[input_data.shape[0] // 2, :]
    sampled_skeleton[0, :] = input_data[0, :]
    sampled_skeleton[-1, :] = input_data[int(input_data.shape[0]) - 1, :]
    # 补充 z 轴信息
    z_coord = np.linspace(-70, 380, sampled_skeleton.shape[0])
    skeleton = np.insert(sampled_skeleton, 2, values=z_coord, axis=1)
    skeleton = skeleton[:, [1, 0, 2]]
    occ_pipe = new_pipe3d(skeleton)
    create_ifc_entity(occ_pipe, 'air pipe', 'IfcFlowSegment', [0.0, 1.0, 0.0])


def create_ifcextrudedareasolid_for_pipe(ifcfile, point, ifcaxis2placement, extrude_dir, extrusion, radius):
    location = create_ifcaxis2placement_2d(ifcfile, point=point, dir1=tuple((1.0, 0.0)))
    ifcclosedprofile = ifcfile.createIfcCircleProfileDef("AREA", "my pipe", location, radius)
    ifcdir = ifcfile.createIfcDirection(extrude_dir)
    ifcextrudedareasolid = ifcfile.createIfcExtrudedAreaSolid(ifcclosedprofile, ifcaxis2placement, ifcdir, extrusion)
    return ifcextrudedareasolid


def create_ifcpipe(_pipe_placement, circle_origin, extrusion_placement_origin, extrusion_placement_z_dir,
                   extrusion_placement_x_dir, extrusion_dir, extrusion_length, radius):
    extrusion_placement = create_ifcaxis2placement(ifcfile, extrusion_placement_origin, extrusion_placement_z_dir,
                                                   extrusion_placement_x_dir)
    solid = create_ifcextrudedareasolid_for_pipe(ifcfile, circle_origin, extrusion_placement, extrusion_dir,
                                                 extrusion_length, radius)
    body_representation = ifcfile.createIfcShapeRepresentation(context, "Body", "SweptSolid", [solid])
    product_shape = ifcfile.createIfcProductDefinitionShape(None, None, [body_representation])
    pipe = ifcfile.createIfcFlowSegment(create_guid(), owner_history, "Pipe", "An awesome pipe", None, _pipe_placement,
                                        product_shape, None)
    return pipe


def create_pipe_by_coords(data, _scale, color):
    my_radius = data['radius']
    pipes = []
    pipe_placement = create_ifclocalplacement(ifcfile, relative_to=storey_placement)
    if 'horizontal' in data.keys():
        h_lines = data['horizontal']
        h_lines = _scale * h_lines
        for idx in range(h_lines.shape[0]):
            h_line = h_lines[idx, :]
            extrusion_length = h_line[2] - h_line[0]
            # x 方向
            pipe = create_ifcpipe(_pipe_placement=pipe_placement,
                                  circle_origin=(0., 0.),
                                  extrusion_placement_origin=O,
                                  extrusion_placement_z_dir=(0.0, 1.0, 0.0),
                                  extrusion_placement_x_dir=(1.0, 0.0, 0.0),
                                  extrusion_dir=(0.0, 0.0, 1.0),
                                  extrusion_length=float(extrusion_length),
                                  radius=my_radius)
            new_pipe = copy_deep(ifcfile, pipe)
            # 平移
            movement = np.array([
                [1.0, 0.0, 0.0, float(h_line[1])],
                [0.0, 1.0, 0.0, float(h_line[0])],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]])
            patcher = Patcher(ifcfile)
            patcher.patch(new_pipe, movement)
            new_pipe.GlobalId = ifcopenshell.guid.new()
            pipes.append(new_pipe)
    if 'vertical' in data.keys():
        v_lines = data['vertical']
        v_lines = _scale * v_lines
        for idx in range(v_lines.shape[0]):
            v_line = v_lines[idx, :]
            extrusion_length = v_line[3] - v_line[1]
            # y 方向
            pipe = create_ifcpipe(_pipe_placement=pipe_placement,
                                  circle_origin=(0., 0.),
                                  extrusion_placement_origin=O,
                                  extrusion_placement_z_dir=(1.0, 0.0, 0.0),
                                  extrusion_placement_x_dir=(0.0, 0.0, 1.0),
                                  extrusion_dir=(0.0, 0.0, 1.0),
                                  extrusion_length=float(extrusion_length),
                                  radius=my_radius)
            new_pipe = copy_deep(ifcfile, pipe)
            # 平移
            movement = np.array([
                [1.0, 0.0, 0.0, float(v_line[1])],
                [0.0, 1.0, 0.0, float(v_line[0])],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]])
            patcher = Patcher(ifcfile)
            patcher.patch(new_pipe, movement)
            new_pipe.GlobalId = ifcopenshell.guid.new()
            pipes.append(new_pipe)
    if 'inclined' in data.keys():
        i_lines = data['inclined']
        i_lines = _scale * i_lines
        k = (i_lines[:, 3] - i_lines[:, 1]) / (i_lines[:, 2] - i_lines[:, 0])

        for idx in range(i_lines.shape[0]):
            i_line = i_lines[idx, :]
            extrusion_length = np.sqrt(np.sum(np.square(i_line[2:] - i_line[:2])))
            angle = math.degrees(math.atan(k[idx]))
            # x 方向
            pipe = create_ifcpipe(_pipe_placement=pipe_placement,
                                  circle_origin=(0., 0.),
                                  extrusion_placement_origin=O,
                                  extrusion_placement_z_dir=(1.0, 0.0, 0.0),
                                  extrusion_placement_x_dir=(0.0, 0.0, 1.0),
                                  extrusion_dir=(0.0, 0.0, 1.0),
                                  extrusion_length=float(extrusion_length),
                                  radius=my_radius)
            new_pipe = copy_deep(ifcfile, pipe)
            # 平移和绕 z 轴旋转
            movement = np.array([
                [math.sin(math.radians(angle)), math.cos(math.radians(angle)), 0.0, float(i_line[1])],
                [math.cos(math.radians(angle)), -math.sin(math.radians(angle)), 0.0, float(i_line[0])],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0]])
            patcher = Patcher(ifcfile)
            patcher.patch(new_pipe, movement)
            new_pipe.GlobalId = ifcopenshell.guid.new()
            pipes.append(new_pipe)
    ifcfile.createIfcRelContainedInSpatialStructure(create_guid(), owner_history,
                                                    "Building Storey Container", None, pipes, building_storey)
    create_pipe_material(pipes)
    colour(pipes, color)


def create_pipe_by_occ(radius, x1_, y1_, z1_, x2_, y2_, z2_):
    loc_ = gp_Ax2(gp_Pnt(x1_, y1_, z1_), gp_Dir(x2_-x1_, y2_-y1_, z2_-z1_))
    dist = math.sqrt((x2_-x1_)**2 + (y2_-y1_)**2 + (z2_-z1_)**2)
    pipe_ = BRepPrimAPI_MakeCylinder(loc_, radius, dist).Shape()
    return pipe_


def create_l_shape_fitting(cx, cy, cz, r1, r2, quadrant):
    if quadrant == 1:
        # 第一象限
        loc = gp_Ax2(gp_Pnt(cx, cy+r1, cz), gp_Dir(0, -1, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx+r1, cy, cz), gp_Dir(-1, 0, 0))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()

        S1 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx+r1, cy+r1, cz), gp_Dir(0, 0, 1)), r1, r2, math.pi).Shape()
        S2 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx+r1, cy+r1, cz), gp_Dir(0, 0, 1)), r1, r2, 3*math.pi/2).Shape()
        S = BRepAlgoAPI_Cut(S2, S1).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 2:
        # 第二象限
        loc = gp_Ax2(gp_Pnt(cx, cy+r1, cz), gp_Dir(0, -1, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx-r1, cy, cz), gp_Dir(-1, 0, 0))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()

        S1 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx-r1, cy+r1, cz), gp_Dir(0, 0, 1)), r1, r2, 3*math.pi/2).Shape()
        S2 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx-r1, cy+r1, cz), gp_Dir(0, 0, 1)), r1, r2, 2*math.pi).Shape()
        S = BRepAlgoAPI_Cut(S2, S1).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 3:
        # 第三象限
        loc = gp_Ax2(gp_Pnt(cx, cy - r1, cz), gp_Dir(0, -1, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx - r1, cy, cz), gp_Dir(-1, 0, 0))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()

        S = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx - r1, cy - r1, cz), gp_Dir(0, 0, 1)), r1, r2,
                                   math.pi / 2).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 4:
        # 第四象限
        loc = gp_Ax2(gp_Pnt(cx, cy - r1, cz), gp_Dir(0, -1, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx + r1, cy, cz), gp_Dir(-1, 0, 0))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()

        S = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx + r1, cy - r1, cz), gp_Dir(0, 0, -1)), r1, r2,
                                   math.pi / 2).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 5:
        # 第五卦限
        loc = gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(-1, 0, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx - r1, cy, cz + r1), gp_Dir(0, 0, -1))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()
        S3 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz + r1), gp_Dir(0, 1, 0)), r1, r2, math.pi).Shape()
        S4 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz + r1), gp_Dir(0, 1, 0)), r1, r2, 3 * math.pi / 2).Shape()
        S = BRepAlgoAPI_Cut(S4, S3).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 6:
        loc = gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(0, -1, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx, cy - r1, cz + r1), gp_Dir(0, 0, -1))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()
        S3 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz + r1), gp_Dir(1, 0, 0)), r1, r2, math.pi / 2).Shape()
        S4 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz + r1), gp_Dir(1, 0, 0)), r1, r2, math.pi).Shape()
        S = BRepAlgoAPI_Cut(S4, S3).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 7:
        loc = gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(1, 0, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx + r1, cy, cz + r1), gp_Dir(0, 0, -1))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()
        S3 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz + r1), gp_Dir(0, 1, 0)), r1, r2, math.pi / 2).Shape()
        S4 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz + r1), gp_Dir(0, 1, 0)), r1, r2, math.pi).Shape()
        S = BRepAlgoAPI_Cut(S4, S3).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 8:
        loc = gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(0, 1, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx, cy + r1, cz + r1), gp_Dir(0, 0, -1))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()
        S3 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz + r1), gp_Dir(1, 0, 0)), r1, r2, math.pi).Shape()
        S4 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz + r1), gp_Dir(1, 0, 0)), r1, r2, 3 * math.pi / 2).Shape()
        S = BRepAlgoAPI_Cut(S4, S3).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 9:
        loc = gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(-1, 0, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx - r1, cy, cz - r1), gp_Dir(0, 0, 1))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()
        S3 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz - r1), gp_Dir(0, 1, 0)), r1, r2, 3 * math.pi / 2).Shape()
        S4 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz - r1), gp_Dir(0, 1, 0)), r1, r2, 2 * math.pi).Shape()
        S = BRepAlgoAPI_Cut(S4, S3).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 10:
        loc = gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(0, -1, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx, cy - r1, cz - r1), gp_Dir(0, 0, 1))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()
        S = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz - r1), gp_Dir(1, 0, 0)), r1, r2, math.pi / 2).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 11:
        loc = gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(1, 0, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx + r1, cy, cz - r1), gp_Dir(0, 0, 1))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()
        S = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz - r1), gp_Dir(0, 1, 0)), r1, r2, math.pi / 2).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    if quadrant == 12:
        loc = gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(0, 1, 0))
        cylinder1 = BRepPrimAPI_MakeCylinder(loc, 1.3 * r2, 10).Shape()
        loc1 = gp_Ax2(gp_Pnt(cx, cy + r1, cz - r1), gp_Dir(0, 0, 1))
        cylinder2 = BRepPrimAPI_MakeCylinder(loc1, 1.3 * r2, 10).Shape()
        S3 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz - r1), gp_Dir(1, 0, 0)), r1, r2, 3 * math.pi / 2).Shape()
        S4 = BRepPrimAPI_MakeTorus(gp_Ax2(gp_Pnt(cx, cy, cz - r1), gp_Dir(1, 0, 0)), r1, r2, 2 * math.pi).Shape()
        S = BRepAlgoAPI_Cut(S4, S3).Shape()
        my_body = BRepAlgoAPI_Fuse(cylinder1, cylinder2).Shape()
        l_shape_fitting = BRepAlgoAPI_Fuse(my_body, S).Shape()

    return l_shape_fitting


def create_v_shape_fitting(cx, cy, cz, radius, dir_x1, dir_y1, dir_z1, dir_x2, dir_y2, dir_z2):
    cylinder_a1 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(dir_x1, dir_y1, dir_z1)), 1.3 * radius,
                                           1.1 * radius).Shape()
    cylinder_a2 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(dir_x1, dir_y1, dir_z1)), 1.3 * radius,
                                           radius).Shape()
    cylinder_a = BRepAlgoAPI_Cut(cylinder_a1, cylinder_a2).Shape()
    cylinder_b1 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(dir_x2, dir_y2, dir_z2)), 1.3 * radius,
                                           1.1 * radius).Shape()
    cylinder_b2 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(dir_x2, dir_y2, dir_z2)), 1.3 * radius,
                                           radius).Shape()
    cylinder_b = BRepAlgoAPI_Cut(cylinder_b1, cylinder_b2).Shape()
    pipe1 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(dir_x1, dir_y1, dir_z1)), radius, radius).Shape()
    pipe2 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(dir_x2, dir_y2, dir_z2)), radius, radius).Shape()
    sphere = create_sphere(cx, cy, cz, radius)
    my_body = BRepAlgoAPI_Fuse(sphere, pipe1).Shape()
    my_body = BRepAlgoAPI_Fuse(my_body, pipe2).Shape()
    my_body = BRepAlgoAPI_Fuse(my_body, cylinder_a).Shape()
    v_shape_fitting = BRepAlgoAPI_Fuse(my_body, cylinder_b).Shape()
    # create_ifc_entity(v_shape_fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])
    return v_shape_fitting


def create_pipe_and_fitting(my_dict):
    radius = my_dict['radius']
    v_lines = my_dict['vertical']
    h_lines = my_dict['horizontal']
    i_lines = my_dict['inclined']
    lines = np.vstack((v_lines, h_lines, i_lines)) * scale
    new_lines = lines[np.argsort(lines[:, 0])]

    for idx, current_line in enumerate(new_lines):
        tmp_line = current_line.copy()
        if not ((tmp_line[0] <= tmp_line[2]) and (tmp_line[1] <= tmp_line[3])):
            new_lines[idx, :2] = tmp_line[2:]
            new_lines[idx, 2:] = tmp_line[:2]
    next_pipe = []
    for idx, line in enumerate(new_lines):
        if idx == 0:
            first_fitting = create_l_shape_fitting(line[1], line[0], 0, 1.5 * radius, radius, 5)
            create_ifc_entity(first_fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

        if idx == new_lines.shape[0] - 1:
            # pipe = create_pipe_by_occ(radius, line[1], line[0], 0, line[3], line[2], 0)
            pipe = next_pipe.pop()
        else:
            fitting_loc = line[2:]
            v1 = line[:2] - fitting_loc
            v2 = new_lines[idx + 1, 2:] - fitting_loc
            flag = np.dot(v1, v2)

            if flag == 0:
                # 说明这两根直线是垂直的，用 L 型弯头
                if v1[0] == 0:
                    v_flag = v1[1]
                    h_flag = v2[0]
                else:
                    h_flag = v1[0]
                    v_flag = v2[1]

                if (h_flag <= 0) and (v_flag >= 0):
                    # print('第四象限')
                    quadrant = 4
                elif (h_flag >= 0) and (v_flag >= 0):
                    # print('第三象限')
                    quadrant = 3
                elif (h_flag >= 0) and (v_flag <= 0):
                    # print('第二象限')
                    quadrant = 2
                elif (h_flag <= 0) and (v_flag <= 0):
                    # print('第一象限')
                    quadrant = 1
                fitting = create_l_shape_fitting(fitting_loc[1], fitting_loc[0], 0, 1.5 * radius, radius, quadrant)
                create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])
                if line[0] == line[2]:
                    pipe = create_pipe_by_occ(radius, line[1], line[0], 0, line[3] - 1.5 * radius, line[2], 0)
                    if v_flag <= 0:
                        new_lines[idx + 1, 0] = new_lines[idx + 1, 0] + 1.5 * radius
                    else:
                        new_lines[idx + 1, 0] = new_lines[idx + 1, 0] - 1.5 * radius
                if line[1] == line[3]:
                    pipe = create_pipe_by_occ(radius, line[1], line[0], 0, line[3], line[2] - 1.5 * radius, 0)
                    if h_flag <= 0:
                        new_lines[idx + 1, 1] = new_lines[idx + 1, 1] + 1.5 * radius
                    else:
                        new_lines[idx + 1, 1] = new_lines[idx + 1, 1] - 1.5 * radius
            else:
                cx, cy, cz = fitting_loc[1], fitting_loc[0], 0.0

                if next_pipe:
                    pipe = next_pipe.pop()
                else:
                    pipe = create_pipe_by_occ(radius, line[1], line[0], 0, line[3], line[2], 0)

                # cylinder_a1 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(v1[1], v1[0], 0)), 1.3*radius,
                #                                  1.1*radius).Shape()
                # cylinder_a2 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(v1[1], v1[0], 0)), 1.3*radius,
                #                                  radius).Shape()
                # cylinder_a = BRepAlgoAPI_Cut(cylinder_a1, cylinder_a2).Shape()
                # cylinder_b1 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(v2[1], v2[0], 0)), 1.3*radius,
                #                                  1.1*radius).Shape()
                # cylinder_b2 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(v2[1], v2[0], 0)), 1.3*radius,
                #                                  radius).Shape()
                # cylinder_b = BRepAlgoAPI_Cut(cylinder_b1, cylinder_b2).Shape()
                # pipe1 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(v1[1], v1[0], 0)), radius, radius).Shape()
                # pipe2 = BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(cx, cy, cz), gp_Dir(v2[1], v2[0], 0)), radius, radius).Shape()
                # my_body = BRepAlgoAPI_Fuse(sphere, pipe1).Shape()
                # my_body = BRepAlgoAPI_Fuse(my_body, pipe2).Shape()
                # my_body = BRepAlgoAPI_Fuse(my_body, cylinder_a).Shape()
                # fitting = BRepAlgoAPI_Fuse(my_body, cylinder_b).Shape()
                # create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])
                fitting = create_v_shape_fitting(cx, cy, cz, radius, v1[1], v1[0], 0, v2[1], v2[0], 0)
                pipe = BRepAlgoAPI_Cut(pipe, fitting).Shape()
                next_line = new_lines[idx + 1, :]
                next_pipe_ = create_pipe_by_occ(radius, next_line[1], next_line[0], 0, next_line[3], next_line[2], 0)
                next_pipe_ = BRepAlgoAPI_Cut(next_pipe_, fitting).Shape()
                next_pipe.append(next_pipe_)
        create_ifc_entity(pipe, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])


def create_t_shape_fitting(radius1, radius2, s, c, e, t, height=0.0):
    if len(s) < 3:
        s = np.hstack((s[1], s[0], height))

    if len(c) < 3:
        c = np.hstack((c[1], c[0], height))

    if len(e) < 3:
        e = np.hstack((e[1], e[0], height))

    cylinder_length = 1.5 * radius1
    x_dir = s - c
    cylinder_x = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(c[0], c[1], c[2]), gp_Dir(x_dir[0], x_dir[1], x_dir[2])),
        radius1, cylinder_length).Shape()
    cylinder_x1 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(c[0], c[1], c[2]), gp_Dir(x_dir[0], x_dir[1], x_dir[2])),
        1.3*radius1, cylinder_length + 0.1*radius1).Shape()
    cylinder_x2 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(c[0], c[1], c[2]), gp_Dir(x_dir[0], x_dir[1], x_dir[2])),
        1.3*radius1, cylinder_length).Shape()
    cylinder_x3 = BRepAlgoAPI_Cut(cylinder_x1, cylinder_x2).Shape()

    y_dir = e - c
    cylinder_y = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(c[0], c[1], c[2]), gp_Dir(y_dir[0], y_dir[1], y_dir[2])),
        radius1, cylinder_length).Shape()
    cylinder_y1 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(c[0], c[1], c[2]), gp_Dir(y_dir[0], y_dir[1], y_dir[2])),
        1.3*radius1, cylinder_length + 0.1*radius1).Shape()
    cylinder_y2 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(c[0], c[1], c[2]), gp_Dir(y_dir[0], y_dir[1], y_dir[2])),
        1.3*radius1, cylinder_length).Shape()
    cylinder_y3 = BRepAlgoAPI_Cut(cylinder_y1, cylinder_y2).Shape()

    z_dir = t
    cylinder_z = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(c[0], c[1], c[2]), gp_Dir(z_dir[0], z_dir[1], z_dir[2])),
        radius2, 1.5 * cylinder_length).Shape()
    cylinder_z1 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(c[0], c[1], c[2]), gp_Dir(z_dir[0], z_dir[1], z_dir[2])),
        1.3*radius2, 1.5 * cylinder_length + 0.1*radius2).Shape()
    cylinder_z2 = BRepPrimAPI_MakeCylinder(
        gp_Ax2(gp_Pnt(c[0], c[1], c[2]), gp_Dir(z_dir[0], z_dir[1], z_dir[2])),
        1.3*radius2, 1.5 * cylinder_length).Shape()
    cylinder_z3 = BRepAlgoAPI_Cut(cylinder_z1, cylinder_z2).Shape()

    my_body = BRepAlgoAPI_Fuse(cylinder_x, cylinder_x3).Shape()
    my_body = BRepAlgoAPI_Fuse(my_body, cylinder_y).Shape()
    my_body = BRepAlgoAPI_Fuse(my_body, cylinder_y3).Shape()
    my_body = BRepAlgoAPI_Fuse(my_body, cylinder_z).Shape()
    fitting = BRepAlgoAPI_Fuse(my_body, cylinder_z3).Shape()
    return fitting


def generate_tube(input_data, r1_, r2_):
    lines = input_data['lines']
    pipe_list, fitting_list = [], []
    for line in lines:
        line = line * scale
        pipe_ = create_pipe_by_occ(r1_, line[1], line[0], 0, line[3], line[2], 0)
        pipe_list.append(pipe_)
        # create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])

    l_shape_fittings = input_data['l_shape_fittings']
    for fitting_info in l_shape_fittings:
        c = fitting_info[0] * scale
        quadrant = fitting_info[1]
        fitting = create_l_shape_fitting(c[1], c[0], 0, 1.5 * r1_, r1_, quadrant)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    v_shape_fittings = input_data['v_shape_fittings']
    for fitting_info in v_shape_fittings:
        c = fitting_info[0] * scale
        v1 = fitting_info[1].astype('float')
        v2 = fitting_info[2].astype('float')
        fitting = create_v_shape_fitting(c[1], c[0], 0, r1_, v1[1], v1[0], 0, v2[1], v2[0], 0)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    t_shape_fittings = input_data['t_shape_fittings']
    for fitting_info in t_shape_fittings:
        s = fitting_info[0] * scale
        c = fitting_info[1] * scale
        e = fitting_info[2] * scale
        fitting = create_t_shape_fitting(r1_, r2_, s, c, e, np.array([0., 0., 1.]))
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    my_body = None
    while fitting_list:
        if my_body is None:
            my_body = fitting_list.pop()
        my_body = BRepAlgoAPI_Fuse(my_body, fitting_list.pop()).Shape()
    # create_ifc_entity(my_body, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    for pipe in pipe_list:
        trimmed_pipe = BRepAlgoAPI_Cut(pipe, my_body).Shape()
        create_ifc_entity(trimmed_pipe, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])


def generate_up_tube(input_data):
    input_lines = input_data['lines']
    for line in input_lines:
        line = line * scale
        pipe_ = create_pipe_by_occ(r2, line[1], line[0], 1.5*1.5*r1+0.1*r2 + 1.5*r2,
                                   line[3], line[2], 1.5*1.5*r1+0.1*r2 + 1.5*r2)
        create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])

    l_shape_fittings = input_data['l_shape_fittings']
    for fitting_info in l_shape_fittings:
        c = fitting_info[0] * scale
        quadrant = fitting_info[1]
        fitting = create_l_shape_fitting(c[1], c[0], 1.5*1.5*r1+0.1*r2 + 1.5*r2, 1.5 * r2, r2, quadrant)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])


if __name__ == '__main__':
    start = time.time()
    scale = 300/58  # 这个系数要重新算一下

    # %% wall
    with open("pkl/wall.pkl", "rb") as tf:
        wall = pickle.load(tf)
    contour = wall['contour']
    for key in contour.keys():
        generate_ifc_polygon('IfcWall', 'wall', contour[key], 0.0, 4500.0, scale, [0.73, 0.73, 0.73])
    void = wall['void']
    generate_ifc_polygon_void('IfcWall', 'wall', void[0], void[1], 0.0, 4500.0, scale, [0.73, 0.73, 0.73])
    # %% duct
    # with open("pkl/duct.pkl", "rb") as tf:
    #     duct = pickle.load(tf)
    # contour = duct['contour']
    # for key in contour.keys():
    #     generate_ifc_polygon('IfcFlowSegment', 'duct', contour[key], 0.0, 50.0, scale, [0.93, 0.75, 0.38])
    # %% slab
    with open("pkl/slab.pkl", "rb") as tf:
        slab = pickle.load(tf)
    contour = slab['contour']
    for key in contour.keys():
        generate_ifc_polygon('IfcSlab', 'slab', contour[key], 0.0, 200.0, scale, [0.73, 0.73, 0.73])
    # %% machine
    with open("pkl/machine.pkl", "rb") as tf:
        machine = pickle.load(tf)
    pump_coords = machine['pump']
    for coord in pump_coords:
        coord = coord * scale
        x1, y1, x2, y2 = coord[0].item(), coord[1].item(), coord[2].item(), coord[3].item()
        create_ifc_rect(y1, x1, y2, x2, 'pump', 200.0, 800.0, [1.0, 0., 0.])
    chiller_coords = machine['chiller']
    for coord in chiller_coords:
        coord = coord * scale
        x1, y1, x2, y2 = coord[0].item(), coord[1].item(), coord[2].item(), coord[3].item()
        create_ifc_rect(y1, x1, y2, x2, 'chiller', 200.0, 1500.0, [1.0, 0.0, 0.])
    # %%
    # with open("pkl/main-pipe.pkl", "rb") as tf:
    #     my_data = pickle.load(tf)
    #     for key in my_data.keys():
    #         # create_pipe_by_coords(my_data[key], scale, [0.0, 1.0, 0.])
    #         create_pipe_and_fitting(my_data[key])
    # %%
    with open("pkl/tube1.pkl", "rb") as tf:
        my_tube1 = pickle.load(tf)
    # with open("pkl/tube2.pkl", "rb") as tf:
    #     my_tube2 = pickle.load(tf)
    with open("pkl/new_tube2.pkl", "rb") as tf:
        my_tube2 = pickle.load(tf)
    r1 = my_tube1['radius1']
    r2 = my_tube1['radius2']
    generate_tube(my_tube1, r1, r2)
    generate_tube(my_tube2, r1, r2)
    # %% upper tube 1
    with open("pkl/purple1.pkl", "rb") as tf:
        purple1 = pickle.load(tf)
    with open("pkl/purple2.pkl", "rb") as tf:
        purple2 = pickle.load(tf)
    generate_up_tube(purple1)
    generate_up_tube(purple2)
    # %% upper tube 2
    with open('pkl/purple3.pkl', 'rb') as tf:
        purple3 = pickle.load(tf)
    lines = purple3['lines']
    pipe_list, fitting_list = [], []
    for line in lines:
        line = line * scale
        pipe_ = create_pipe_by_occ(r2, line[1], line[0], 1.5*1.5*r1+0.1*r2 + 1.5*r2, line[3], line[2], 1.5*1.5*r1+0.1*r2 + 1.5*r2)
        pipe_list.append(pipe_)
        # create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
    l_shape_fittings = purple3['l_shape_fittings']
    for fitting_info in l_shape_fittings:
        c = fitting_info[0] * scale
        quadrant = fitting_info[1]
        fitting = create_l_shape_fitting(c[1], c[0], 1.5*1.5*r1+0.1*r2 + 1.5*r2, 1.5 * r2, r2, quadrant)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])
    t_shape_fittings = purple3['t_shape_fittings']
    for fitting_info in t_shape_fittings:
        s = fitting_info[0] * scale
        c = fitting_info[1] * scale
        e = fitting_info[2] * scale
        t = np.array([fitting_info[3][1], fitting_info[3][0], 0.0])
        fitting = create_t_shape_fitting(r2, r2, s, c, e, t, height=1.5*1.5*r1+0.1*r2 + 1.5*r2)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])
    my_body = None
    while fitting_list:
        if my_body is None:
            my_body = fitting_list.pop()
        my_body = BRepAlgoAPI_Fuse(my_body, fitting_list.pop()).Shape()

    for pipe in pipe_list:
        trimmed_pipe = BRepAlgoAPI_Cut(pipe, my_body).Shape()
        create_ifc_entity(trimmed_pipe, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
    # %% lower tube 1
    with open('pkl/green2.pkl', 'rb') as tf:
        green2 = pickle.load(tf)
    green2_height = -250.0
    lines = green2['lines']
    pipe_list, fitting_list = [], []
    for line in lines:
        line = line * scale
        pipe_ = create_pipe_by_occ(r2, line[1], line[0], green2_height, line[3], line[2], green2_height)
        pipe_list.append(pipe_)
        # create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])

    l_shape_fittings = green2['l_shape_fittings']
    for fitting_info in l_shape_fittings:
        c = fitting_info[0] * scale
        quadrant = fitting_info[1]
        fitting = create_l_shape_fitting(c[1], c[0], green2_height, 1.5 * r2, r2, quadrant)
        if quadrant >= 5:
            pipe_ = create_pipe_by_occ(r2, c[1], c[0]+1.5*r2, green2_height+1.5*r2, c[1], c[0]+1.5*r2, 0)
            create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    t_shape_fittings = green2['t_shape_fittings']
    for fitting_info in t_shape_fittings:
        s = fitting_info[0] * scale
        c = fitting_info[1] * scale
        e = fitting_info[2] * scale
        t = np.array([fitting_info[3][1], fitting_info[3][0], 0.0])
        fitting = create_t_shape_fitting(r2, r2, s, c, e, t, height=green2_height)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    my_body = None
    while fitting_list:
        if my_body is None:
            my_body = fitting_list.pop()
        my_body = BRepAlgoAPI_Fuse(my_body, fitting_list.pop()).Shape()

    for pipe in pipe_list:
        trimmed_pipe = BRepAlgoAPI_Cut(pipe, my_body).Shape()
        create_ifc_entity(trimmed_pipe, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
    # %% green
    with open('pkl/green.pkl', 'rb') as tf:
        green = pickle.load(tf)
    green_height = 1.5*r2
    lines = green['lines']
    pipe_list, fitting_list = [], []
    for line in lines:
        line = line * scale
        pipe_ = create_pipe_by_occ(r2, line[1], line[0], green_height, line[3], line[2], green_height)
        pipe_list.append(pipe_)
        # create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])

    l_shape_fittings = green['l_shape_fittings']
    for fitting_info in l_shape_fittings:
        c = fitting_info[0] * scale
        quadrant = fitting_info[1]
        fitting = create_l_shape_fitting(c[1], c[0], green_height, 1.5 * r2, r2, quadrant)
        # if quadrant >= 5:
        #     pipe_ = create_pipe_by_occ(r2, c[1], c[0]+1.5*r2, green2_height+1.5*r2, c[1], c[0]+1.5*r2, 0)
        #     create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    t_shape_fittings = green['t_shape_fittings']
    for fitting_info in t_shape_fittings:
        s = fitting_info[0] * scale
        c = fitting_info[1] * scale
        e = fitting_info[2] * scale
        t = np.array([fitting_info[3][1], fitting_info[3][0], 0.0])
        fitting = create_t_shape_fitting(r2, r2, s, c, e, t, height=green_height)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    v_shape_fittings = green['v_shape_fittings']
    for fitting_info in v_shape_fittings:
        c = fitting_info[0] * scale
        v1 = fitting_info[1].astype('float')
        v2 = fitting_info[2].astype('float')
        fitting = create_v_shape_fitting(c[1], c[0], green_height, r2,
                                         v1[1], v1[0], 0.0,
                                         v2[1], v2[0], 0.0)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    my_body = None
    while fitting_list:
        if my_body is None:
            my_body = fitting_list.pop()
        my_body = BRepAlgoAPI_Fuse(my_body, fitting_list.pop()).Shape()

    for pipe in pipe_list:
        trimmed_pipe = BRepAlgoAPI_Cut(pipe, my_body).Shape()
        create_ifc_entity(trimmed_pipe, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
    # %% deep blue
    with open('pkl/deep_blue.pkl', 'rb') as tf:
        deep_blue = pickle.load(tf)
    deep_blue_height = -200.0
    lines = deep_blue['lines']
    pipe_list, fitting_list = [], []
    for line in lines:
        line = line * scale
        pipe_ = create_pipe_by_occ(r2, line[1], line[0], deep_blue_height,
                                   line[3], line[2], deep_blue_height)
        pipe_list.append(pipe_)
        # create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])

    l_shape_fittings = deep_blue['l_shape_fittings']
    for fitting_info in l_shape_fittings:
        c = fitting_info[0] * scale
        quadrant = fitting_info[1]
        if quadrant == 8:
            pipe_ = create_pipe_by_occ(r2, c[1], c[0]+1.5*r2, deep_blue_height+1.5*r2, c[1], c[0]+1.5*r2, 0)
            create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
        if quadrant == 7:
            pipe_ = create_pipe_by_occ(r2, c[1]+1.5*r2, c[0], deep_blue_height+1.5*r2, c[1]+1.5*r2, c[0], 0)
            create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
        fitting = create_l_shape_fitting(c[1], c[0], deep_blue_height, 1.5 * r2, r2, quadrant)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    t_shape_fittings = deep_blue['t_shape_fittings']
    for fitting_info in t_shape_fittings:
        s = fitting_info[0] * scale
        c = fitting_info[1] * scale
        e = fitting_info[2] * scale
        t = np.array([fitting_info[3][1], fitting_info[3][0], 0.0])
        fitting = create_t_shape_fitting(r2, r2, s, c, e, t, height=deep_blue_height)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    my_body = None
    while fitting_list:
        if my_body is None:
            my_body = fitting_list.pop()
        my_body = BRepAlgoAPI_Fuse(my_body, fitting_list.pop()).Shape()

    for pipe in pipe_list:
        trimmed_pipe = BRepAlgoAPI_Cut(pipe, my_body).Shape()
        create_ifc_entity(trimmed_pipe, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
    # %% green 3
    with open('pkl/green3.pkl', 'rb') as tf:
        green3 = pickle.load(tf)
    green3_height = -200.0
    lines = green3['lines']
    pipe_list, fitting_list = [], []
    for line in lines:
        line = line * scale
        pipe_ = create_pipe_by_occ(r2, line[1], line[0], green3_height,
                                   line[3], line[2], green3_height)
        pipe_list.append(pipe_)
        # create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])

    l_shape_fittings = green3['l_shape_fittings']
    for fitting_info in l_shape_fittings:
        c = fitting_info[0] * scale
        quadrant = fitting_info[1]
        if quadrant == 8:
            pipe_ = create_pipe_by_occ(r2, c[1], c[0]+1.5*r2, green3_height+1.5*r2, c[1], c[0]+1.5*r2, 0)
            create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
        if quadrant == 5:
            pipe_ = create_pipe_by_occ(r2, c[1] - 1.5*r2, c[0], green3_height+1.5*r2, c[1] - 1.5*r2, c[0], 0)
            create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
        fitting = create_l_shape_fitting(c[1], c[0], green3_height, 1.5 * r2, r2, quadrant)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    t_shape_fittings = green3['t_shape_fittings']
    for fitting_info in t_shape_fittings:
        s = fitting_info[0] * scale
        c = fitting_info[1] * scale
        e = fitting_info[2] * scale
        t = np.array([fitting_info[3][1], fitting_info[3][0], 0.0])
        fitting = create_t_shape_fitting(r2, r2, s, c, e, t, height=green3_height)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    my_body = None
    while fitting_list:
        if my_body is None:
            my_body = fitting_list.pop()
        my_body = BRepAlgoAPI_Fuse(my_body, fitting_list.pop()).Shape()

    for pipe in pipe_list:
        trimmed_pipe = BRepAlgoAPI_Cut(pipe, my_body).Shape()
        create_ifc_entity(trimmed_pipe, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])
    # %% green 1
    with open('pkl/green1.pkl', 'rb') as tf:
        green1 = pickle.load(tf)
    green1_height = -500.0
    lines = green1['lines']
    pipe_list, fitting_list = [], []
    for line in lines:
        line = line * scale
        pipe_ = create_pipe_by_occ(r2, line[1], line[0], green1_height,
                                   line[3], line[2], green1_height)
        pipe_list.append(pipe_)
        # create_ifc_entity(pipe_, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])

    l_shape_fittings = green1['l_shape_fittings']
    for fitting_info in l_shape_fittings:
        c = fitting_info[0] * scale
        quadrant = fitting_info[1]
        fitting = create_l_shape_fitting(c[1], c[0], green1_height, 1.5 * r2, r2, quadrant)
        fitting_list.append(fitting)
        create_ifc_entity(fitting, 'fitting', 'IfcFlowFitting', [0.0, 0.0, 1.])

    my_body = None
    while fitting_list:
        if my_body is None:
            my_body = fitting_list.pop()
        my_body = BRepAlgoAPI_Fuse(my_body, fitting_list.pop()).Shape()

    for pipe in pipe_list:
        trimmed_pipe = BRepAlgoAPI_Cut(pipe, my_body).Shape()
        create_ifc_entity(trimmed_pipe, 'pipe', 'IfcFlowSegment', [0.0, 1.0, 0.])

    # %%
    ifcfile.write('demo3.ifc')
    end = time.time()
    print('used time: ', end-start)


