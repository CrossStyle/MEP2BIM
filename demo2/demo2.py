import uuid
import time
import tempfile
import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util
import pickle


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


create_guid = lambda: ifcopenshell.guid.compress(uuid.uuid1().hex)
# IFC template creation
filename = "demo2.ifc"
timestamp = round(time.time())
timestring = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(timestamp))
creator = "Zaolin Pan"
organization = "HKUST"
application, application_version = "IfcOpenShell", "0.7"
project_globalid, project_name = create_guid(), "the MEP system 2"

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


def main():
    start = time.time()
    scale = 10600/467

    with open("pkl/duct.pkl", "rb") as tf:
        duct = pickle.load(tf)
    duct_contour = duct['contour']
    for key in duct_contour.keys():
        generate_ifc_polygon('IfcFlowSegment', 'duct', duct_contour[key], 0.0, 500.0, scale, [0.93, 0.75, 0.38])

    with open("pkl/muffler.pkl", "rb") as tf:
        muffler = pickle.load(tf)
    muffler_contour = muffler['contour']
    for key in muffler_contour.keys():
        generate_ifc_polygon('IfcAirTerminalBox', 'muffler', muffler_contour[key], 0.0, 500.0, scale, [0.0, 1., 0.0])

    with open("pkl/grille.pkl", "rb") as tf:
        muffler = pickle.load(tf)
    muffler_contour = muffler['contour']
    for key in muffler_contour.keys():
        generate_ifc_polygon('IfcAirTerminalBox', 'vav_5', muffler_contour[key], 0.0, 500.0, scale,
                             [0.0, 0., 1.0])

    with open("pkl/air_outlet.pkl", "rb") as tf:
        air_outlet = pickle.load(tf)
    air_outlet_contour = air_outlet['contour']
    for key in air_outlet_contour.keys():
        generate_ifc_polygon('IfcFlowTerminal', 'grille_1', air_outlet_contour[key], 0.0, 500.0, scale,
                             [0.0, 1., 1.0])

    with open("pkl/duct1.pkl", "rb") as tf:
        air_outlet = pickle.load(tf)
    air_outlet_contour = air_outlet['contour']
    for key in air_outlet_contour.keys():
        generate_ifc_polygon('IfcFlowSegment', 'duct', air_outlet_contour[key], 200.0, 200.0, scale,
                             [0.0, 0.5, 1.0])

    ifcfile.write(filename)
    end = time.time()
    print('used time: ', end-start)


if __name__ == '__main__':
    main()
