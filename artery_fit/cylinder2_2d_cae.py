# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2022 replay file
# Run by Joan Laubrie on Tue Mar 12 11:05 2025
#
from abaqus import *
from abaqusConstants import *
from caeModules import *

import numpy as np
#=============================================================================#
#====================     GEOMETRIC  PARAMETERS        =======================#
#=============================================================================#
length = 5.800
stretch_z = 1.5
n_lamellar_layers = 5

# unloaded geometry
do_ini = 0.827
h_ini = 0.105
di_ini = do_ini - 2.0*h_ini
ratio_media = 0.735
vol_ini = 0.25*np.pi*(do_ini**2 - di_ini**2)*length
lamellar_thickness_ini = 0.002
muscle_thickness_ini = (h_ini*ratio_media-float(n_lamellar_layers+1)*lamellar_thickness_ini)/float(n_lamellar_layers)
print("initial volume: " + str(vol_ini))
print("outer diameter (tf) = " + str(do_ini))
print("inner diameter (tf) = " + str(di_ini))
print("thickness (tf) = " + str(h_ini))
print("lamellar thickness (tf) = " + str(lamellar_thickness_ini))
print("muscle layer thickness (tf) = " + str(muscle_thickness_ini))
#lambda_z = 1.5
#lambda_t = 1.0

# geometry at 40mmHg
do_40 = 0.945
di_40 = np.sqrt(do_40**2 - 4.0*vol_ini/(np.pi*length*stretch_z))
h_40 = 0.5*(do_40 - di_40)
lamellar_thickness_40 = lamellar_thickness_ini*h_40/h_ini
muscle_thickness_40 = muscle_thickness_ini*h_40/h_ini
print("outer diameter (P=40-mmHg) = " + str(do_40))
print("inner diameter (P=40-mmHg) = " + str(di_40))
print("thickness (P=40-mmHg) = " + str(h_40))
print("lamellar thickness (P=40-mmHg) = " + str(lamellar_thickness_40))
print("muscle layer thickness (P=40-mmHg) = " + str(muscle_thickness_40))
#lambda_z = 1.00
#lambda_t = 0.87
#lambda_f = 0.92

# geometry at 90mmHg
do_90 = 1.297
di_90 = np.sqrt(do_90**2 - 4.0*vol_ini/(np.pi*length*stretch_z))
h_90 = 0.5*(do_90 - di_90)
lamellar_thickness_90 = lamellar_thickness_ini*h_90/h_ini
muscle_thickness_90 = muscle_thickness_ini*h_90/h_ini
print("outer diameter (P=90-mmHg) = " + str(do_90))
print("inner diameter (P=90-mmHg) = " + str(di_90))
print("thickness (P=90-mmHg) = " + str(h_90))
print("lamellar thickness (P=90-mmHg) = " + str(lamellar_thickness_90))
print("muscle layer thickness (P=90-mmHg) = " + str(muscle_thickness_90))
lambda_z = 1.00
lambda_t = 0.56
lambda_f = 0.90

#----------------------------------------------------------------#
# overall geometry
r_outer = 0.5*do_90
r_inner = 0.5*di_90
lamellar_thickness = lamellar_thickness_90
muscle_thickness = muscle_thickness_90
angle = np.pi/2.0

# interlayer radius, interface media-adventitia
r_interlayer = np.sqrt(ratio_media*r_outer**2 + (1.0-ratio_media)*r_inner**2)

# thickness by layer
t_artery = r_outer - r_inner
t_media = r_interlayer - r_inner
t_adventitia = r_outer - r_outer

# elastin lamellas
unit_step = (t_media-lamellar_thickness)/n_lamellar_layers      # lamellar unit thickness
r_lamella = np.zeros((n_lamellar_layers+1,2), dtype=np.float64)
r_muscle = np.zeros((n_lamellar_layers), dtype=np.float64)
r_lamella[0,0] = r_inner                       #inner in first lamella
r_lamella[-1,1] = r_outer                      #outer in last lamella
for i in range(n_lamellar_layers):
    # inner radius for a elastin lamellar layer
    r_lamella[i+1,0] = r_inner + float(i+1)*unit_step
    # outer radius for a elastin lamellar layer
    r_lamella[i,1] = r_inner + lamellar_thickness + float(i)*unit_step
    # outer radius for a elastin lamellar layer
    r_muscle[i] = 0.5*(r_lamella[i,1] + r_lamella[i+1,0])

print("Muscle radii = " + str(r_muscle))

#=============================================================================#
#====================     SIMULATION  PARAMETERS        ======================#
#=============================================================================#
# mesh parameters
elem_shape1 = QUAD   # QUAD, TRI
elem_type1 = CPE4      # shell: S4, S3; membrane: M3D4, M3D3
elem_size1 = 0.004

# reinforcement mesh parameters
elem_type2 = T2D3      # 3D continuous elements: T3D2, 2D continuous elements: T2D2, T2D3
elem_size2 = 0.10

# material parameters. total density of tissue = 1050.0 [kg/m3]; stiffness = 100 [N.m/kg]
# media layer. total_density = 1050.0e-9
fraction_mat1 = 0.6
fraction_fib1 = 0.4
# matrix properties
E_gnd1 = 0.25
nu_gnd1 = 0.49
prestretch_circ1 = lambda_t
prestretch_long1 = lambda_z
# fibre properties
k1_fib1 = 0.40
k2_fib1 = 3.5
angle_fib1 = 42.0
prestretch_fib1 = lambda_f
# remodeling properties
tau_fibre1 = 12000.0      # x3 -> 90-days
growth_fibre1 = 0.0

# adventitia layer. total_density = 1050.0e-9
fraction_mat2 = 0.35
fraction_fib2 = 0.65
# matrix properties
E_gnd2 = 0.25
nu_gnd2 = 0.49
prestretch_circ2 = lambda_t
prestretch_long2 = lambda_z
# fibre properties
k1_fib2 = 0.40
k2_fib2 = 3.5
angle_fib2 = 53.0
prestretch_fib2 = lambda_f
# remodeling properties
tau_fibre2 = 12000.0      # x3 -> 90-days
growth_fibre2 = 0.0

# cell
E_cell = 0.10
s_max = 0.10
l_min = 0.65
l_max = 1.40
# remodeling properties
tau_cell = 12000.0        # x3 -> 90-days
growth_cell = 0.0

inner_pressure1 = 0.020                 # lumen pressure
period1 = 1.0                           # prestress step
period2 = 1.0                           # pressure-load step
time1 = period1
time2 = time1 + period2

frame_factor = 1.0
visco_tol = 1.0
#=============================================================================#
#================       CELL GEOMETRY CONSTRUCTION        ====================#
#=============================================================================#
# dimensions of a cell
cell_semilength = 0.020  #*h_ini/h_40
cell_semiwidth = 0.5*muscle_thickness
cell_radius = cell_semiwidth
cell_area = np.pi*cell_radius**2

# cell tilting and number of cells per interlamellar layer
tilt_angle = 0.0*np.pi/180.0
n_cells_per_layer = np.pi*r_muscle/cell_semilength
n_cells_per_layer = n_cells_per_layer.astype(int) #- 5

# define coordinates for the center of the cells
cell_to_cell_angle = 2.0*np.pi/n_cells_per_layer.astype(float)
layer_phase = 0.6*cell_to_cell_angle
theta0 = np.zeros((r_muscle.shape[0],np.max(n_cells_per_layer)), dtype=float)
theta0 += np.nan
for ir in range(r_muscle.shape[0]):
    factor = ir % 2
    for itheta in range(n_cells_per_layer[ir]):
        theta0[ir,itheta] = float(itheta)*cell_to_cell_angle[ir] + factor*layer_phase[ir]

# coordinates of cells extremes
x0_c = np.zeros((r_muscle.shape[0],theta0.shape[1]), dtype=np.float64)
y0_c = np.zeros((r_muscle.shape[0],theta0.shape[1]), dtype=np.float64)
length_array = np.zeros((r_muscle.shape[0],theta0.shape[1],2), dtype=np.float64)
for ir in range(r_muscle.shape[0]):
    for itheta in range(theta0.shape[1]):
        # cell center
        x0_c[ir,itheta] = r_muscle[ir]*np.cos(theta0[ir,itheta])
        y0_c[ir,itheta] = r_muscle[ir]*np.sin(theta0[ir,itheta])
        # cell semi-length
        length_array[ir,itheta,0] = cell_semilength*np.sin(theta0[ir,itheta] + tilt_angle)
        length_array[ir,itheta,1] = cell_semilength*np.cos(theta0[ir,itheta] + tilt_angle)

#=============================================================================#
#========================      MODEL CREATION       ==========================#
#=============================================================================#
# change name to the model
#mdb.models.changeKey(fromName='Model-1',toName='valve')
mdb.Model(name='aortic_wall')
model1 = mdb.models['aortic_wall']

#=============================================================================#
#======================     ARTERIAL WALL GEOMETRY     =======================#
#=============================================================================#
# initial center plane. feature #1
model1.ConstrainedSketch(name='__profile__', sheetSize=1.0)
model1.sketches['__profile__'].CircleByCenterPerimeter(center=(0.0, 0.0), point1=(r_outer, 0.0))
model1.sketches['__profile__'].CircleByCenterPerimeter(center=(0.0, 0.0), point1=(r_inner, 0.0))
part1 = model1.Part(dimensionality=TWO_D_PLANAR, name='ground', type=DEFORMABLE_BODY)
part1.BaseShell(sketch=model1.sketches['__profile__'])
del model1.sketches['__profile__']

# split the arteial wall into layers
print("split the wall into the layers")
f0, e0 = part1.faces, part1.edges
print(f0)
for i in f0:
    print(i)
print(e0)
for i in e0:
    print(i)
t0 = part1.MakeSketchTransform(sketchPlane=f0[0], sketchUpEdge=e0[0],
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s0 = model1.ConstrainedSketch(name='__profile__',
    sheetSize=0.5, gridSpacing=0.1, transform=t0)
s0.setPrimaryObject(option=SUPERIMPOSE)
part1.projectReferencesOntoSketch(sketch=s0, filter=COPLANAR_EDGES)
s0.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(r_interlayer, 0.0))
pickedFaces = part1.faces.getByBoundingBox(
                xMin=-1.1*r_outer, yMin=-1.1*r_outer, zMin=-1.0,
                xMax=1.1*r_outer, yMax=1.1*r_outer, zMax=1.0)
part1.PartitionFaceBySketch(sketchUpEdge=part1.edges[0], faces=pickedFaces, sketch=s0)
s0.unsetPrimaryObject()
del model1.sketches['__profile__']

# take just the portion of circumference of interest
# draw the lines to cut the geometry
print("cut the wall to the region of interest")
f0, e0 = part1.faces, part1.edges
print(f0)
for i in f0:
    print(i)
print(e0)
for i in e0:
    print(i)
t0 = part1.MakeSketchTransform(sketchPlane=f0[0], sketchUpEdge=e0[1],
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
s0 = model1.ConstrainedSketch(name='__profile__', sheetSize=2.0, 
    gridSpacing=0.1, transform=t0)
s0.setPrimaryObject(option=SUPERIMPOSE)
part1.projectReferencesOntoSketch(sketch=s0, filter=COPLANAR_EDGES)
s0.Line(point1=(0.0, 0.0), point2=(1.1*r_outer, 0.0))
s0.Line(point1=(0.0, 0.0), point2=(1.1*r_outer*np.cos(angle), 1.1*r_outer*np.sin(angle)))
pickedFaces = part1.faces.getByBoundingBox(
                xMin=-1.1*r_outer, yMin=-1.1*r_outer, zMin=-1.0, 
                xMax=1.1*r_outer, yMax=1.1*r_outer, zMax=1.0)
part1.PartitionFaceBySketch(sketchUpEdge=part1.edges[1], faces=pickedFaces, sketch=s0)
s0.unsetPrimaryObject()
del model1.sketches['__profile__']
# make a list with the faces to remove
all_faces = part1.faces.getByBoundingBox(
                xMin=-1.1*r_outer, yMin=-1.1*r_outer, zMin=-1.0, 
                xMax=1.1*r_outer, yMax=1.1*r_outer, zMax=1.0)
not_picked_faces = part1.faces.getByBoundingBox(
                xMin=-0.1*r_outer, yMin=-0.1*r_outer, zMin=-1.0, 
                xMax=1.1*r_outer, yMax=1.1*r_outer*np.sin(angle), zMax=1.0)
not_picked_faces_index = []
for iface in not_picked_faces:
    not_picked_faces_index.append(iface.index)
picked_faces = []
for iface in all_faces:
    if not iface.index in not_picked_faces_index:
        picked_faces.append(iface)
part1.RemoveFaces(faceList = picked_faces, deleteCells=False)

#----------------------------------------------------------------#
# update validity
part1.checkGeometry()

#----------------------------------------------------------------#
# sets for the whole solid
all_faces = part1.faces.getByBoundingCylinder(
                center1=(0,0,-1.0), center2=(0,0,1.0), radius=1.001*r_outer)
# sets for the media layer
media_faces = part1.faces.getByBoundingCylinder(
                center1=(0,0,-1.0), center2=(0,0,1.0), radius=1.001*r_interlayer)

# body set
part1.Set(faces=all_faces, name='body')
# media layer set
part1.Set(faces=media_faces, name='media')
# adventitia layer set
part1.SetByBoolean(
        name='adventitia',
        operation=DIFFERENCE,
        sets=(
            part1.sets['body'],
            part1.sets['media'], )
        )

#----------------------------------------------------------------#
all_edges = part1.edges.getByBoundingCylinder(
                center1=(0,0,-1.0), center2=(0,0,1.0), radius=1.001*r_outer)
not_outer_edges = part1.edges.getByBoundingCylinder(
                center1=(0,0,-1.0), center2=(0,0,1.0), radius=1.001*r_interlayer)
inner_edge = part1.edges.getByBoundingCylinder(
                center1=(0,0,-1.0), center2=(0,0,1.0), radius=1.001*r_inner)
top_edge = part1.edges.getByBoundingBox(
                xMin=0.9*r_inner*np.cos(angle), yMin=0.9*r_inner*np.sin(angle), zMin=-1.0, 
                xMax=1.1*r_outer*np.cos(angle), yMax=1.1*r_outer*np.sin(angle), zMax=1.0)
bottom_edge = part1.edges.getByBoundingBox(
                xMin=0.9*r_inner, yMin=-0.1*r_inner, zMin=-1.0, 
                xMax=1.1*r_outer, yMax=0.1*r_inner, zMax=1.0)

# sets for the borders
# all_edge set
part1.Set(edges=all_edges, name='all_edges')
# not_outer_edge set
part1.Set(edges=not_outer_edges, name='not_outer_edges')
# top_edge set
part1.Set(edges=top_edge, name='top_edge')
# bottom_edge set
part1.Set(edges=bottom_edge, name='bottom_edge')
# inner set
part1.Set(edges=inner_edge, name='inner_edge')
# interlayer edge set
part1.SetByBoolean(
        name='interlayer_edge',
        operation=DIFFERENCE,
        sets=(
            part1.sets['not_outer_edges'],
            part1.sets['inner_edge'],
            part1.sets['top_edge'],
            part1.sets['bottom_edge'], )
        )
# outer edge set
part1.SetByBoolean(
        name='outer_edge',
        operation=DIFFERENCE, 
        sets=(
            part1.sets['all_edges'], 
            part1.sets['interlayer_edge'],
            part1.sets['inner_edge'],
            part1.sets['top_edge'],
            part1.sets['bottom_edge'], )
        )

#----------------------------------------------------------------#
# surface
part1.Surface(side1Edges=all_edges, name='all_surfs')
part1.Surface(side1Edges=not_outer_edges, name='not_outer_surfs')
part1.Surface(side1Edges=top_edge, name='top_surf')
part1.Surface(side1Edges=bottom_edge, name='bottom_surf')
part1.Surface(side1Edges=inner_edge, name='inner_surf')
part1.SurfaceByBoolean(
        name='interlayer_surf',
        operation=DIFFERENCE,
        surfaces=(
            part1.surfaces['not_outer_surfs'],
            part1.surfaces['inner_surf'],
            part1.surfaces['top_surf'],
            part1.surfaces['bottom_surf'], )
        )
part1.SurfaceByBoolean(
        name='outer_surf',
        operation=DIFFERENCE, 
        surfaces=(
            part1.surfaces['all_surfs'],
            part1.surfaces['interlayer_surf'],
            part1.surfaces['inner_surf'],
            part1.surfaces['top_surf'],
            part1.surfaces['bottom_surf'], )
        )

#----------------------------------------------------------------#
outer_point = part1.vertices.getByBoundingBox(
                xMin=0.99*r_outer, yMin=-0.01*r_outer, zMin=-1.0,
                xMax=1.01*r_outer, yMax=0.01*r_outer, zMax=1.0)
inner_point = part1.vertices.getByBoundingBox(
                xMin=0.99*r_inner, yMin=-0.01*r_inner, zMin=-1.0,
                xMax=1.01*r_inner, yMax=0.01*r_inner, zMax=1.0)

# outer_point set
part1.Set(vertices=outer_point, name='outer_point')
# inner_point set
part1.Set(vertices=inner_point, name='inner_point')

#----------------------------------------------------------------#
# material orientation
part1.MaterialOrientation(
        region=part1.sets['body'], orientationType=DISCRETE, axis=AXIS_3,
        normalAxisDefinition=VECTOR, normalAxisVector=(0.0, 0.0, 1.0),
        flipNormalDirection=False, normalAxisDirection=AXIS_3,
        primaryAxisDefinition=EDGE, primaryAxisRegion=part1.sets['inner_edge'],
        primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)

#=============================================================================#
#=====================     REINFORCEMENT GEOMETRY     ========================#
#=============================================================================#
# initial center plane. feature #1
sketch2 = model1.ConstrainedSketch(name='__profile__', sheetSize=1.0)
sketch2.setPrimaryObject(option=STANDALONE)
for ir in range(r_muscle.shape[0]):
    for itheta in range(theta0.shape[1]):
        if np.isnan(theta0[ir,itheta]): continue
        sketch2.Line(
            point1=(x0_c[ir,itheta]+length_array[ir,itheta,0], y0_c[ir,itheta]-length_array[ir,itheta,1]),
            point2=(x0_c[ir,itheta]-length_array[ir,itheta,0], y0_c[ir,itheta]+length_array[ir,itheta,1]) )
part2 = model1.Part(name='cell', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
part2.BaseWire(sketch=sketch2)
sketch2.unsetPrimaryObject()
del model1.sketches['__profile__']

#----------------------------------------------------------------#
# cells are made all along the circumference initially, then we have to cut
print("cut the region of interest for the cells ...")
part2.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
part2.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
pickedEdges = part2.edges.getByBoundingCylinder(
                center1=(0,0,-1.0), center2=(0,0,1.0), radius=1.001*r_outer)
part2.PartitionEdgeByDatumPlane(datumPlane=part2.datums[3], edges=pickedEdges)
pickedEdges = part2.edges.getByBoundingCylinder(
                center1=(0,0,-1.0), center2=(0,0,1.0), radius=1.001*r_outer)
part2.PartitionEdgeByDatumPlane(datumPlane=part2.datums[2], edges=pickedEdges)

# make a list with the faces to remove
all_edges = part2.edges
not_picked_edges = part2.edges.getByBoundingBox(
        xMin=-0.001*r_inner, yMin=-0.001*r_inner, zMin=-0.1,
        xMax= 1.001*r_outer, yMax= 1.001*r_outer, zMax= 0.1)
not_picked_edges_index = []
for iedge in not_picked_edges:
    not_picked_edges_index.append(iedge.index)
picked_edges = []
for iedge in all_edges:
    if not iedge.index in not_picked_edges_index:
        picked_edges.append(iedge)
part2.RemoveWireEdges(wireEdgeList=picked_edges)

#----------------------------------------------------------------#
# update validity
part2.checkGeometry()

#----------------------------------------------------------------#
# sets for the solid
verts2, edges2 = part2.vertices, part2.edges
part2.Set(edges=edges2, name='body')
part2.Set(vertices=verts2, name='points')

#=============================================================================#
#========================       MESH        ==================================#
#=============================================================================#
# meshing
print("building a mesh for the geometry")
part1.setMeshControls(
        regions=part1.faces.getByBoundingCylinder(
                center1=(0,0,-1.0), center2=(0,0,1.0), radius=1.001*r_outer), 
        elemShape=elem_shape1)
part1.seedPart(size=elem_size1, deviationFactor=0.1, minSizeFactor=0.1)
part1.generateMesh()
elemType1 = mesh.ElemType(elemCode=elem_type1, elemLibrary=STANDARD, secondOrderAccuracy=OFF)
part1.setElementType(regions=part1.sets['body'], elemTypes=(elemType1, ))

#-----------------------------------------------------------------------------#
# truss meshing
part2.seedPart(size=elem_size2, deviationFactor=0.5, minSizeFactor=0.5)
part2.generateMesh()
elemType2 = mesh.ElemType(elemCode=elem_type2, elemLibrary=STANDARD)
part2.setElementType(regions=part2.sets['body'], elemTypes=(elemType2, ))

#=============================================================================#
#========================       MATERIALS        =============================#
#=============================================================================#
print("assigning a material to the object")
# material host
material1_name = 'ground_media'
model1.Material(name=material1_name)
material1 = model1.materials[material1_name]
material1.UserMaterial(mechanicalConstants=(E_gnd1, nu_gnd1, k1_fib1, k2_fib1, angle_fib1,
                        prestretch_circ1, prestretch_long1, prestretch_fib1,
                        fraction_mat1, fraction_fib1, tau_fibre1, growth_fibre1))
material1.Depvar(n=23)

# section membrane
model1.HomogeneousSolidSection(name='media_section', material=material1_name)

# section assignment, BOTTOM_SURFACE, TOP_SURFACE
part1.SectionAssignment(region=part1.sets['media'], sectionName='media_section',
        offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

#-----------------------------------------------------------------------------#
material2_name = 'ground_adven'
model1.Material(name=material2_name)
material2 = model1.materials[material2_name]
material2.UserMaterial(mechanicalConstants=(E_gnd2, nu_gnd2, k1_fib2, k2_fib2, angle_fib2,
                        prestretch_circ2, prestretch_long2, prestretch_fib2,
                        fraction_mat2, fraction_fib2, tau_fibre2, growth_fibre2))
material2.Depvar(n=23)

# section membrane
model1.HomogeneousSolidSection(name='adven_section', material=material2_name)

# section assignment, BOTTOM_SURFACE, TOP_SURFACE
part1.SectionAssignment(region=part1.sets['adventitia'], sectionName='adven_section',
        offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
        thicknessAssignment=FROM_SECTION)

#-----------------------------------------------------------------------------#
# material guest
name_guest = 'cell'
model1.Material(name=name_guest)
material5 = model1.materials[name_guest]
material5.UserMaterial(mechanicalConstants=(E_cell, s_max, l_min, l_max, tau_cell, growth_cell))
material5.Depvar(n=23)

model1.TrussSection(name='section_cell', material=name_guest, area=cell_area)

# section assignment, BOTTOM_SURFACE, TOP_SURFACE
part2.SectionAssignment(region=part2.sets['body'], sectionName='section_cell',
        offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
        thicknessAssignment=FROM_SECTION)

#=============================================================================#
#===================           ASSEMBLY               ========================#
#=============================================================================#
# assembly
print("making the assembly of the problem")
assembly = model1.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name='ground-1', part=part1, dependent=ON)
assembly.Instance(name='cell-1', part=part2, dependent=ON)

# creation of a cylindrical coordinates system. more appropriate for BC
assembly.DatumCsysByThreePoints(name='sys_cylindrical', coordSysType=CYLINDRICAL, 
    origin=(0.0, 0.0, 0.0), point1=(1.0, 0.0, 0.0), point2=(0.0, 1.0, 0.0))
print(assembly.datums)
for i in assembly.datums.keys():
    print(assembly.datums[i])
sys_cyl = assembly.datums[6]

#=============================================================================#
#===================          INTERACTION             ========================#
#=============================================================================#
model1.EmbeddedRegion(
        name='inclusion_interaction',
        embeddedRegion=assembly.instances['cell-1'].sets['body'],
        hostRegion=assembly.instances['ground-1'].sets['body'],
        weightFactorTolerance=1e-06, absoluteTolerance=0.0,
        fractionalTolerance=0.05, toleranceMethod=BOTH)

#=============================================================================#
#=================       STEPS AND OUTPUT REQUESTS        ====================#
#=============================================================================#
# steps
initial_inc1 = period1/20.0
max_inc1 = period1/4.0
model1.ViscoStep(name='Step-1', previous='Initial', timePeriod=period1,
                  timeIncrementationMethod=AUTOMATIC,
                  stabilizationMethod=DISSIPATED_ENERGY_FRACTION, adaptiveDampingRatio = 0.05,
                  cetol=visco_tol, integration=EXPLICIT_ONLY,
                  maxNumInc=100, initialInc=initial_inc1, minInc=1e-5, maxInc=max_inc1,
                  nlgeom=ON)
initial_inc2 = period2/20.0
max_inc2 = period2/4.0
model1.ViscoStep(name='Step-2', previous='Step-1', timePeriod=period2,
                  timeIncrementationMethod=AUTOMATIC,
                  stabilizationMethod=DISSIPATED_ENERGY_FRACTION, adaptiveDampingRatio = 0.05,
                  cetol=visco_tol, integration=EXPLICIT_ONLY,
                  maxNumInc=100, initialInc=initial_inc2, minInc=1e-5, maxInc=max_inc2,
                  nlgeom=ON)

# field output for both bodies
model1.fieldOutputRequests['F-Output-1'].setValues(
        frequency=1, variables = ('S', 'NE', 'LE', 'SDV', 'U', 'RF', 'EVOL'),
        region=model1.rootAssembly.allInstances['ground-1'].sets['body'])
model1.FieldOutputRequest(name='F-Output-2', createStepName='Step-1')
model1.fieldOutputRequests['F-Output-2'].setValues(
        frequency=1, variables = ('S', 'NE', 'LE', 'SDV', 'U', 'RF', 'EVOL'),
        region=model1.rootAssembly.allInstances['cell-1'].sets['body'])

# history output for SAMPLE (GROUND)
del model1.historyOutputRequests['H-Output-1']
model1.HistoryOutputRequest(name='H-Output-1g', createStepName='Step-1',
                            region=model1.rootAssembly.allInstances['ground-1'].sets['body'],
                            variables=PRESELECT, numIntervals=20)
model1.historyOutputRequests['H-Output-1g'].deactivate('Step-2')
model1.HistoryOutputRequest(name='H-Output-2g', createStepName='Step-2',
                            region=model1.rootAssembly.allInstances['ground-1'].sets['body'],
                            variables=PRESELECT, numIntervals=20)

# history output for SAMPLE (CELL)
model1.HistoryOutputRequest(name='H-Output-1c', createStepName='Step-1',
                            region=model1.rootAssembly.allInstances['cell-1'].sets['body'],
                            variables=PRESELECT, numIntervals=20)
model1.historyOutputRequests['H-Output-1c'].deactivate('Step-2')
model1.HistoryOutputRequest(name='H-Output-2c', createStepName='Step-2',
                            region=model1.rootAssembly.allInstances['cell-1'].sets['body'],
                            variables=PRESELECT, numIntervals=20)

#=============================================================================#
#================       BOUNDARY CONDITIONS AND LOADS        =================#
#=============================================================================#
# amplitude
model1.TabularAmplitude(name='load_amp1', timeSpan=TOTAL,
        smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (time1, 0.0), (time2, 1.0) ))

# boundary conditions
model1.YsymmBC(name='top_edge_bc', createStepName='Initial',
        region=assembly.instances['ground-1'].sets['top_edge'],
        localCsys=sys_cyl)
model1.YsymmBC(name='bottom_edge_bc', createStepName='Initial',
        region=assembly.instances['ground-1'].sets['bottom_edge'],
        localCsys=sys_cyl)

# loads
model1.Pressure(name='pressure1', createStepName='Step-2',
        region=assembly.instances['ground-1'].surfaces['inner_surf'],
        distributionType=UNIFORM, field='', magnitude=inner_pressure1, amplitude='load_amp1')

#=============================================================================#
#======================       JOB AND OUTPUT        ==========================#
#=============================================================================#

# create job
mdb.Job(name='cylinder2_2d', model='aortic_wall', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
        userSubroutine='', 
        scratch='', resultsFormat=ODB, 
        multiprocessingMode=THREADS, numCpus=1, numDomains=1, numGPUs=0, 
        parallelizationMethodExplicit=DOMAIN)

# add lines into the input file
model1.keywordBlock.synchVersions(storeNodesAndElements=False)
model1.keywordBlock.insert(88, """** \n** PREDEFINED FIELDS\n** \n** Name: initial_stretch Type: Solution""")
model1.keywordBlock.insert(89, """*Initial Conditions, type=SOLUTION, user""")

# write input file
mdb.jobs['cylinder2_2d'].writeInput()


