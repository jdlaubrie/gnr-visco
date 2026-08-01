# abaqus cae noGUI=test_2d_cae.py
#
# Abaqus/CAE Release 2022 replay file
# Run by Joan Laubrie on Mon Mar  10 12:30 2025
#
from abaqus import *
from abaqusConstants import *
from caeModules import *

#=============================================================================#
#========================       PARAMETERS        ============================#
#=============================================================================#
# geometry
a_side = 1.0
b_side = 1.0

# mesh parameters
elem_shape = QUAD   # QUAD, TRI
elem_type = CPE4      # shell: S4, S3; membrane: M3D4, M3D3; plane strain: CPE3, CPE4(S); plane stress: CPS3, CPS4(S)
elem_size = 0.5

# material parameters. total density of tissue = 1050.0 [kg/m3]
#total_density = 1050.0e-9
nu = 0.495
# ground substance
E_gnd = 0.220

disp_x = 0.5*a_side
period = 1.0

#=============================================================================#
#========================      MODEL CREATION       ==========================#
#=============================================================================#
# change name to the model
#mdb.models.changeKey(fromName='Model-1',toName='valve')
mdb.Model(name='test')
model1 = mdb.models['test']

#=============================================================================#
#===========================        GEOMETRY       ===========================#
#=============================================================================#
# initial center plane. feature #1
model1.ConstrainedSketch(name='__profile__', sheetSize=10.0)
model1.sketches['__profile__'].rectangle(point1=(0.0, 0.0), point2=(a_side, b_side))
part1 = model1.Part(dimensionality=TWO_D_PLANAR, name='sample', type=DEFORMABLE_BODY)
#part1 = model1.Part(dimensionality=THREE_D, name='sample', type=DEFORMABLE_BODY)
part1.BaseShell(sketch=model1.sketches['__profile__'])
del model1.sketches['__profile__']

# update validity
part1.checkGeometry()

# body set
part1.Set(faces=part1.faces.getByBoundingBox(
                xMin=-1.0, yMin=-1.0, zMin=-1.0, 
                xMax=a_side+1.0, yMax=b_side+1.0, zMax=1.0), 
          name='body')
# bottom edge set
part1.Set(edges=part1.edges.getByBoundingBox(
                xMin=-0.01*a_side, yMin=-0.01*b_side, zMin=-0.1,
                xMax=1.01*a_side, yMax=0.01*b_side, zMax=0.1),
          name='bottom')
# top edge set
part1.Set(edges=part1.edges.getByBoundingBox(
                xMin=-0.01*a_side, yMin=0.99*b_side, zMin=-0.1,
                xMax=1.01*a_side, yMax=1.01*b_side, zMax=0.1),
          name='top')
# left edge set
part1.Set(edges=part1.edges.getByBoundingBox(
                xMin=-0.01*a_side, yMin=-0.01*b_side, zMin=-0.1,
                xMax=0.01*a_side, yMax=1.01*b_side, zMax=0.1),
          name='left')
# right edge set
part1.Set(edges=part1.edges.getByBoundingBox(
                xMin=0.99*a_side, yMin=-0.01*b_side, zMin=-0.1,
                xMax=1.01*a_side, yMax=1.01*b_side, zMax=0.1),
          name='right')
# surface
#part1.Surface(side1Faces=part1.faces.getByBoundingBox(
#                xMin=-1.0, yMin=-1.0, zMin=-1.0,
#                xMax=a_side+1.0, yMax=b_side+1.0, zMax=1.0),
#              name='surf')

# material orientation
part1.MaterialOrientation(
        region=part1.sets['body'], orientationType=DISCRETE, axis=AXIS_3,
        normalAxisDefinition=VECTOR, normalAxisVector=(0.0, 0.0, 1.0),
        flipNormalDirection=False, normalAxisDirection=AXIS_3,
        primaryAxisDefinition=EDGE, primaryAxisRegion=part1.sets['bottom'],
        primaryAxisDirection=AXIS_1, flipPrimaryDirection=False)
#
#=============================================================================#
#========================       MESH        ==================================#
#=============================================================================#
# meshing
part1.setMeshControls(
        regions=part1.faces.getByBoundingBox(
            xMin=-1.0, yMin=-1.0, zMin=-1.0, 
            xMax=a_side+1.0, yMax=b_side+1.0, zMax=1.0), 
        elemShape=elem_shape)
part1.seedPart(size=elem_size, deviationFactor=0.5, minSizeFactor=0.5)
part1.generateMesh()
elemType1 = mesh.ElemType(elemCode=elem_type, elemLibrary=STANDARD, secondOrderAccuracy=OFF)
part1.setElementType(regions=part1.sets['body'], elemTypes=(elemType1, ))

#=============================================================================#
#========================       MATERIALS        =============================#
#=============================================================================#
# material
material_name = 'ground'
model1.Material(name=material_name)
material1 = model1.materials[material_name]
material1.Density(table=((1050.0e-09, ), ))
material1.UserMaterial(mechanicalConstants=(E_gnd, nu))
#material1.UserMaterial(mechanicalConstants=(E_gnd, nu, k1, k2, theta_c, dispe_c, stretch_t, stretch_z, stretch_f, 1.0, 1.0))
#material1.UserMaterial(mechanicalConstants=(E_myo, stretch_m, s_max, lambda_min, lambda_max))
#material1.UserMaterial(mechanicalConstants=(a_gs, b_gs, ratio_gs))
#material1.Elastic(type=ISOTROPIC, table=((elasticity, poisson),))
# C10, C01, C20, C11, C02
#material1.Hyperelastic(
#        materialType=ISOTROPIC,
#        testData=OFF,
#        type=NEO_HOOKE,
#        #type=POLYNOMIAL,
#        volumetricResponse=VOLUMETRIC_DATA,
#        n=2,
#        table = ((0.42, 1.0/(0.42*24.6)),) )  #intima
#        #table = ((0.58, 1.0/(0.58*24.6)),) ) #adventitia
#        #table = ((-1.1373, 1.206, 6.5364, -17.819, 12.870, 0.2911, 0.0), ) ) # adventitia
#        #table = ((-0.7699, 0.8235, 2.623, -7.5097, 5.8136, 0.3731, 0.0),) ) # intima
#material1.Damping(alpha=1.0)
material1.Depvar(n=4)

model1.HomogeneousSolidSection(name='section_1', material=material_name)

# section assignment, BOTTOM_SURFACE, TOP_SURFACE
part1.SectionAssignment(region=part1.sets['body'], sectionName='section_1', 
        offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

#=============================================================================#
#===================    ASSEMBLY AND INTERACTION      ========================#
#=============================================================================#
# assembly
assembly = model1.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name='sample-1', part=part1, dependent=ON)

#=============================================================================#
#======================       STEPS AND LOADS        =========================#
#=============================================================================#
# steps
initial_inc = period/20.0
max_inc = period/4.0
#model1.ImplicitDynamicsStep(name='Step-1', previous='Initial', timePeriod=period,
#        nlgeom=ON, timeIncrementationMethod=AUTOMATIC,
#        application=QUASI_STATIC,
#        maxNumInc=100, initialInc=initial_inc, minInc=1e-5,maxInc=max_inc)
model1.StaticStep(name='Step-1', previous='Initial', timePeriod=period,
        nlgeom=ON, timeIncrementationMethod=AUTOMATIC,
        maxNumInc=100, initialInc=initial_inc, minInc=1e-5, maxInc=max_inc)
#model1.StaticStep(name='Step-2', previous='Step-1', timePeriod=1.0,
#        nlgeom=ON, timeIncrementationMethod=AUTOMATIC,
#        maxNumInc=100, initialInc=0.05, minInc=1e-5, maxInc=0.5)

# amplitude
model1.TabularAmplitude(name='load_amp', timeSpan=STEP, 
        smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (period, 1.0)))

# boundary conditions
model1.XsymmBC(name='left_bc', createStepName='Initial',
        region=assembly.instances['sample-1'].sets['left'], localCsys=None)
model1.YsymmBC(name='bottom_bc', createStepName='Initial',
        region=assembly.instances['sample-1'].sets['bottom'], localCsys=None)

# loads
model1.DisplacementBC(name='load1_bc', createStepName='Step-1',
        region=assembly.instances['sample-1'].sets['right'],
        u1=disp_x, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
        amplitude='load_amp', fixed=OFF, distributionType=UNIFORM, fieldName='',
        localCsys=None)

#model1.boundaryConditions['load1_bc'].deactivate(stepName='Step-2')
#model1.DisplacementBC(name='load2_bc', createStepName='Step-2',
#        region=assembly.instances['sample-1'].sets['right'],
#        u1=disp_x, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET,
#        amplitude='load_amp', fixed=OFF, distributionType=UNIFORM, fieldName='',
#        localCsys=None)

#=============================================================================#
#======================       JOB AND OUTPUT        ==========================#
#=============================================================================#
# output
model1.fieldOutputRequests['F-Output-1'].setValues(
        #frequency=1, variables = ('S', 'NE', 'LE', 'SDV', 'U', 'RF', 'EVOL', 'A', 'V'))
        frequency = 1, variables = ('S', 'NE', 'LE', 'SDV', 'U', 'RF', 'EVOL'))
model1.historyOutputRequests['H-Output-1'].setValues(numIntervals=20)

# create job
mdb.Job(name='test_2d', model='test', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
        #userSubroutine='exp_linear_2d.f',
        scratch='', resultsFormat=ODB, 
        multiprocessingMode=THREADS, numCpus=1, numDomains=1, numGPUs=0)

# write input file
mdb.jobs['test_2d'].writeInput()
