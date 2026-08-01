# abaqus python
from odbAccess import openOdb
import numpy as np

step_first = 'Step-1'

#===================================================================#
def stress_strain(odb, name):

    part1 = odb.rootAssembly.instances['GROUND-1']
    #part2 = odb.rootAssembly.instances['CELL-1']
    # label of outer_point node set. to follow the diameter displacement
    outer_point = part1.nodeSets['OUTER_POINT'].nodes[0].label
    inner_point = part1.nodeSets['INNER_POINT'].nodes[0].label
    outer_point_coord = part1.nodeSets['OUTER_POINT'].nodes[0].coordinates[0]
    inner_point_coord = part1.nodeSets['INNER_POINT'].nodes[0].coordinates[0]

    # set a frame for reference
    frame_obj0 = odb.steps[step_first].frames[0]

    # number of integration points in the simulation
    ngauss = len(frame_obj0.fieldOutputs['S'].values)
    ngauss_gnd = len(frame_obj0.fieldOutputs['S'].bulkDataBlocks[0].data)
    ngauss_cell = len(frame_obj0.fieldOutputs['S'].bulkDataBlocks[1].data)
    print("Number of integration points = "+str(ngauss))
    print("Number of integration points (ground) = "+str(ngauss_gnd))
    print("Number of integration points (cells) = "+str(ngauss_cell))

    # number of elements in the simulation
    nelem = len(frame_obj0.fieldOutputs['EVOL'].values)
    nelem_gnd = len(frame_obj0.fieldOutputs['EVOL'].bulkDataBlocks[0].data)
    nelem_cell = len(frame_obj0.fieldOutputs['EVOL'].bulkDataBlocks[1].data)
    print("Number of elements = "+str(nelem))
    print("Number of elements (ground) = "+str(nelem_gnd))
    print("Number of elements (cells) = "+str(nelem_cell))

    # seek the value id for the outer_point label
    nnodes = len(frame_obj0.fieldOutputs['U'].values)
    for i in range(nnodes):
        if (frame_obj0.fieldOutputs['U'].values[i].nodeLabel == outer_point) and \
                (frame_obj0.fieldOutputs['U'].values[i].instance.name == 'GROUND-1'):
            outer_point_id = i
        if (frame_obj0.fieldOutputs['U'].values[i].nodeLabel == inner_point) and \
                (frame_obj0.fieldOutputs['U'].values[i].instance.name == 'GROUND-1'):
            inner_point_id = i

    nframes = 0
    for istep in odb.steps.keys():
        nframes += len(odb.steps[istep].frames)
    time = np.zeros((nframes), dtype=np.float64)

    data_out = np.zeros((25,nframes), dtype=np.float64)
    jframe = 0
    initial_time = 0.0
    for istep in odb.steps.keys():
        for iframe in range(len(odb.steps[istep].frames)):
            # frame values
            frame_obj = odb.steps[istep].frames[iframe]
            time[jframe] = frame_obj.frameValue + initial_time
            # value given by node
            displacement_outer = frame_obj.fieldOutputs['U'].values[outer_point_id].data[0]
            displacement_inner = frame_obj.fieldOutputs['U'].values[inner_point_id].data[0]
            # value given by element in the ground material
            volume = sum(frame_obj.fieldOutputs['EVOL'].bulkDataBlocks[0].data)
            # These values are given by integration point
            stress_vm = sum(frame_obj.fieldOutputs['S'].bulkDataBlocks[0].mises)
            stress_tot = sum(frame_obj.fieldOutputs['SDV19'].bulkDataBlocks[0].data)
            stress_vol = sum(frame_obj.fieldOutputs['SDV20'].bulkDataBlocks[0].data)
            stress_mat = sum(frame_obj.fieldOutputs['SDV21'].bulkDataBlocks[0].data)
            stress_fib = sum(frame_obj.fieldOutputs['SDV22'].bulkDataBlocks[0].data)
            stress_cell = sum(frame_obj.fieldOutputs['S'].bulkDataBlocks[1].data[:,0])
            stress_pas = sum(frame_obj.fieldOutputs['SDV21'].bulkDataBlocks[1].data)
            stress_act = sum(frame_obj.fieldOutputs['SDV22'].bulkDataBlocks[1].data)
            rstretch_11 = sum(frame_obj.fieldOutputs['SDV1'].bulkDataBlocks[0].data)
            rstretch_22 = sum(frame_obj.fieldOutputs['SDV2'].bulkDataBlocks[0].data)
            rstretch_33 = sum(frame_obj.fieldOutputs['SDV3'].bulkDataBlocks[0].data)
            rstretch_fib = sum(frame_obj.fieldOutputs['SDV7'].bulkDataBlocks[0].data)
            rstretch_cell = sum(frame_obj.fieldOutputs['SDV7'].bulkDataBlocks[1].data)
            det_pl8 = sum(frame_obj.fieldOutputs['SDV8'].bulkDataBlocks[0].data)
            det_pl9 = sum(frame_obj.fieldOutputs['SDV9'].bulkDataBlocks[0].data)
            det_ela = sum(frame_obj.fieldOutputs['SDV10'].bulkDataBlocks[0].data)
            det_tot = sum(frame_obj.fieldOutputs['SDV11'].bulkDataBlocks[0].data)
            rho_tot = sum(frame_obj.fieldOutputs['SDV12'].bulkDataBlocks[0].data)
            phi_mat = sum(frame_obj.fieldOutputs['SDV13'].bulkDataBlocks[0].data)
            phi_fib = sum(frame_obj.fieldOutputs['SDV14'].bulkDataBlocks[0].data)
            det_pl9_cell = sum(frame_obj.fieldOutputs['SDV9'].bulkDataBlocks[1].data)
            rho_cell = sum(frame_obj.fieldOutputs['SDV12'].bulkDataBlocks[1].data)

            # collection to export
            data_out[0,jframe] = outer_point_coord + displacement_outer
            data_out[1,jframe] = inner_point_coord + displacement_inner
            data_out[2,jframe] = stress_tot/float(ngauss_gnd)
            data_out[3,jframe] = stress_vol/float(ngauss_gnd)
            data_out[4,jframe] = stress_mat/float(ngauss_gnd)
            data_out[5,jframe] = stress_fib/float(ngauss_gnd)
            data_out[6,jframe] = stress_cell/float(ngauss_cell)
            data_out[7,jframe] = rstretch_11/float(ngauss_gnd)
            data_out[8,jframe] = rstretch_22/float(ngauss_gnd)
            data_out[9,jframe] = rstretch_33/float(ngauss_gnd)
            data_out[10,jframe] = rstretch_fib/float(ngauss_gnd)
            data_out[11,jframe] = rstretch_cell/float(ngauss_cell)
            data_out[12,jframe] = det_tot/float(ngauss_gnd)
            data_out[13,jframe] = det_ela/float(ngauss_gnd)
            data_out[14,jframe] = det_pl9/float(ngauss_gnd)
            data_out[15,jframe] = det_pl8/float(ngauss_gnd)
            data_out[16,jframe] = phi_mat/float(ngauss_gnd)
            data_out[17,jframe] = phi_fib/float(ngauss_gnd)
            data_out[18,jframe] = rho_tot/float(ngauss_gnd)
            data_out[19,jframe] = det_pl9_cell/float(ngauss_cell)
            data_out[20,jframe] = rho_cell/float(ngauss_cell)
            data_out[21,jframe] = stress_vm/float(ngauss_gnd)
            data_out[22,jframe] = stress_pas/float(ngauss_cell)
            data_out[23,jframe] = stress_act/float(ngauss_cell)
            data_out[24,jframe] = volume
            jframe += 1

        initial_time += frame_obj.frameValue

    output_data = np.transpose(np.vstack((time,data_out)))
    np.savetxt('data_'+name+'.dat',output_data,fmt='%f')

#===================================================================#
odb_file = ['cylinder3_2d.odb'] #

for idb in odb_file:
    odb = openOdb(idb)
    name=idb[:9]
    print("Processing simulation case: " + name)
    stress_strain(odb, name)
