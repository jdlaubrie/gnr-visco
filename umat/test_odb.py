# abaqus python test_odb.py
#
from odbAccess import openOdb
import numpy as np

#===================================================================#
def stress_strain(odb, name):

    nframes = 0
    for istep in odb.steps.keys():
        nframes += len(odb.steps[istep].frames)
    time = np.zeros((nframes), dtype=np.float64)
    data_out = np.zeros((5,nframes), dtype=np.float64)
    jframe = 0
    initial_time = 0.0
    for istep in odb.steps.keys():
        for iframe in range(len(odb.steps[istep].frames)):
            frame_obj = odb.steps[istep].frames[iframe]
            integration_points = len(frame_obj.fieldOutputs['S'].values)
            time[jframe] = frame_obj.frameValue + initial_time
            stress_11 = 0.0
            stress_22 = 0.0
            strain_11 = 0.0
            strain_22 = 0.0
            # value given by element
            volume = frame_obj.fieldOutputs['EVOL'].values[0].data
            # These values are given by integration point
            for i in range(integration_points):
                stress_11 += frame_obj.fieldOutputs['S'].values[i].data[0]
                stress_22 += frame_obj.fieldOutputs['S'].values[i].data[1]
                strain_11 += frame_obj.fieldOutputs['NE'].values[i].data[0]
                strain_22 += frame_obj.fieldOutputs['NE'].values[i].data[1]
            data_out[0,jframe] = stress_11/float(integration_points)
            data_out[1,jframe] = stress_22/float(integration_points)
            data_out[2,jframe] = strain_11/float(integration_points)
            data_out[3,jframe] = strain_22/float(integration_points)
            data_out[4,jframe] = volume
            jframe += 1

        initial_time = frame_obj.frameValue

    output_data = np.transpose(np.vstack((time,data_out)))
    np.savetxt('data_'+name+'.dat',output_data,fmt='%f')

#===================================================================#
odb_file = ['test_2d.odb'] #

for idb in odb_file:
    odb = openOdb(idb)
    name=idb[:-4]
    print("Processing simulation case: " + name)
    stress_strain(odb, name)
