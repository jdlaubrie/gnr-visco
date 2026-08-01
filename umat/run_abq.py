# python3 run_abq.py
#
import os
import time

# =============================================================================#
# ==============       ABAQUS COMMAND EXECUTION          ======================#
# =============================================================================#
def abq_execution(job_name, u_routine, mode='analysis', cpus=1, mp_mode='threads'):
    # mode = {'analysis', 'interactive', 'datacheck'}
    # processors distribution load
    parallel = ' cpus=' + str(cpus) + ' mp_mode=' + mp_mode

    print("Executing job= " + job_name)
    command = 'abaqus job=' + job_name + ' user=' + u_routine + ' double ' + mode + parallel + ' ask_delete=OFF'
    os.system(command)  # OPERATION TO ABAQUS
    time.sleep(5)
    # call the script for postprocessing
    file_path = job_name + '.lck'
    while os.path.exists(file_path):
        print("still computing. wait 5 more seconds ...")
        time.sleep(5)
    print("ABAQUS simulation finished.")

    #print("Simulation postprocessing. Collecting outer radius ...")
    #os.system(
    #    'abaqus python outer_radius_odb.py ' + job_name)  # OPERATION TO ABAQUS PYTHON
    #radius_data = np.loadtxt(job_name + "_pd.dat", dtype=float)

#=============================================================================#
#==============       SPECIAL INPUT FILE EDITING        ======================#
#=============================================================================#
case_name = 'test_2d'
user_routine = 'imp_neo_hooke_2d.f'
# processors distribution load
cpus = 1

abq_execution(case_name, user_routine, cpus=cpus)