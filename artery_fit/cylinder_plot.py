# python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

name = "data_cylinder3"
#-------------------------------------------------------------------#
def error_curve(ref_data, radius_sim):
    eps = 0.05
    determination_90 = [-1.0]
    rmse_90 = [-1.0]
    error_90 = -1.0

    # assign pressure "0" to prestretch frames. frames of "Step-1"
#    diameter_sim = []
#    for i in range(radius_sim.shape[0]):
#        if radius_sim[i, 0] >= 0.0:
#            diameter_sim.append(radius_sim[i,:])
    diameter_sim = np.copy(radius_sim)
    diameter_sim[:, 0] -= 1.0
    diameter_sim[:, 0] *= 150.0
    diameter_sim[:, 1] *= 2000.0
    diameter_avg = np.mean(ref_data[:,0])

    # check abaqus convergence
    end_pressure = 150.0
    if radius_sim[-1,0] < (2.0-eps):
        end_pressure = diameter_sim[-1,0]

    # estimate error between simulation output and experimental data
    error_sim = 0.0
    count = 0
    diameter_pred = np.zeros((ref_data.shape[0]), dtype=float)
    for i in range(ref_data.shape[0]):
        count = i
        if ref_data[i,1] > end_pressure: break
        if np.min(np.abs(diameter_sim[:,0] - ref_data[i,1])) < 1.0:
            idx = np.argmin(np.abs(diameter_sim[:,0] - ref_data[i,1]))
            diameter_pred[i] = diameter_sim[idx,1]
        else:
            idx = np.argmin(np.abs(diameter_sim[:, 0] - ref_data[i, 1]))
            if diameter_sim[idx, 0] < ref_data[i, 1]:
                diameter_pred[i] = diameter_sim[idx, 1] + (ref_data[i, 1] - diameter_sim[idx, 0])\
                *(diameter_sim[idx+1, 1]-diameter_sim[idx, 1])/(diameter_sim[idx+1, 0]-diameter_sim[idx, 0])
            else:
                diameter_pred[i] = diameter_sim[idx, 1] + (ref_data[i, 1] - diameter_sim[idx, 0])\
                *(diameter_sim[idx, 1]-diameter_sim[idx-1, 1])/(diameter_sim[idx, 0]-diameter_sim[idx-1, 0])
        error_sim += ((ref_data[i, 0] - diameter_pred[i]) / diameter_avg) ** 2

    from sklearn.metrics import r2_score, root_mean_squared_error
    determination = r2_score(ref_data[:count,0], diameter_pred[:count], multioutput='raw_values')
    rmse = root_mean_squared_error(ref_data[:count, 0], diameter_pred[:count], multioutput='raw_values')
    print("\n->R2={}".format(determination))
    print("->RMSE={}".format(rmse))
    print("->error={}\n".format(error_sim))
    if count>=8:
        error_90 = np.abs(ref_data[8,0]-diameter_pred[8])
        determination_90 = r2_score(ref_data[8:count,0], diameter_pred[8:count], multioutput='raw_values')
        rmse_90 = root_mean_squared_error(ref_data[8:count, 0], diameter_pred[8:count], multioutput='raw_values')
        print("->error_90={}-um".format(error_90))
        print("->R2_90={}".format(determination_90))
        print("->RMSE_90={}\n".format(rmse_90))

    return determination[0], rmse[0], error_90, determination_90[0], rmse_90[0]

#-------------------------------------------------------------------#
# reference data for the fitting
data_path = "./"
df = pd.read_excel(data_path + "data_ferruzzi.ods", sheet_name="20_weeks")
df_dta = df.loc[2:15, ['Unnamed: 6', 'Unnamed: 7']]
dta_data = df_dta.to_numpy(dtype=float)

#-------------------------------------------------------------------#
folder = "./"
stress_data = np.loadtxt(folder + name + ".dat", dtype=float)
r2, rmse, eps_90, r2_90, rmse_90 = error_curve(dta_data, stress_data)

pressure_array = 150.0*(stress_data[:, 0] - 1.0)
for i in range(pressure_array.shape[0]):
    if pressure_array[i] < 0.0:
        pressure_array[i] = 0.0

# plot:
plt.rcParams.update({'font.size': 12})
fig, ax = plt.subplots(figsize=plt.figaspect(1.0))
ax.grid()
ax.scatter(dta_data[:,0], dta_data[:,1], s=100, marker='*', c='k')
ax.plot(2.0e3*stress_data[:, 1], pressure_array, 'k-.', lw=2, label='overall')
ax.set_ylim(bottom=-5,top=150)
#ax.set_xlim(left=0.97,right=1.15)
ax.set_xlabel(r'outer diameter [$\mu$m]')
ax.set_ylabel('pressure [mmHg]')

# place a text box in bottom right
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr_r2 = ( f'$R^2$=' f'{r2:.3f}' '\n'
               r'$RMSE$=' f'{rmse:.1f}' r'-$\mu$m' )
# place a text box in upper left in axes coords
plt.text(0.55, 0.15, textstr_r2, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props) #

# place a text box in bottom right
textstr_r2_90 = ( r'$R^2_{ph}$=' f'{r2_90:.3f}' '\n'
                  r'$RMSE_{ph}$=' f'{rmse_90:.1f}' r'-$\mu$m' '\n'
                  r'$\epsilon_{90}$=' f'{eps_90:.1f}' r'-$\mu$m' )
# place a text box in upper left in axes coords
plt.text(0.05, 0.80, textstr_r2_90, transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props) #

fig.tight_layout()
plt.show
FIGURENAME = folder + name + "_pd.pdf"
plt.savefig(FIGURENAME)
plt.close('all')
