# python3
import numpy as np
import matplotlib.pyplot as plt

name = "cylinder3"
#-------------------------------------------------------------------#
data = np.loadtxt("data_" + name + ".dat", dtype=float)
u_outer_error = 100.0*(data[20, 1]-data[0, 1])/data[0, 1]
u_inner_error = 100.0*(data[20, 2]-data[0, 2])/data[0, 2]
print("displacement outer-reference = {}-%".format(u_outer_error))
print("displacement inner-reference = {}-%".format(u_inner_error))

#-------------------------------------------------------------------#
plt.rcParams.update({'font.size':12})
fig, ax = plt.subplots(2, 3, figsize=(15, 9)) #,
fig.subplots_adjust(wspace=0.25, hspace=0.25)

# time-stress plot
ax[0,0].plot(data[:, 0], data[:, 3]*1.0e3, lw=1, label=r'$p^{mat}$')
ax[0,0].plot(data[:, 0], data[:, 4]*1.0e3, lw=1, label=r'$p^{fib}$')
ax[0,0].plot(data[:, 0], data[:, 5]*1.0e3, lw=1, label=r'$s^{mat}$')
ax[0,0].plot(data[:, 0], data[:, 6]*1.0e3, lw=1, label=r'$s^{fib}$')
ax[0,0].plot(data[:, 0], data[:, 22]*1.0e3, lw=1, label=r'$\sigma_{VM}$')
ax[0,0].axhline(y=data[20, 4]*1.0e3, color='k', linestyle='--', lw=0.5)
ax[0,0].axhline(y=data[20, 5]*1.0e3, color='k', linestyle='--', lw=0.5)
ax[0,0].axhline(y=data[20, 6]*1.0e3, color='k', linestyle='--', lw=0.5)
ax[0,0].axvline(x=data[20, 0], color='k', linestyle='--', lw=0.5)
#ax[0,0].axvline(x=data[72, 0], color='k', linestyle='--', lw=0.5)
#ax[0,0].axvline(x=data[92, 0], color='k', linestyle='--', lw=0.5)
ax[0,0].set_xlabel(r'time [3-days]')
ax[0,0].set_ylabel(r'stress $\sigma$ [kPa]')
ax[0,0].legend(loc='lower center')
ax[0,0].set_xlim(-2.0,202.0)
ax[0,0].set_ylim(0.0,300.0)

# time-displacement plot
ax[1,2].plot(data[:, 0], data[:, 1], lw=1, label=r'$u_{outer}$')
ax[1,2].plot(data[:, 0], data[:, 2], lw=1, label=r'$u_{inner}$')
ax[1,2].axvline(x=data[20, 0], color='k', linestyle='--', lw=0.5)
#ax[1,0].axvline(x=data[72, 0], color='k', linestyle='--', lw=0.5)
#ax[1,0].axvline(x=data[92, 0], color='k', linestyle='--', lw=0.5)
ax[1,2].set_xlabel(r'time [3-days]')
ax[1,2].set_ylabel(r'displacement $u$ [mm]')
ax[1,2].legend()
ax[1,2].set_xlim(-2.0,202.0)
ax[1,2].set_ylim(0.60,0.75)

# time-stretch plot
ax[0,1].plot(data[:, 0], 1.0/np.sqrt(data[:, 8]), lw=1, label=r'$\lambda^{mat}_{r11}$')
ax[0,1].plot(data[:, 0], 1.0/np.sqrt(data[:, 9]), lw=1, label=r'$\lambda^{mat}_{r22}$')
ax[0,1].plot(data[:, 0], 1.0/np.sqrt(data[:, 10]), lw=1, label=r'$\lambda^{mat}_{r33}$')
ax[0,1].plot(data[:, 0], data[:, 11], lw=1, label=r'$\lambda^{fib}_r$')
ax[0,1].plot(data[:, 0], data[:, 12], lw=1, label=r'$\lambda^{cell}_r$')
ax[0,1].axvline(x=data[20, 0], color='k', linestyle='--', lw=0.5)
#ax[0,1].axvline(x=data[72, 0], color='k', linestyle='--', lw=0.5)
#ax[0,1].axvline(x=data[92, 0], color='k', linestyle='--', lw=0.5)
ax[0,1].set_xlabel(r'time [3-days]')
ax[0,1].set_ylabel(r'remodeling stretch $\lambda$ [-]')
ax[0,1].legend()
ax[0,1].set_xlim(-2.0,202.0)
ax[0,1].set_ylim(0.5,1.9)

# mass-fractions
ax[1,1].plot(data[:, 0], data[:, 15], lw=1, label=r'$\varrho_R/\varrho_{R0}$')
ax[1,1].plot(data[:, 0], data[:, 17], lw=1, label=r'$\phi^{mat}$')
ax[1,1].plot(data[:, 0], data[:, 18], lw=1, label=r'$\phi^{fib}$')
ax[1,1].plot(data[:, 0], data[:, 20], lw=1, label=r'$\varrho_R^{cell}/\varrho_{R0}^{cell}$')
ax[1,1].axvline(x=data[20, 0], color='k', linestyle='--', lw=0.5)
#ax[1,1].axvline(x=data[72, 0], color='k', linestyle='--', lw=0.5)
#ax[1,1].axvline(x=data[92, 0], color='k', linestyle='--', lw=0.5)
ax[1,1].set_xlabel(r'time [3-days]')
ax[1,1].set_ylabel(r'mass fraction $\phi$ [-]')
ax[1,1].legend()
ax[1,1].set_xlim(-2.0,202.0)
ax[1,1].set_ylim(0.35,1.25)

# cell-stress
ax[0,2].plot(data[:, 0], data[:, 7]*1.0e3, lw=1, label=r'$\sigma^{cell}$')
ax[0,2].plot(data[:, 0], data[:, 23]*1.0e3, lw=1, label=r'$\sigma_{pas}^{cell}$')
ax[0,2].plot(data[:, 0], data[:, 24]*1.0e3, lw=1, label=r'$\sigma_{act}^{cell}$')
ax[0,2].axhline(y=data[20, 7]*1.0e3, color='k', linestyle='--', lw=0.5)
ax[0,2].axvline(x=data[20, 0], color='k', linestyle='--', lw=0.5)
#ax[0,2].axvline(x=data[72, 0], color='k', linestyle='--', lw=0.5)
#ax[0,2].axvline(x=data[92, 0], color='k', linestyle='--', lw=0.5)
ax[0,2].set_xlabel(r'time [3-days]')
ax[0,2].set_ylabel(r'stress $\sigma$ [kPa]')
ax[0,2].legend() #loc='lower right'
ax[0,2].set_xlim(-2.0,202.0)
ax[0,2].set_ylim(-5.0,100.0)

# jacobians
ax[1,0].plot(data[:, 0], data[:, 13], lw=1, label=r'$J^{gnd}$')
ax[1,0].plot(data[:, 0], data[:, 14], lw=1, label=r'$J_e^{gnd}$')
ax[1,0].plot(data[:, 0], data[:, 15], lw=1, label=r'$J_g^{gnd}$')
ax[1,0].plot(data[:, 0], data[:, 20], lw=1, label=r'$J_g^{cell}$')
#ax[1,2].plot(data[:, 0], data[:, 16], lw=1, label=r'$J_g^*$')
ax[1,0].axvline(x=data[20, 0], color='k', linestyle='--', lw=0.5)
#ax[1,2].axvline(x=data[72, 0], color='k', linestyle='--', lw=0.5)
#ax[1,2].axvline(x=data[92, 0], color='k', linestyle='--', lw=0.5)
ax[1,0].set_xlabel(r'time [3-days]')
ax[1,0].set_ylabel(r'jacobian [-]')
ax[1,0].legend()
ax[1,0].set_xlim(-2.0,202.0)
ax[1,0].set_ylim(0.99,1.25)

# place a text box in bottom right
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr_u_outer = ( fr'error[$u_o$] = {u_outer_error:.2f}-%' )
textstr_u_inner = ( fr'error[$u_i$] = {u_inner_error:.2f}-%' )
# place a text box in upper left in axes coords, cylinder0={y=0.83,y=0.16}; cylinder1={y=0.60,y=0.08}; cylinder2={y=0.87,y=0.12}; cylinder3={y=0.92,y=0.34}
plt.text(0.10, 0.42, textstr_u_outer, transform=ax[1,2].transAxes, fontsize=10,
        verticalalignment='top', bbox=props) #
plt.text(0.10, 0.14, textstr_u_inner, transform=ax[1,2].transAxes, fontsize=10,
        verticalalignment='top', bbox=props) #

#fig.tight_layout()
plt.show
FIGURENAME = name + '_graph.pdf'
plt.savefig(FIGURENAME)
plt.close('all')
