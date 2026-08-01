# python3 test_plot.py
#
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------------------#
data_input = np.loadtxt("data_test_2d.dat", dtype=float)

#-------------------------------------------------------------------#
plt.rcParams.update({'font.size':10})
fig = plt.figure(figsize=plt.figaspect(0.5))
gs = fig.add_gridspec(1, 2, hspace=0.25, wspace=0.25)
ax = gs.subplots(sharex=False, sharey=False)

#ax.axhline(y=0.0, color='k', linestyle='--', linewidth=1)
ax[0].plot(data_input[:, 3] + 1.0, data_input[:, 1], lw=1, label='stress-strain')
ax[1].plot(data_input[:, 3] + 1.0, 100.0*(data_input[:, 5]/data_input[0, 5]-1.0), lw=1, label='volume')

ax[0].set_xlabel(r'stretch $\lambda$ [-]')
ax[0].set_ylabel(r'stress $\sigma$ [MPa]')
ax[0].legend()
ax[1].set_xlabel(r'stretch $\lambda$ [-]')
ax[1].set_ylabel(r'volume change $V/V_R$ [%]')
ax[1].legend()

#fig.tight_layout()
plt.show
FIGURENAME = 'test_2d_graph.pdf'
plt.savefig(FIGURENAME)
plt.close('all')
