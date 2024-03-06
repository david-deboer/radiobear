# Benchmark radiobear against Cheng Li radio code (inherited from M. Janssen) 
import radiobear as rb 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
cmap = cm.magma 


# Initialize Jupiter 
config_file = 'config_benchmark.par' 
j = rb.planet.Planet('jupiter', plot_atm=False, plot_bright=False, verbose=True, config_file = config_file) 

# Run radiobear at a set of emission angles at the equator 
freqs = [0.6, 1.25, 2.6, 5.2, 10, 21.9] 
b = [[0.0,0.0],[np.sin(np.radians(15)),0.0],[np.sin(np.radians(30)),0.0],[np.sin(np.radians(45)),0.0]]
rv = j.run(freqs,b=b, reuse_override='False')


# Expexted results 
T_rb = np.array(
	   [[805.38385, 424.58337, 298.74915, 227.0338 , 179.9148 , 147.3627 ],
       [790.66736, 420.24728, 296.62033, 225.69223, 179.22531, 147.24104],
       [745.55975, 407.28714, 290.07492, 221.56863, 177.13437, 146.8625 ],
       [668.4528 , 385.59015, 278.51672, 214.30847, 173.56506, 146.17966]]) 


# Comparison values from Cheng Li's radio code 
T_CL = np.array(
	   [[823.63, 424.45, 299.05, 228.63, 181.26, 142.84],
       [807.45, 420.09, 296.96, 227.3 , 180.55, 142.67],
       [758.75, 407.09, 290.56, 223.24, 178.39, 142.16],
       [677.81, 385.46, 279.27, 216.06, 174.69, 141.22]])

 

# Plot the results 
fig, ax = plt.subplots(2, 1,figsize=(10,8))
plt.subplots_adjust(hspace=0)
labels = [] 
for i in range(len(b)): 
	ax[0].plot(freqs,rv.Tb[i,:],linestyle='-',  color=cmap(i/(len(b)-1)))
	ax[0].plot(freqs,T_rb[i,:], linestyle='--', color=cmap(i/(len(b)-1)))
	ax[0].plot(freqs,T_CL[i,:], linestyle=':',  color=cmap(i/(len(b)-1)))
	ax[1].plot(freqs,rv.Tb[i,:]-T_rb[i,:],linestyle='--',  color=cmap(i/(len(b)-1)))
	ax[1].plot(freqs,rv.Tb[i,:]-T_CL[i,:], linestyle=':', color=cmap(i/(len(b)-1)))
	labels.append(f'{np.degrees(np.arcsin(b[i][0])):0.0f}deg')
ax[0].plot([],[],label='This installation',color='gray',linestyle='-')
ax[0].plot([],[],label='radiobear benchmark',color='gray',linestyle='--')
ax[0].plot([],[],label='CL benchmark',color='gray',linestyle=':')
legend = ax[0].legend() 
ax1 = ax[0].add_artist(legend) 

line1, = plt.gca().plot([],[],color=cmap(0/(len(b)-1)),linestyle='-')
line2, = plt.gca().plot([],[],color=cmap(1/(len(b)-1)),linestyle='-')
line3, = plt.gca().plot([],[],color=cmap(2/(len(b)-1)),linestyle='-')
line4, = plt.gca().plot([],[],color=cmap(3/(len(b)-1)),linestyle='-')



legend2 = plt.legend([line1,line2,line3,line4],labels)

ax[0].set_title('Spectrum comparison')
ax[0].set_ylabel('Brightness temperature (K)')
ax[0].set_xticks([])

ax[1].set_ylabel(r'$\Delta$ Brightness temperature (K)')
ax[1].set_xlabel('Frequency (GHz)')


# Plot the profile 
# Convert from ppm to solar abundances (based on Asplund Solar abundances)
N2H = 1.48e-4 # N/H2 ratio 
O2H = 1.07e-3 # O/H2 ratio 
S2H = 2.89e-5 # S/H2 ratio 

fig, ax = plt.subplots(1, 1,figsize=(10,8))
i = 0 
ax.plot( j.atmos[0].gas[j.atmos[0].config.C['NH3']]/N2H, j.atmos[0].gas[j.atmos[0].config.C['P']], color='firebrick',label='NH3',linestyle= '-')
ax.plot( j.atmos[0].gas[j.atmos[0].config.C['H2O']]/O2H, j.atmos[0].gas[j.atmos[0].config.C['P']], color='navy',label='H2O',linestyle= '-')
ax.plot( j.atmos[0].gas[j.atmos[0].config.C['H2S']]/O2H, j.atmos[0].gas[j.atmos[0].config.C['P']], color='lightsalmon',label='H2S',linestyle= '-')
ax2 = ax.twiny() 
ax2.plot(j.atmos[0].gas[j.atmos[0].config.C['T']], j.atmos[0].gas[j.atmos[0].config.C['P']],color='gray',linestyle='--',label='Temperature')
ax2.set_xlabel('Tempeature (K)')
ax2.xaxis.label.set_color('gray')
ax2.tick_params(axis='x',colors = 'gray') 
ax2.spines['top'].set_color('gray')
ax.invert_yaxis()
ax.set_yscale('log')
ax.set_ylim([500,0.1])

ax.set_ylabel('P (bar)')
ax.set_xlabel(r'Abundance (solar)')
ax.legend(loc='lower left')
plt.show()





