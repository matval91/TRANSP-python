# MODULE TO CREATE HIGH-QUALITY PLOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

#=====================================================================================
# SET TEXT FONT AND SIZE
#=====================================================================================
plt.rc('font', family='serif', serif='Palatino')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
#=====================================================================================


#x=np.linspace(0,6.28,100)
#y=np.cos(x)

fig=plt.figure()
ax=fig.add_subplot(111)


#=====================================================================================
# SET LINE GRAPHIC
#=====================================================================================
ax.plot(rho,jnbi,'r', linewidth=3.3, label='NEMO/SPOT')
ax.plot(rho_ascot, j_ascot, linewidth=3.3, color='b', label='BBNBI/ASCOT')
ax.legend()
#=====================================================================================

#=====================================================================================
# SET AXIS LABELS
#=====================================================================================
ax.set_xlabel(r'x [mm] $\rho$')
ax.set_ylabel(r'\TeX\ style! $\sum e^{-i\pi}$')
#=====================================================================================


#=====================================================================================
# ADJUST SUBPLOT IN FRAME
#=====================================================================================
plt.subplots_adjust(top=0.95,bottom=0.2,left=0.2,right=0.95)
#=====================================================================================

#=====================================================================================
# SET TICK LOCATION
#=====================================================================================

# Create your ticker object with M ticks
M = 4
yticks = ticker.MaxNLocator(M)
yticks_m=ticker.MaxNLocator(M*2)
xticks = ticker.MaxNLocator(M)
# Set the yaxis major locator using your ticker object. You can also choose the minor
# tick positions with set_minor_locator.
ax.yaxis.set_major_locator(yticks)
#ax.yaxis.set_minor_locator(yticks_m)
ax.xaxis.set_major_locator(xticks)
#=====================================================================================



#=====================================================================================
# ADD ANOTHER SUBGRAPH (with subaxes) in the plot
#=====================================================================================
a = plt.axes([0.3, 0.6, .2, .2])#, facecolor='y')
plt.plot(x, np.sin(x))
plt.title('sin(x)')
plt.xlim(0, 6.28)
plt.xticks([])
plt.yticks([])
#=====================================================================================

#=====================================================================================
# Change ticks parameters
#=====================================================================================
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
#=====================================================================================






# SHOW
plt.show()
