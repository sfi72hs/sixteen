import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.use('agg')
mpl.rcParams['text.usetex']=True

mpl.rcParams.update(mpl.rcParamsDefault)

# Sets a few parameters regarding the plot.
plt.rcParams['text.usetex']       = True
plt.rcParams['font.size']         = 34*.5
plt.rcParams['font.family']       = 'serif'
plt.rcParams['font.serif']        = 'Computer Modern Roman'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.major.size']  = 8
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.major.size']  = 8
