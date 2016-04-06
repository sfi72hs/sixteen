import networkx as nx
import numpy as np

from rewire import *


N = 500  # number of nodes 
graph_density = 0.1
mx_init = init_graph(N, graph_density)
is_infected_init = (np.random.random(N) < .05)

def bf(self_state, neighs_infected, neighs_noninfected, global_infected):
    return self_state + neighs_infected

rdict = run_rewire(mx_init, is_infected_init, benefit_function=bf, opts=dict(NUM_ITERS=1000))

print 'Last 10 iterations infected (out of 1000 iterations):',rdict['num_infected'][-10:]
import matplotlib.pyplot as plt
plt.plot(rdict['num_infected'])