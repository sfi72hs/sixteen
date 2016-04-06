import numpy as np
import networkx as nx
import scipy.stats

def init_graph(N, graph_density):
    mx_init = np.random.random( (N,N) )
    mx_init = mx_init < (graph_density)
    mx_init = np.tril(mx_init)
    mx_init = mx_init + mx_init.T
    np.fill_diagonal(mx_init, 0)
    return mx_init

def get_inf_assortativity(mx, is_infected):
    G=nx.from_numpy_matrix(mx)
    nx.set_node_attributes(G, "infected", dict(enumerate(is_infected)))
    return nx.attribute_assortativity_coefficient(G, 'infected')

def run_rewire(mx_init, is_infected_init, benefit_function, opts={}):
    """
    # Parameters:
        mx_init          - NxN numpy array adjacency matrix to start with
        is_infected_init - N bool vector, who is infected at beginning
        benefit_function - benefit function with signature f(self_state, neighs_infected, neighs_noninfected, global_infected)
        opts             - Options dictionary, see code below
    
    # Returns dictionary with entries:
        num_infected  - list of number of infected nodes over iterations 
        num_rewired   - list of number of rewirings (cumulative) over iterations
        assortativity - infected-or-not-attribute assortativity, over iterations (downsampled)
        is_infected   - final iteration vector of which nodes are infected
        mx            - final iteration adjacency matrix
        rewire_types  - types of rewiring tracked by rewire_counts_list
        rewite_counts_list = rewite_counts_list,
    )

    """
    beta       = opts.get('beta', 1)          # strength of beta (influence of benefit on rewiring choices)
    NUM_ITERS = opts.get('NUM_ITERS', 10000)  # number of iterations to run
    
    p_recovery = opts.get('p_recovery', 0.)   # probability of node recovering
    p_transmit = opts.get('p_transmit', 0.001)  # probability of transmission across link, per each node's iteration
    p_rewire   = opts.get('p_rewire'  , 0.5)  # probability of rewiring/adding, per node's iteration

    do_add     = opts.get('do_add'   , False)      # add links adaptive network change?
    do_rewire  = opts.get('do_rewire', True)       # rewire links during adaptive network change?
    do_null    = opts.get('do_null'  , False)  # allow null (benefit-neutral) rewiring.adding? 
    do_only_beneficial = opts.get('do_only_beneficial'  , True)  # only add/rewire if its benefit-enhancing
    
    iterate_only_infected  = opts.get('iterate_only_infected', False)       # do spreading / adaptive network change only for infected noedes
    
    Ksteps    = None # 3
    do_only_beneficial = True
    
    def expected_benefit(is_node_infected, num_neighs_infected, num_neighs_noninfected, num_global_infected):
        num_neighbors = num_neighs_infected + num_neighs_noninfected
        if not is_node_infected:
            prob_willbe_infected = 1.0 - (1.0 - p_transmit)**num_neighs_infected
        else:
            prob_willbe_infected = 1.0 - p_recovery
            
        prob_Nnewneighbors_willbe_infected = np.zeros(num_neighs_noninfected+1)
        
        # don't consider their recovery
        if not is_node_infected:
            prob_Nnewneighbors_willbe_infected[0]=1.0
        else:
            for i in range(num_neighs_noninfected+1):
                prob_Nnewneighbors_willbe_infected[i] = scipy.stats.binom.pmf(i, num_neighs_noninfected, p_transmit)

        exp_benefit = (1-prob_willbe_infected) * \
            sum([p*benefit_function(0, num_neighs_infected+newinfs, num_neighs_noninfected-newinfs, num_global_infected+newinfs-is_node_infected)
            for newinfs, p in enumerate(prob_Nnewneighbors_willbe_infected) ])
        
        exp_benefit += prob_willbe_infected * \
            sum([p*benefit_function(1, num_neighs_infected+newinfs, num_neighs_noninfected-newinfs, num_global_infected+newinfs+(1-is_node_infected))
            for newinfs, p in enumerate(prob_Nnewneighbors_willbe_infected) ])
        return exp_benefit

    N = len(is_infected_init)
    assert(N == mx_init.shape[0])

    is_infected = is_infected_init.copy()
    infection_time = np.zeros(N, dtype='int') - 1
    infection_time[is_infected] = 0

    debug     = False

    mx = mx_init.copy()
    mxPlusI = mx.copy()
    np.fill_diagonal(mxPlusI, True)
    mxPlusIint = mxPlusI.astype('int')

    meanK = mx.sum(axis=0).mean()
    #print "MEAN K: %0.5f" % meanK
    #print "Inf chance: %0.5f" % (1 - (1-p_transmit)**meanK)

    c_rewired = 0
    rewire_types = ['i-i->n','n-i->n','i-n->i','n-n->i',]
    if do_only_beneficial:
        rewire_types += ['i-none','n-none']
    if do_null:
        rewire_types += 'i-i->i','n-i->i','i-n->n','n-n->n'
    if do_add:
        rewire_types += ['i+i','i+n','n+i','n+n']
    rewire_counts = {}
    for rt in rewire_types:
        rewire_counts[rt] = 0
    rewire_counts_list = []
    num_infected = []
    num_rewired  = []
    mean_degree  = []
    mean_degree_inf = []
    assortativity = []

    for iteration in range(int(NUM_ITERS)):
        if iteration % 100 == 0:
            assortativity.append(get_inf_assortativity(mx, is_infected))
            
        rewire_counts_list.append([rewire_counts[rt] for rt in rewire_types])
        mean_degree.append(mx.sum(axis=0).mean())
        mean_degree_inf.append(mx[is_infected,:].sum(axis=1).mean())
        num_global_infected = is_infected.sum()
        num_infected.append(num_global_infected)
        num_rewired.append(c_rewired)
        
        # choose "ego" node to rewire
        if iterate_only_infected:
            node = np.random.choice(np.flatnonzero(is_infected)) 
        else:
            node = np.random.choice(N)

        neighbors = mx[node,:]
        num_neighbors = neighbors.sum()
        neighbor_status = is_infected[neighbors]
        num_neighs_infected = neighbor_status.sum()
        num_neighs_noninfected = (~neighbor_status).sum()
        
        if (np.random.random() < p_rewire):
            #radiusK = (np.linalg.matrix_power(mxPlusIint, Ksteps) > 0).astype('int')
            #possible_switch_mx = radiusK - mxPlusIint
            if Ksteps is None:
                possible_new_neighs = np.ones(N, dtype='bool')
            else:
                possible_new_neighs = np.zeros(N, dtype='bool')
                possible_new_neighs[node] = 1
                for i in range(Ksteps):
                    for j in np.flatnonzero(possible_new_neighs):
                        possible_new_neighs[mx[j,:]] = 1
            possible_new_neighs = possible_new_neighs & ~neighbors
            #continue

            #possible_new_neighs = possible_switch_mx[node,:]
            num_radiusK = possible_new_neighs.sum()
            num_radiusK_infected = is_infected[possible_new_neighs].sum()
            num_radiusK_noninfected = num_radiusK - num_radiusK_infected
            #prop_radiusK_infected = float(num_radiusK_infected)/num_radiusK

            #if debug and np.logical_and(possible_switch_mx != 0, possible_switch_mx != 1 ).any():
            #    raise Exception('error')
            if debug and num_radiusK_noninfected < 0 or num_radiusK_infected < 0:
                raise Exception('error -- radius both noninfected and infected less than 0')


            curB = expected_benefit(is_infected[node], num_neighs_infected, num_neighs_noninfected, num_global_infected)

            newBinf2noninf, newBnoninf2inf, newBnewinfedge, newBnewnoninfedge = -np.inf, -np.inf, -np.inf, -np.inf
            # given we swap infected neighbor for non-infected neighbor
            if num_neighs_infected > 0 and num_radiusK_noninfected > 0:
                newBinf2noninf = expected_benefit(is_infected[node], num_neighs_infected-1, num_neighs_noninfected+1, num_global_infected)
            # given we swap non-infected neighbor for infected neighbor
            if num_neighs_noninfected > 0 and num_radiusK_infected > 0:
                newBnoninf2inf = expected_benefit(is_infected[node], num_neighs_infected+1, num_neighs_noninfected-1, num_global_infected)
            if num_radiusK_infected > 0:
                newBnewinfedge = expected_benefit(is_infected[node], num_neighs_infected+1, num_neighs_noninfected, num_global_infected)
            if num_radiusK_noninfected > 0:
                newBnewnoninfedge = expected_benefit(is_infected[node], num_neighs_infected, num_neighs_noninfected+1, num_global_infected)
                
            if debug:
                print '#n-inf: %d, #n-ninf: %d, #r-inf: %d, #r-ninf: %d' % (
                    num_neighs_infected, num_neighs_noninfected, num_radiusK_infected, num_radiusK_noninfected)

            if num_radiusK and (do_add or (do_rewire and num_neighbors>0)):
                #print curB , newBinf2noninf, newBnoninf2inf
                #log_ps = np.array([curB , newBinf2noninf, newBnoninf2inf])
                
                log_ps = []
                if do_null:
                    if do_only_beneficial:
                        raise Exception('do not use do_only_beneficial and do_null both')
                    log_ps.append(0)
                else:
                    log_ps.append(-np.inf)
                    
                if do_rewire and num_neighbors>0:
                    log_ps += [newBinf2noninf-curB, newBnoninf2inf-curB]
                else:
                    log_ps += [-np.inf, -np.inf]
                    
                if do_add:
                    log_ps += [newBnewinfedge-curB, newBnewnoninfedge-curB]
                else:
                    log_ps += [-np.inf, -np.inf]

                log_ps = np.array(log_ps)
                if beta == 0:
                    log_ps[:] = 0
                else:
                    log_ps *= beta

                if do_only_beneficial:
                    log_ps[log_ps<=0] = -np.inf
                    
                ctype = 'i' if is_infected[node] else 'n'
                if np.all(np.isinf(-log_ps)):
                    rewire_counts['%s-none' % ctype] += 1
                else:
                    #log_ps = np.array([0 , newBinf2noninf-curB, newBnoninf2inf-curB, newBnewinfedge-curB, newBnewnoninfedge-curB])
                    #log_ps = np.array([0 , -10e4, -10e4, newBnewinfedge-curB, newBnewnoninfedge-curB])
                    prob_rewire = np.exp(log_ps)
                    prob_rewire /= prob_rewire.sum()

                    rewire_action = np.random.choice(len(prob_rewire), p=prob_rewire)

                    assert(len(prob_rewire)==5)
                    if rewire_action == 0:
                        prop_neigh_infected = (float(num_neighs_infected)/num_neighbors) if num_neighbors > 0 else 1
                        prop_radiusK_infected = (float(num_radiusK_infected)/num_radiusK) if num_radiusK > 0 else 1 # TODO 
                        #p_switch_inf2noninf = prop_neigh_infected*(1-prop_radiusK_infected)
                        #p_switch_noninf2inf = (1-prop_neigh_infected)*prop_radiusK_infected
                        #p_switch_same = 1 - p_switch_inf2noninf - p_switch_noninf2inf

                        prop_both_infected = (prop_neigh_infected*prop_radiusK_infected)
                        prop_both_infected = prop_both_infected/(prop_both_infected+(1-prop_neigh_infected)*(1-prop_radiusK_infected))
                        if np.random.random() < prop_both_infected:
                            node_to_rewire_to = np.random.choice(np.flatnonzero(is_infected & possible_new_neighs))
                            node_to_unwire_from = np.random.choice(np.flatnonzero(is_infected & neighbors))
                            rewire_counts['%s-i->i' % ctype] += 1
                        else:
                            node_to_rewire_to = np.random.choice(np.flatnonzero(~is_infected & possible_new_neighs))
                            node_to_unwire_from = np.random.choice(np.flatnonzero(~is_infected & neighbors))
                            rewire_counts['%s-n->n' % ctype] += 1
                    if rewire_action == 1:
                        node_to_unwire_from = np.random.choice(np.flatnonzero(is_infected & neighbors))
                        node_to_rewire_to = np.random.choice(np.flatnonzero(~is_infected & possible_new_neighs))
                        rewire_counts['%s-i->n' % ctype] += 1
                    if rewire_action == 2:
                        node_to_unwire_from = np.random.choice(np.flatnonzero(~is_infected & neighbors))
                        node_to_rewire_to = np.random.choice(np.flatnonzero(is_infected & possible_new_neighs))
                        rewire_counts['%s-n->i' % ctype] += 1
                    if rewire_action == 3:
                        node_to_unwire_from = None
                        node_to_rewire_to = np.random.choice(np.flatnonzero(is_infected & possible_new_neighs))
                        rewire_counts['%s+i' % ctype] += 1
                    if rewire_action == 4:
                        node_to_unwire_from = None
                        node_to_rewire_to = np.random.choice(np.flatnonzero(~is_infected & possible_new_neighs))
                        rewire_counts['%s+n' % ctype] += 1

                    if debug and iteration % 20 == 0:
                        print is_infected[node], prob_rewire, num_neighs_noninfected, num_neighs_infected, num_radiusK_noninfected, num_radiusK_infected
                        print newBnewinfedge, newBnewnoninfedge
                
                    mx[node, node_to_rewire_to] = 1
                    mx[node_to_rewire_to, node] = 1
                    if node_to_unwire_from is not None:
                        mx[node, node_to_unwire_from] = 0
                        mx[node_to_unwire_from, node] = 0

                    c_rewired+=1
                    if debug:
                        print "Rewiring %d to %d" % (node_to_unwire_from, node_to_rewire_to)

                    mxPlusI = mx.copy()
                    np.fill_diagonal(mxPlusI, True)
                    mxPlusIint = mxPlusI.astype('int')

                    neighbors = mx[node,:]
                    num_neighbors = neighbors.sum()
                    neighbor_status = is_infected[neighbors]
                    num_neighs_infected = neighbor_status.sum()
                    num_neighs_noninfected = (~neighbor_status).sum()     
            
            
        total_links = mx.sum()

        SIlinks = mx[is_infected,:][:,~is_infected]
        edge_targets = np.nonzero(SIlinks)[1]
        if len(edge_targets):
            edges_to_infect = np.random.random(len(edge_targets)) < p_transmit
            #cedge = edges[np.random.choice(len(edges))]
            not_infected_ixs = np.flatnonzero(~is_infected)
            #print "Setting ",is_infected[not_infected_ixs[cedge[1]]]
            #print is_infected[not_infected_ixs[edge_targets[edges_to_infect]]]
            is_infected[not_infected_ixs[edge_targets[edges_to_infect]]] = True

        """
        if is_infected[node] and num_neighs_noninfected > 0:
            #print "HER", num_neighs_noninfected
            #oldB = benefit_function(is_infected[i], num_neighs_infected, num_neighs_noninfected, num_global_infected)
            #newB = benefit_function(0, num_neighs_infected + 1, num_neighs_noninfected - 1, num_global_infected + 1)

            #prob_transmit = np.exp(beta * (newB - oldB - cost))
            #prob_transmit = prob_transmit / (1. + prob_transmit)
            
            prob_transmit = p_transmit # * (num_neighbors / float(total_links))

            neighs_to_infect = np.random.random(num_neighs_noninfected) < prob_transmit
            if iteration % 100 == 0: # debug:
                #print iteration, "Infected %d -> %s" % (node, neighs_noninfected_ixs[neighs_to_infect])
                print iteration,  "degree=", num_neighbors, "num_infected=", neighs_to_infect.sum()
                print num_neighs_noninfected, neighs_to_infect
                print 'mean inf degree:', mx[is_infected,:].sum(axis=0).mean()

            if neighs_to_infect.any():
                neighbor_ixs = np.flatnonzero(neighbors)
                neighs_noninfected_ixs = neighbor_ixs[~neighbor_status]
                is_infected[neighs_noninfected_ixs[neighs_to_infect]] = True
        """

        if np.random.random() < p_recovery:
            is_infected[i] = False

    return dict(
        num_infected = num_infected, 
        num_rewired  = num_rewired,
        mean_degree = mean_degree,
        mean_degree_inf = mean_degree_inf,
        assortativity = assortativity,
        is_infected = is_infected,
        mx = mx,
        rewire_types=rewire_types,
        rewire_counts_list = rewire_counts_list,
    )


            