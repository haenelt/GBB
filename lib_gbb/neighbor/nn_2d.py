def nn_2d(ind, adjm, n_iter):
    """
    This function finds all nearest neigbors on the surface mesh which show a connection within the
    numbers of iterations steps.
    Inputs:
        *ind: vertex index.
        *adjm: sparse adjacency matrix (in csr_matrix format).
        *n_iter: number of iterations.
    Outputs:
        *nn: array of neighbor indices.
        
    created by Daniel Haenelt
    Date created: 21-12-2019     
    Last modified: 22-12-2019
    """
    import numpy as np

    # get first order neighbours
    nn = adjm[ind,:].indices

    i = 0
    while i < n_iter - 1:
        nn_temp = nn.copy()
        for j in range(len(nn_temp)):
            nn = np.append(nn, adjm[nn_temp[j],:].indices)
            nn = np.unique(nn)

        i += 1
        
    return nn
