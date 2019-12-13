def cost_BBR(A=None):
    # cost function
        
    # Qv: percent contrast measure
    # Q0: offset parameter
    # Mv: slope parameter
    # hv: weight for each vertex
    # gv: gray matter intensity at distance d away from the surface mesh along a chosen axis
    # wv: white matter intensity at distance d away from the surface mesh along a chosen axis
    # J: cost function
    # B. subset of vertices

    Q0 = 0
    Mv = 0.5
    hv = 1

    #Qv = 100 * ( gv -wv ) / ( 0.5 * ( gv + wv ) )

    # J = 1 / N sum_{v in B} hv * ( 1 + tanh(Mv * (Qv - Q0)) )
