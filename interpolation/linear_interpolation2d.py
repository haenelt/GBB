def linear_interpolation2d(x, y, x_s, y_s, f_s):
    # bilinear interpolation
    import numpy as np
    
    #x_s = np.array([12, 24])
    #y_s = np.array([3, 10])
    #f_s = np.array([[2, 5],[10, 20]])
    #x = 15
    #y = 10
    
    Delta_x = x_s[1] - x_s[0]
    Delta_y = y_s[1] - y_s[0]
    
    A = 1 / (Delta_x*Delta_y)
    B = np.array([x_s[1] - x, x - x_s[0]])
    C = [y_s[1] - y, y - y_s[0]]
    
    f = A * np.dot(B,np.dot(f_s,C))
    
    return f