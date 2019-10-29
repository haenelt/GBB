def linear_interpolation2d(x, y, x_s, y_s, f_s):
    # bilinear interpolation
    import numpy as np
    
    x_s = np.array(x_s)
    y_s = np.array(y_s)
    f_s = np.array(f_s)
    
    Delta_x = x_s[1] - x_s[0]
    Delta_y = y_s[1] - y_s[0]
    
    A = 1 / (Delta_x*Delta_y)
    B = np.array([x_s[1] - x, x - x_s[0]])
    C = [y_s[1] - y, y - y_s[0]]
    
    f = A * np.dot(B,np.dot(f_s,C))
    
    return f