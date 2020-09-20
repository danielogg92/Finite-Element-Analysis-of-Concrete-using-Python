# -*- coding: utf-8 -*-
def rect_prop(x_dim,y_dim):
    """
    @params x_dim is width
    @params y_dim is depth
    
    @returns A,Ixx,Iyy,j
    """
    A=x_dim*y_dim;
    Ixx=1/12*x_dim*y_dim**3;
    return A,Ixx;