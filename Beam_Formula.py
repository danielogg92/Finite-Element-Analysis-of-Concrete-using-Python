# -*- coding: utf-8 -*-

def bar_formula(P1,U1,E,A,x_ratio,length,wx):
    """
    calculates bar force and displacement.

    Parameters
    ----------
    P0 : float
        Axial force at End 1. Negative is compression (in N)
    U0 : float
        Displacement at end 1 (in mm)
    E : float
        Youngs modulus in MPa.    
    A : float
        Cross-sectional area in sq.mm.
    x_ratio : float
        ratio of position relative to end 1. Can vary from 0 to 1.
    length : float
        length in mm.
    wx : float
        traction. Postive goes from end 1 to end 2. In N/mm
    Returns
    -------
    Px : float
        axial force relative to end 1
    Ux : float
        displacement mm relative to end 1
    """
    #Get the x-position
    x=x_ratio*length;
    #Calculate the axial force
    Px=P1-wx*x;
    #Calculate the axial displacement
    Ux=U1+1/(A*E)*(P1*x-1/2*wx*x**2);
    return Px,Ux;

def beam_formula(V1,M1,U1,R1,E,Ix,x_ratio,length,wy):
    """
    Calculates the flexural deformations and forces, rotations and
    displacements

    Parameters
    ----------
    V1 : float
        shear force at end 1 (in N).
    M1 : float
        moment at end 1 (in N-mm).
    U1 : float
        displacement at end 1 (in mm).
    R1 : float
        rotation at end 1 (in rads (mm/mm)).
    E : float
        youngs modulus (in MPa or N/mm**2).
    Ix : float
        second moment of area (in mm**4).
    x_ratio : float
        ratio along the element where the section cut is taken. Value 
        varies from 0 to 1.
    length : float
        length of the element (in mm).
    wy : float
        applied loading (in N/mm).

    Returns
    -------
    Vx : float
        shear force at section cut (in N).
    Mx : float
        moment at section cut (in N-mm).
    Rx : float
        rotation at section cut (in rad or mm/mm).
    Ux : float
        displacement at section cut (in mm).
    """
    #Get the x-position
    x=x_ratio*length;
    #Get the shear at section x
    Vx=V1+wy*x;
    #Get the moment at section x
    Mx=M1+V1*x+1/2*wy*x**2;
    #Get the rotation at section x
    Rx=R1+1/(E*Ix)*(M1*x+1/2*V1*x**2+1/6*wy*x**3);
    #Get the displacement at section x
    Ux=U1+R1*x+1/(E*Ix)*(1/2*M1*x**2+1/6*V1*x**3+1/24*wy*x**4);
    return Vx,Mx,Rx,Ux;