# -*- coding: utf-8 -*-
import Beam_Formula as bf;
import math as m;
import matplotlib.pyplot as plt; 
import numpy as np;
import sys;
import time;

def beam_element_stiff(A,Ix,E,L):
    """
    Parameters
    ----------
    A : float
        Cross-sectional area.
    Ix : float
        Second moment of inertia.
    E : float
        Youngs modulus.
    L : float
        Length.

    Returns
    -------
    ele_mat : numpy.ndarray
        Elemental stiffness matrix.
    """
    ele_mat=np.array([[A*E/L,0,0,-A*E/L,0,0],
                      [0,12*E*Ix/L**3,6*E*Ix/L**2,0,-12*E*Ix/L**3,6*E*Ix/L**2],
                      [0,6*E*Ix/L**2,4*E*Ix/L,0,-6*E*Ix/L**2,2*E*Ix/L],
                      [-A*E/L,0,0,A*E/L,0,0],
                      [0,-12*E*Ix/L**3,-6*E*Ix/L**2,0,12*E*Ix/L**3,-6*E*Ix/L**2],
                      [0,6*E*Ix/L**2,2*E*Ix/L,0,-6*E*Ix/L**2,4*E*Ix/L]]);
    return ele_mat;

def transformation_6(cos_x,cos_y):
    """
    Beam transformation matrix. Elements has 6 DOF (Ux,Uy,Rz).

    Parameters
    ----------
    cos_x : float
        x-cosine of the frame element. Element orientated in the positive x-
        direction cos_x=1.0.
    cos_y : float
        y-cosine of the frame element. Element orientated in the positive y-
        direction cos_y=1.0.

    Returns
    -------
    transform : numpy array
        Numpy array of size (6,6).

    """
    transform=np.array([[cos_x,cos_y,0,0,0,0],
                        [-cos_y,cos_x,0,0,0,0],
                        [0,0,1,0,0,0],
                        [0,0,0,cos_x,cos_y,0],
                        [0,0,0,-cos_y,cos_x,0],
                        [0,0,0,0,0,1]]);
    return transform;

def rect_prop(x_dim,y_dim):
    """
    @params x_dim is width
    @params y_dim is depth
    
    @returns A,Ixx,Iyy,j
    """
    A=x_dim*y_dim;
    Ixx=1/12*x_dim*y_dim**3;
    return A,Ixx;

def frame_eq_node_loads(udl_x,udl_y,ele_len):
    """
    
    """
    shape_mat=np.array([[1/2,0],
                        [0,1/2],
                        [0,ele_len/12],
                        [1/2,0],
                        [0,1/2],
                        [0,-ele_len/12]]);
    load_mat=np.array([[udl_x],
                       [udl_y]]);
    frame_node_loads=shape_mat @ load_mat * ele_len;
    return frame_node_loads;

def print_nodal_data(node_no,nodes,node_loads,p_global,u_global):
    """
    Prints the nodal data post analysis
    
    Parameters
    ----------
    node_no : int
        Node number.
    nodes : dict
        The nodal dictionary containing coordinates, fixity and DOF numbers.
    node_loads : dict
        The node loads dictionary containined point loads.
    p_global : Numpy array
        The global loads matrix post-analysis.
    u_global : Numpy array
        The global displacements matrix post-analysis.

    Returns
    -------
    None.
    """
    #Get the nodal coordinates
    x,y=round(nodes[node_no]['x'],1),round(nodes[node_no]['y'],1);
    #Get the nodal restraint and DOF Id
    ux_fix,ux_dof=nodes[node_no]['fix_x'],nodes[node_no]['Ux_dof'];
    uy_fix,uy_dof=nodes[node_no]['fix_y'],nodes[node_no]['Uy_dof'];
    rz_fix,rz_dof=nodes[node_no]['fix_mz'],nodes[node_no]['Rz_dof'];
    print('NODE: {} \nX: {} mm    Y: {} mm'.format(node_no,x,y));
    print('X-restraint: {}    Y-restraint: {}    RZ-restraint: {}\n'.format(
        ux_fix,uy_fix,rz_fix));
    #Get the node loads defined by user
    if(node_loads.get(node_no,False)!=False):
        x_load=round(node_loads[node_no][0]*10**-3,3);
        y_load=round(node_loads[node_no][1]*10**-3,3);
        mz_moment_applied=round(node_loads[node_no][2]*10**-6,1);
    else:
        x_load,y_load,mz_moment=0.0,0.0,0.0;
    print('Applied Loads:');
    print('X-force: {} kN    Y-force: {} kN    MZ-moment: {} kN-m\n'.format(
        x_load,y_load,mz_moment_applied));
    #Get the nodal forces
    x_force=round(p_global[ux_dof-1][0]*10**-3,3);
    y_force=round(p_global[uy_dof-1][0]*10**-3,3);
    mz_moment=round(p_global[rz_dof-1][0]*10**-6,1);
    #Get the node reactions if restraint is True for any DOF
    if(ux_fix or uy_fix or rz_fix):
        reaction_string='NODE REACTIONS:\n';
        if(ux_fix):
            x_reaction=round(x_force-x_load,3);
            reaction_string+='RX = {} kN    '.format(x_reaction);
        if(uy_fix):
            y_reaction=round(y_force-y_load,3);
            reaction_string+='RY = {} kN    '.format(y_reaction);
        if(rz_fix):
            mz_reaction=round(mz_moment-mz_moment_applied,1);
            reaction_string+='RMZ = {} kN-m    '.format(mz_reaction);
        print(reaction_string+'\n');
    #Get the nodal displacements
    ux=round(u_global[ux_dof-1][0],3);
    uy=round(u_global[uy_dof-1][0],3);
    rz=round(u_global[rz_dof-1][0],5);
    print('NODAL DISPLACEMENTS:');
    print(('X-disp: {} mm    Y-displ: {} mm    Z-rot: {} rads\n'.format(
        ux,uy,rz)));
    return None;

def print_element_data(ele_no,elements,nodes,p_global,u_global):
    """
    Prints the element data post analysis

    Parameters
    ----------
    ele_no : int
        Element Number.
    elements : dict
        The element dictionary containing nodes and sections.
    p_global : numpy array
        The global loads matrix post-analysis.
    u_global : numpy array
        The global displacement matrix post-analysis

    Returns
    -------
    None.
    """
    #Get the nodes
    node_i=elements[ele_no]['Node_1'];
    node_j=elements[ele_no]['Node_2'];
    #Get the elements global DOF numbers
    element_dof=[nodes[node_i]['Ux_dof'],nodes[node_i]['Uy_dof'],
                 nodes[node_i]['Rz_dof'],nodes[node_j]['Ux_dof'],
                 nodes[node_j]['Uy_dof'],nodes[node_j]['Rz_dof']];
    #Get frame properties
    section=elements[ele_no]['Section'];
    A=sections[section]['A'];
    Ixx=sections[section]['Ixx'];
    mat=sections[section]['Material'];
    E=materials[mat]['E'];
    length=elements[ele_no]['Length'];
    #Get the element stiffness matrix
    k_local=elements[ele_no]['k_local'];
    #Get the local displacement from the u_global matrix
    u_local=np.zeros((6,1));
    for i in range(len(element_dof)):
        u_local[i][0]=u_global[element_dof[i]-1][0];
    #Calculate the elemental internal forces
    p_local=k_local @ u_local;
    #Get the displacements and forces and multiply to units mm,rads,kN,kN-m
    ux_1,ux_2=round(u_local[0][0],3),round(u_local[3][0],3);
    uy_1,uy_2=round(u_local[1][0],3),round(u_local[4][0],3);
    rz_1,rz_2=round(u_local[2][0],6),round(u_local[5][0],6);
    fx_1,fx_2=round(p_local[0][0]*10**-3,1),round(p_local[3][0]*10**-3,1);
    fy_1,fy_2=round(p_local[1][0]*10**-3,1),round(p_local[4][0]*10**-3,1);
    mz_1,mz_2=round(p_local[2][0]*10**-6,1),round(p_local[5][0]*10**-6,1);
    #Print the data to the console
    print('Element No: {}'.format(ele_no));
    print('Node i: {}    Node j: {}    Length: {} mm'.format(
        node_i,node_j,length));
    print('Section: {}    Material: {}'.format(section,mat));
    print('A: {}e3 mm2    Ixx: {}e6 mm4    E: {} MPa'.format(
        round(A*10**-3,1),round(Ixx*10**-6,1),round(E)));
    print('Ux.i: {} mm    Uy.i: {} mm    Rz.i: {} rads'.format(
        ux_1,uy_1,rz_1));
    print('Ux.j: {} mm    Uy.j: {} mm    Rz.j: {} rads'.format(
        ux_2,uy_2,rz_2));
    print('Fx.i: {} kN    Fy.i: {} kN    Mz.i: {} kN-m'.format(
        fx_1,fy_1,mz_1));
    print('Fx.j: {} kN    Fy.j: {} kN    Mz.j: {} kN-m\n'.format(
        fx_2,fy_2,mz_2));
    return None;

def get_element_section_cuts(ele_no,elements,cuts=5):
    """
    """
    section=elements[ele_no]['Section'];
    A=sections[section]['A'];
    Ixx=sections[section]['Ixx'];
    mat=sections[section]['Material'];
    E=materials[mat]['E'];
    length=elements[ele_no]['Length'];
    k_local=elements[ele_no]['k_local'];
    u_local=elements[ele_no]['u_local'];
    p_local=elements[ele_no]['p_local'];
    ux_1=round(u_local[0][0],3);
    uy_1=round(u_local[1][0],3);
    rz_1=round(u_local[2][0],6);
    ux_2=round(u_local[3][0],3);
    uy_2=round(u_local[4][0],3);
    rz_2=round(u_local[5][0],6);
    fx_1=round(p_local[0][0],1);
    fy_1=round(p_local[1][0],1);
    mz_1=round(p_local[2][0],1);
    fx_2=round(p_local[3][0],1);
    fy_2=round(p_local[4][0],1);
    mz_2=round(p_local[5][0],1);
    del_rz=rz_1-(bf.beam_moment_end1_deflection(-mz_1,length,E,Ixx,0)[1]+
            bf.beam_moment_end2_deflection(mz_2,length,E,Ixx,0)[1])
    for i in range(cuts+1):
        cut_x=round(i/cuts*length,1);
        fx_cut=round(-fx_1+i/cuts*(fx_2+fx_1),1);
        fy_cut=round(-fy_1+i/cuts*(fy_2+fy_1),1);
        mz_cut=round(-mz_1+i/cuts*(mz_2+mz_1),1);
        ux_cut=ux_1+i/cuts*(ux_2-ux_1);
        uy_cut=round(uy_1+i/cuts*(uy_2-uy_1)+
                bf.beam_moment_end1_deflection(-mz_1,length,E,Ixx,i/5)[0]+
                bf.beam_moment_end2_deflection(-mz_2,length,E,Ixx,i/5)[0],3);
        rz_cut=round(del_rz+bf.beam_moment_end1_deflection(-mz_1,length,E,Ixx,i/5)[1]+
                     bf.beam_moment_end2_deflection(mz_2,length,E,Ixx,i/5)[1],5);
        print('Cut {}: {} mm'.format(i+1,cut_x));
        print('Fx= {} kN    Fy= {} kN    Mz= {} kN-m'.format(
            round(fx_cut*10**-3,1),round(fy_cut*10**-3,1),
            round(mz_cut*10**-6,1)));
        print('Ux= {} mm    Uy= {} mm    Rz= {} rads'.format(
            round(ux_cut,3),round(uy_cut,3),
            round(rz_cut,6)));        
    return None;

def plot_structure(nodes,elements):
    """
    Plots the structure using the nodes and elements dictionary

    Parameters
    ----------
    nodes : dict
        nodes database.
    elements : dict
        elements database.

    Returns
    -------
    fig,ax
    """
    fig,ax=plt.subplots(); #Create the plot
    #Loop thru elements dictionary and add frames to plot
    for ele in elements:
        node_i=elements[ele]['Node_1'];
        node_j=elements[ele]['Node_2'];
        x_1,y_1=nodes[node_i]['x'],nodes[node_i]['y'];
        x_2,y_2=nodes[node_j]['x'],nodes[node_j]['y'];
        ax.plot([x_1,x_2],[y_1,y_2],
                color='black',
                linewidth=2);
    #Loop thru nodes dictionary and add nodes to the plot
    for node in nodes:
        x_1,y_1=nodes[node]['x'],nodes[node]['y'];
        rx,ry=nodes[node]['fix_x'],nodes[node]['fix_y'];
        ax.scatter(x_1,y_1,
                   c='black',
                   s=20);
        #If the node is restrained, a restraint identifier will be shown
        if(rx==True):
            ax.plot([x_1,x_1-200,x_1-200,x_1],
                    [y_1,y_1+100,y_1-100,y_1],
                    color='black',
                    linewidth=1);        
        if(ry==True):
            ax.plot([x_1,x_1-100,x_1+100,x_1],
                    [y_1,y_1-200,y_1-200,y_1],
                    color='black',
                    linewidth=1);
    ax.set_xlabel('x-Axis');
    ax.set_ylabel('y-Axis');
    ax.set_title('Structure');
    plt.grid(b=True);
    return fig,ax;

def plot_element_forces(ele_no,case_no,elements,load_cases,action='Mz'):
    """
    Plot the element forces

    Parameters
    ----------
    ele_no : int
        Element number.
    case_no : int
        Load case number.
    elements : dict
        Elements database.
    load_cases : dict
        Load Case database.
    action : str, optional
        Which action to be plotted, can be 'Px', 'Vy' or 'Mz'.
        The default is 'Mz'.

    Returns
    -------
    fig : TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    """
    action_keys={
        'Px':{'axis':'Axial Force: Px (kN)'},
        'Vy':{'axis':'Shear Force: Vy (kN)'},
        'Mz':{'axis':'Bending Moment: Mz (kN-m)'},
        }
    title='Element {} : Load Case {}'.format(
        ele_no,load_cases[case_no]['Case Name']);
    #Find the element rotation
    rot=elements[ele_no]['Rotation'];
    if(rot>-1 and rot<1):
        orientation='Beam';
    elif(rot>89 and rot<91):
        orientation='Column';
    if(rot>-1 and rot<1):
        orientation='Diagonal';
    #Create the plot
    fig,ax=plt.subplots(); 
    x_pos,ele_actions=[],[];
    #Loop the element section cuts
    for cut in elements[ele_no]['sect_cuts'][case_no]:
        x=elements[ele_no]['sect_cuts'][case_no][cut]['Offset'];
        f=elements[ele_no]['sect_cuts'][case_no][cut][action];
        x_pos+=[x];
        ele_actions+=[f];
        #Plot the section cut lines
        if(orientation=='Column'):
            ax.plot([0,f],[x,x],
                    color='red',
                    linewidth=0.5);
        else:
            ax.plot([x,x],[0,f],
                    color='red',
                    linewidth=0.5);
    #Plot the element forces
    if(orientation=='Column'):
        ax.plot(ele_actions,x_pos,
                color='red',
                linewidth=2);
        ax.plot([0,0],[x_pos[0],x_pos[-1]],
                color='black',
                linewidth=2);
        ax.scatter([0,0],[x_pos[0],x_pos[-1]],
                   c='black',
                   s=20);
        ax.set_ylabel('Element x-axis (mm)');
        ax.set_xlabel(action_keys[action]['axis']);
        if(action=='Mz'):
            ax.invert_xaxis();
        plt.grid(b=True,
                 axis='x',
                 linestyle=':');         
    else:
        ax.plot(x_pos,ele_actions,
                color='red',
                linewidth=2);
        ax.plot([x_pos[0],x_pos[-1]],[0,0],
                color='black',
                linewidth=2);
        ax.scatter([x_pos[0],x_pos[-1]],[0,0],
                   c='black',
                   s=20);
        ax.set_xlabel('Element x-axis (mm)');
        ax.set_ylabel(action_keys[action]['axis']);
        if(action=='Mz'):
            ax.invert_yaxis();
        plt.grid(b=True,
                 axis='y',
                 linestyle=':');
    #Set the plot title
    ax.set_title(title);
    return fig,ax;

def plot_element_displacement(ele_no,case_no,elements,load_cases,action='Uy'):
    """
    Plot the element displacement

    Parameters
    ----------
    ele_no : int
        Element number.
    case_no : int
        Load case number.
    elements : dict
        Elements database.
    load_cases : dict
        Load Case database.
    action : str, optional
        Which displacement is to be plotted, can be 'Uy'.
        The default is 'Uy'.

    Returns
    -------
    fig : TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    """
    action_keys={
        'Ux':{'axis':'X-Displacement: Ux (mm)'},
        'Uy':{'axis':'Y-Displacement: Uy (mm)'},
        'Rz':{'axis':'Z-Rotation: Rz (rad)'},
        }
    title='Element {} : Load Case {}'.format(
        ele_no,load_cases[case_no]['Case Name']);
    #Find the element rotation
    rot=elements[ele_no]['Rotation'];
    if(rot>-1 and rot<1):
        orientation='Beam';
    elif(rot>89 and rot<91):
        orientation='Column';
    if(rot>-1 and rot<1):
        orientation='Diagonal';
    #Create the plot
    fig,ax=plt.subplots(); 
    x_pos,ele_disp=[],[];
    #Loop the element section cuts
    for cut in elements[ele_no]['sect_cuts'][case_no]:
        x=elements[ele_no]['sect_cuts'][case_no][cut]['Offset'];
        u=elements[ele_no]['sect_cuts'][case_no][cut][action];
        x_pos+=[x];
        ele_disp+=[u];
        #Plot the section cut lines
        if(orientation=='Column'):
            ax.plot([0,u],[x,x],
                    color='red',
                    linewidth=0.5);
        else:
            ax.plot([x,x],[0,u],
                    color='red',
                    linewidth=0.5);
    #Plot the element forces
    if(orientation=='Column'):
        ax.plot(ele_disp,x_pos,
                color='red',
                linewidth=2);
        ax.plot([0,0],[x_pos[0],x_pos[-1]],
                color='black',
                linewidth=2);
        ax.scatter([0,0],[x_pos[0],x_pos[-1]],
                   c='black',
                   s=20);
        ax.set_ylabel('Element x-axis (mm)');
        ax.set_xlabel(action_keys[action]['axis']);
        plt.grid(b=True,
                 axis='x',
                 linestyle=':');         
    else:
        ax.plot(x_pos,ele_disp,
                color='red',
                linewidth=2);
        ax.plot([x_pos[0],x_pos[-1]],[0,0],
                color='black',
                linewidth=2);
        ax.scatter([x_pos[0],x_pos[-1]],[0,0],
                   c='black',
                   s=20);
        ax.set_xlabel('Element x-axis (mm)');
        ax.set_ylabel(action_keys[action]['axis']);
        plt.grid(b=True,
                 axis='y',
                 linestyle=':');
    #Set the plot title
    ax.set_title(title);
    return fig,ax;

def print_element_node_loads(ele_no,elements,case_no):
    """
    Print the element node loads to the console to easily check the node 
    loads for an element and load case
    
    Parameters
    ----------
    ele_no : int
        The element numbers.
    elements : dict
        The elements database stored in a dictionary. Analysis must be 
        performed so that node loads have been calculated and stored
    case_no : int
        The load case number.

    Returns
    -------
    None.
    """
    #Get elements transformation matrix
    transform=elements[ele_no]['T_local'];
    #Get the elements length
    length=elements[ele_no]['Length'];
    #Get the elements uniform load in the local x and y-direction
    w_x,w_y=elements[ele_no]['loads'][case_no];
    #Get the global node loads
    u_local_g=transform.T @ frame_eq_node_loads(w_x,w_y,length)
    #Get the loads for elements individual DOF
    px_i=round(u_local_g[0][0]*10**-3,1);
    py_i=round(u_local_g[1][0]*10**-3,1);
    mz_i=round(u_local_g[2][0]*10**-6,1);
    px_j=round(u_local_g[3][0]*10**-3,1);
    py_j=round(u_local_g[4][0]*10**-3,1);
    mz_j=round(u_local_g[5][0]*10**-6,1);
    #Print the output
    print('ELEMENT {} : CASE {}'.format(ele_no,case_no));
    print('Px.i = {} kN'.format(px_i));
    print('Py.i = {} kN'.format(py_i));
    print('Mz.i = {} kN-m'.format(mz_i));
    print('Px.j = {} kN'.format(px_j));
    print('Py.j = {} kN'.format(py_j));
    print('Mz.j = {} kN-m'.format(mz_j));
    return None;



# #Define nodes
# nodes={
#     1:{'x':0.0,   'y':0.0,'fix_x':True,'fix_y':True,'fix_mz':False},
#     2:{'x':1000.0,'y':0.0,'fix_x':False,'fix_y':False,'fix_mz':False},
#     3:{'x':2000.0,'y':0.0,'fix_x':False,'fix_y':False,'fix_mz':False},
#     4:{'x':3000.0,'y':0.0,'fix_x':False,'fix_y':False,'fix_mz':False},
#     5:{'x':4000.0,'y':0.0,'fix_x':False,'fix_y':False,'fix_mz':False},
#     6:{'x':5000.0,'y':0.0,'fix_x':False,'fix_y':False,'fix_mz':False},
#     7:{'x':6000.0,'y':0.0,'fix_x':True,'fix_y':True,'fix_mz':False}
#     };

# #Define the materials       
# materials={
#     'CONC_32':{'E':32000.0,'Density':2500.0}
#     };

# #Define the sections
# sections={
#     'B200x1000G32':{
#         'A':200000.0,
#         'Ixx':666666666.667,
#         'Material':'CONC_32'
#         }
#     };

# #Define the elements
# elements={
#     1:{'Node_1':1,'Node_2':2,'Section':'B200x1000G32'},
#     2:{'Node_1':2,'Node_2':3,'Section':'B200x1000G32'},
#     3:{'Node_1':3,'Node_2':4,'Section':'B200x1000G32'},
#     4:{'Node_1':4,'Node_2':5,'Section':'B200x1000G32'},
#     5:{'Node_1':5,'Node_2':6,'Section':'B200x1000G32'},
#     6:{'Node_1':6,'Node_2':7,'Section':'B200x1000G32'}
#     };

# #Define the global freedoms
# global_freedoms={
#     'x-translation':False,
#     'y-translation':True,
#     'z-rotation':True
#     };

# #Define the load cases. Gravity is in m/s2
# load_cases={
#     1:{'Case Name':'LOAD CASE'},
#     }

# #Define the frame loads
# frame_loads={};

# #Define the loads. Node loads are in kN
# node_loads={
#     1:{'LOAD CASE':[0,-5]},
#     2:{'LOAD CASE':[0,-10]},
#     3:{'LOAD CASE':[0,-10]},
#     4:{'LOAD CASE':[0,-10]},
#     5:{'LOAD CASE':[0,-10]},
#     6:{'LOAD CASE':[0,-10]},
#     7:{'LOAD CASE':[0,-5]},
#     };

#Define nodes. Cartesian co-ordinates in mm
nodes={
    1:{'x':400.0,   'y':0.0,'fix_x':True,'fix_y':True,'fix_mz':False},
    2:{'x':7600.0,'y':0.0,'fix_x':True,'fix_y':True,'fix_mz':False},
    3:{'x':0.0,'y':3390.0,'fix_x':False,'fix_y':False,'fix_mz':False},
    4:{'x':400.0,'y':3390.0,'fix_x':False,'fix_y':False,'fix_mz':False},
    5:{'x':7600.0,'y':3390.0,'fix_x':False,'fix_y':False,'fix_mz':False},
    6:{'x':8000.0,'y':3390.0,'fix_x':False,'fix_y':False,'fix_mz':False},
    7:{'x':0.0,'y':6890.0,'fix_x':False,'fix_y':False,'fix_mz':False},
    8:{'x':400.0,'y':6890.0,'fix_x':False,'fix_y':False,'fix_mz':False},
    9:{'x':7600.0,'y':6890.0,'fix_x':False,'fix_y':False,'fix_mz':False},
    10:{'x':8000.0,'y':6890.0,'fix_x':False,'fix_y':False,'fix_mz':False},
    };

#Define the materials       
materials={
    'CONC_32':{'E':32000.0,'Density':2500.0},
    'CONC_50':{'E':34500.0,'Density':2500.0}
    };

#Define the sections
sections={
    'B200x4000G32':{
        'A':800000.0,
        'Ixx':2666666666.6,
        'Material':'CONC_32'},
    'C800x350G50':{
        'A':280000.0,
        'Ixx':14933333333.3,
        'Material':'CONC_50'}
    };

#Define the elements
elements={
    1:{'Node_1':1,'Node_2':4,'Section':'C800x350G50'},
    2:{'Node_1':2,'Node_2':5,'Section':'C800x350G50'},
    3:{'Node_1':3,'Node_2':4,'Section':'B200x4000G32'},
    4:{'Node_1':4,'Node_2':5,'Section':'B200x4000G32'},
    5:{'Node_1':5,'Node_2':6,'Section':'B200x4000G32'},
    6:{'Node_1':4,'Node_2':8,'Section':'C800x350G50'},
    7:{'Node_1':5,'Node_2':9,'Section':'C800x350G50'},
    8:{'Node_1':7,'Node_2':8,'Section':'B200x4000G32'},
    9:{'Node_1':8,'Node_2':9,'Section':'B200x4000G32'},
    10:{'Node_1':9,'Node_2':10,'Section':'B200x4000G32'},
    };

#Define the global freedoms
global_freedoms={
    'x-translation':True,
    'y-translation':True,
    'z-rotation':True
    };

#Define the load cases. Gravity is in m/s2
load_cases={
    1:{'Case Name':'DL','gravity':[0,-9.81]},
    2:{'Case Name':'SDL'},
    3:{'Case Name':'LL'},
    4:{'Case Name':'EQ'}
    }

#Define the frame area loads. Frame loads are applied to the local axis of
#of the element
frame_loads={
    3:{'SDL':[0,-4],'LL':[0,-12]},
    4:{'SDL':[0,-4],'LL':[0,-12]},
    5:{'SDL':[0,-4],'LL':[0,-12]},
    8:{'SDL':[0,-4],'LL':[0,-12]},
    9:{'SDL':[0,-4],'LL':[0,-12]},
    10:{'SDL':[0,-4],'LL':[0,-12]},
    };

#Define the loads. Node loads are in kN
node_loads={
    3:{'EQ':[50,0]},
    5:{'EQ':[50,0]},
    8:{'EQ':[50,0]},
    10:{'EQ':[50,0]},
    };

#Starting Analysis
start_time=time.time();
print('BEGINNING ANALYSIS');

#build the dof dictionary and add the dof to the nodes matrix
dof_db={};
dof=0;
for node in nodes:
    x_t_restraint=not(global_freedoms['x-translation']) or nodes[node]['fix_x'];
    y_t_restraint=not(global_freedoms['y-translation']) or nodes[node]['fix_y'];
    z_r_restraint=not(global_freedoms['z-rotation']) or nodes[node]['fix_mz'];
    dof+=1;
    dof_db[dof]={
        'node':node,
        'freedom':'x-translation',
        'restraint':x_t_restraint,
        #'load':Fx_load};
        }
    nodes[node]['Ux_dof']=dof;
    dof+=1;
    dof_db[dof]={
        'node':node,
        'freedom':'y-translation',
        'restraint':y_t_restraint,
        #'load':Fy_load};
        };
    nodes[node]['Uy_dof']=dof;
    dof+=1;
    dof_db[dof]={
        'node':node,
        'freedom':'z-rotation',
        'restraint':z_r_restraint,
        #'load':Mz_load};
        };
    nodes[node]['Rz_dof']=dof;
print('Structure consists {} Nodes, {} Elements, {} Degrees of Freedom'.format(
    len(nodes),len(elements),len(dof_db)));

#Delete the unused variables
del(dof,node,x_t_restraint,y_t_restraint,z_r_restraint,)
# calculate the elemental stiffness matrix and add to the elements matrix
for ele in elements:
    node_1=elements[ele]['Node_1'];
    node_2=elements[ele]['Node_2'];
    #Get the cartesian coordinates of the elements nodes
    x_1,y_1=nodes[node_1]['x'],nodes[node_1]['y'];
    x_2,y_2=nodes[node_2]['x'],nodes[node_2]['y'];
    #Get the element section
    section=elements[ele]['Section'];
    #Calculate the element length
    length=round(((x_2-x_1)**2+(y_2-y_1)**2)**0.5,3);
    #Calculate the elements orientation
    cos_x=(x_2-x_1)/length;
    cos_y=(y_2-y_1)/length;
    rot=m.degrees(m.acos(cos_x));
    #Get the elements parameters
    A=sections[section]['A'];
    Ix=sections[section]['Ixx'];
    material=sections[section]['Material'];
    E=materials[material]['E'];
    #Calculate the element stiffness matrix and transformation matrix
    ele_transform=transformation_6(cos_x,cos_y);
    k_local_0=beam_element_stiff(A,Ix,E,length);
    k_local=ele_transform.T @ k_local_0 @ ele_transform;
    #Add the element length and stiffness matrix to the elements dictionary
    elements[ele]['Length']=length;
    elements[ele]['Rotation']=rot;
    elements[ele]['T_local']=ele_transform;
    elements[ele]['k_local']=k_local;
    elements[ele]['loads']={};
print('Processed Elemental Properties');
#Delete the unused variables
del(A,E,Ix,cos_x,cos_y,ele,ele_transform,k_local,k_local_0,length,material,
    node_1,node_2,section,x_1,x_2,y_1,y_2);

#Initiate the global stiffness matrix
k_global=np.zeros((len(dof_db),len(dof_db)));
#Populate the local and global stiffness matrix
for ele in elements:
    #Get the elements nodes
    node_1=elements[ele]['Node_1'];
    node_2=elements[ele]['Node_2'];
    #Get the degrees of freedoms for the node
    element_dof=[nodes[node_1]['Ux_dof'],
                 nodes[node_1]['Uy_dof'],
                 nodes[node_1]['Rz_dof'],
                 nodes[node_2]['Ux_dof'],
                 nodes[node_2]['Uy_dof'],
                 nodes[node_2]['Rz_dof']];
    #Get the elemental stiffness matrix
    k_local=elements[ele]['k_local'];
    #Add the element stiffness components to the global stiffness matrix
    for i in range(6):
        i_global=element_dof[i]-1;
        for j in range(6):
            j_global=element_dof[j]-1;
            k_global[i_global][j_global]+=k_local[i][j];
print('Built Global stiffness matrix Size: {} kB'.format(
    round(sys.getsizeof(k_global)/1024,1)));

#Delete the unused variables
del(ele,k_local,element_dof,i,i_global,j,j_global,node_1,node_2)

#Create a list of unrestrained degrees of freedom
dofs=[];
for freedom in dof_db:
    if(dof_db[freedom]['restraint']==False):
        dofs+=[freedom];

#Initiate the reduced stiffness matrix
k_global_red=np.zeros((len(dofs),len(dofs)));
#Get the reduced stiffness matrix
for i in range(len(dofs)):
    global_i=dofs[i];
    for j in range(len(dofs)):
        global_j=dofs[j];
        k_global_red[i][j]=k_global[global_i-1][global_j-1];

#Initiate the global loads matrix
# p_global=np.zeros((len(dof_db),1));
p_global={};
for case in load_cases:
    p_global[case]=np.zeros((len(dof_db),1));

#Add the gravity load to p_global
for case in load_cases:
    if(load_cases[case].get('gravity',False)!=False):
        g_x,g_y=load_cases[case]['gravity'];
        for ele in elements:
            #Get nodes
            node_i=elements[ele]['Node_1'];
            node_j=elements[ele]['Node_2'];
            #Get the global dofs
            ele_dofs=[nodes[node_i]['Ux_dof'],nodes[node_i]['Uy_dof'],
                      nodes[node_i]['Rz_dof'],nodes[node_j]['Ux_dof'],
                      nodes[node_j]['Uy_dof'],nodes[node_j]['Rz_dof']];
            #Get the section properties
            section=elements[ele]['Section'];
            rot=elements[ele]['Rotation'];
            length=elements[ele]['Length'];
            mat=sections[section]['Material'];
            A=sections[section]['A'];
            density=materials[mat]['Density'];
            cos_x=elements[ele]['T_local'][0][0];
            cos_y=elements[ele]['T_local'][0][1];
            transform_mat=elements[ele]['T_local'];
            #Get the section global loads
            w_xg=g_x*A*density*10**-9;
            w_yg=g_y*A*density*10**-9;
            #local uniform distributed loads
            w_x=cos_x*w_xg+cos_y*w_yg;
            w_y=-cos_y*w_xg+cos_x*w_yg;
            #get the elements node loads
            n_loads=transform_mat.T @ frame_eq_node_loads(w_x,w_y,length);
            #Frame loads
            if(elements[ele].get('loads',False)==False):
                elements[ele]['loads']={case:[round(w_x,6),round(w_y,6)]};
            elif(elements[ele]['loads'].get(case,False)==False):
                elements[ele]['loads'][case]=[round(w_x,6),round(w_y,6)];
            else:
                elements[ele]['loads'][case][0]+=round(w_x,6);
                elements[ele]['loads'][case][1]+=round(w_y,6);
            for i in range(len(ele_dofs)):
                p_global[case][ele_dofs[i]-1][0]+=round(n_loads[i][0],6);

#Add the frame loads to the p_global matrix
for ele in frame_loads:
    for case_name in frame_loads[ele]:
        case=False;
        for case_nm in load_cases:
            if(load_cases[case_nm]['Case Name']==case_name):
                case=case_nm;
        if(case!=False):
            w_x=frame_loads[ele][case_name][0];
            w_y=frame_loads[ele][case_name][1];
            length=elements[ele]['Length'];
            transform_mat=elements[ele]['T_local'];
            cos_x=elements[ele]['T_local'][0][0];
            cos_y=elements[ele]['T_local'][0][1];
            n_loads=transform_mat.T @ frame_eq_node_loads(w_x,w_y,length);
            #Get nodes
            node_i=elements[ele]['Node_1'];
            node_j=elements[ele]['Node_2'];
            #Get the global dofs
            ele_dofs=[nodes[node_i]['Ux_dof'],nodes[node_i]['Uy_dof'],
                      nodes[node_i]['Rz_dof'],nodes[node_j]['Ux_dof'],
                      nodes[node_j]['Uy_dof'],nodes[node_j]['Rz_dof']];
            for i in range(len(ele_dofs)):
                p_global[case][ele_dofs[i]-1][0]+=round(n_loads[i][0],6);
            if(elements[ele].get('loads',False)==False):
                elements[ele]['loads']={case:[round(w_x,6),round(w_y,6)]};
            elif(elements[ele]['loads'].get(case,False)==False):
                elements[ele]['loads'][case]=[round(w_x,6),round(w_y,6)];
            else:
                elements[ele]['loads'][case][0]+=round(w_x,6);
                elements[ele]['loads'][case][1]+=round(w_y,6);

#Put the node loads into the p_global matrix
for node in node_loads:
    for case_name in node_loads[node]:
        case=False;
        for case_nm in load_cases:
            if(load_cases[case_nm]['Case Name']==case_name):
                case=case_nm;
        if(case!=False):
            f_x=node_loads[node][case_name][0]*10**3;
            f_y=node_loads[node][case_name][1]*10**3;
            #Get the global dofs
            node_dofs=[nodes[node]['Ux_dof'],nodes[node]['Uy_dof']];
            #Put the node loads into the p_global
            p_global[case][node_dofs[0]-1][0]+=round(f_x,6);
            p_global[case][node_dofs[1]-1][0]+=round(f_y,6);

#Initiate the reduced loads matrix
p_global_red={};
for case in load_cases:
    p_global_red[case]=np.zeros((len(dofs),1));

#Get the reduced loads matrix
for case in load_cases:
    for i in range(len(dofs)):
        global_i=dofs[i];
        p_global_red[case][i][0]=p_global[case][global_i-1][0];

print('Built Reduced global loads matrix Size: {} kB'.format(
    round(sys.getsizeof(p_global_red)/1024,1)));

#Calculate the local displacements
u_global_red={};
for case in load_cases:
    u_global_red[case]=np.linalg.inv(k_global_red) @ p_global_red[case];
print('Solved Reduced global displacement matrix Size: {} kB'.format(
    round(sys.getsizeof(u_global_red)/1024,1)));

#Populate the global displacements matrix
u_global={};
for case in load_cases:
    u_global[case]=np.zeros((len(dof_db),1));
    for i in range(len(dofs)):
        u_global[case][dofs[i]-1]=u_global_red[case][i];

#Generate the global loads matrix
p_global={};
for case in load_cases:
    p_global[case]=k_global@u_global[case];

#Post-process Elements
for ele in elements:
    node_i=elements[ele]['Node_1'];
    node_j=elements[ele]['Node_2'];
    k_local=elements[ele]['k_local'];
    ele_dofs=[nodes[node_i]['Ux_dof'],nodes[node_i]['Uy_dof'],nodes[node_i]['Rz_dof'],
              nodes[node_j]['Ux_dof'],nodes[node_j]['Uy_dof'],nodes[node_j]['Rz_dof']];
    u_local={};
    for case in load_cases:
        u_local[case]=np.zeros((6,1));
        for i in range(len(ele_dofs)):
            u_local[case][i]=u_global[case][ele_dofs[i]-1][0];
    p_local={};
    for case in load_cases:
        p_local[case]=k_local@u_local[case];
    elements[ele]['p_local']=p_local;
    elements[ele]['u_local']=u_local;
del(ele,node_i,node_j,k_local,dofs,u_local,i,p_local);

#Post-process elements force
for ele in elements:
    section=elements[ele]['Section'];
    length=elements[ele]['Length'];
    A=sections[section]['A'];
    Ix=sections[section]['Ixx'];
    mat=sections[section]['Material'];
    E=materials[mat]['E'];
    transform=elements[ele]['T_local'];
    for case in load_cases:
        #Get the loads and the local forces and displacements
        w_loads=elements[ele]['loads'].get(case,[0,0]);
        p_local=transform @ elements[ele]['p_local'][case];
        u_local=transform @ elements[ele]['u_local'][case];
        #Get the applied equivalent node loads
        p_applied=frame_eq_node_loads(w_loads[0],w_loads[1],length);
        #Reapply equivalent node loads to p_local
        p_local=p_local-p_applied;
        #Update the local element forces to the elements database
        elements[ele]['p_local'][case]=p_local;
        elements[ele]['u_local'][case]=u_local;
        #Get the local forces
        pi,pj=-p_local[0][0],p_local[3][0];
        vi,vj=p_local[1][0],-p_local[4][0];
        mi,mj=-p_local[2][0],p_local[5][0];
        uxi,uxj=u_local[0][0],u_local[3][0];
        uyi,uyj=u_local[1][0],u_local[4][0];
        rzi,rzj=u_local[2][0],u_local[5][0];
        wx,wy=w_loads[0],w_loads[1];
        section_cuts={};
        for i in range(9):
            x_ratio=round(i/8,3);
            x_cut=round(i*length/8,1);
            Px,Ux=bf.bar_formula(pi,uxi,E,A,x_ratio,length,wx);
            Vy,Mz,Rz,Uy=bf.beam_formula(vi,mi,uyi,rzi,E,Ix,x_ratio,length,wy);
            section_cuts[i]={
                'Offset':x_cut,
                'Px':round(Px*10**-3,4),
                'Vy':round(Vy*10**-3,4),
                'Mz':round(Mz*10**-6,4),
                'Ux':round(Ux,4),
                'Uy':round(Uy,4),
                'Rz':round(Rz,6),
                };
        if(elements[ele].get('sect_cuts',False)==False):
            elements[ele]['sect_cuts']={};
        elements[ele]['sect_cuts'][case]=section_cuts;
        

#Print results
end_time=time.time();
print('COMPLETED ANALYSIS in {} secs\n'.format(round(end_time-start_time,3)));

#Delete unused varibales
# del(A,case,case_name,cos_x,cos_y,density,ele_dofs,f_x,f_y,freedom,global_i,
#     global_j,j,length,mat,n_loads,node,node_dofs,rot,start_time,end_time,
#     transform_mat,w_x,w_xg,w_y,w_yg);

#Print summary to panel
for case in load_cases:
    print('LOAD CASE: {}- {}\n'.format(case,load_cases[case]['Case Name']));
    for node in nodes:
        node_x,node_y=round(nodes[node]['x'],1),round(nodes[node]['y'],1);
        node_x_fix,node_x_dof=nodes[node]['fix_x'],nodes[node]['Ux_dof'];
        node_y_fix,node_y_dof=nodes[node]['fix_y'],nodes[node]['Uy_dof'];
        node_mz_fix,node_mz_dof=nodes[node]['fix_mz'],nodes[node]['Rz_dof'];
        print('NODE: {} \nx: {} mm\ny: {} mm'.format(node,node_x,node_y));
        if(node_x_fix):
            node_reaction=round(p_global[case][node_x_dof-1][0]*10**-3,3);
            print('Reaction RX: {} kN'.format(node_reaction));
        if(node_y_fix):
            node_reaction=round(p_global[case][node_y_dof-1][0]*10**-3,3);
            print('Reaction RY: {} kN'.format(node_reaction));
        if(node_mz_fix):
            node_reaction=round(p_global[case][node_mz_dof-1][0]*10**-6,1);
            print('Reaction MZ: {} kN-m'.format(node_reaction));
        node_disp=round(u_global[case][node_x_dof-1][0],3);
        print('Displacement x: {} mm'.format(node_disp));
        node_disp=round(u_global[case][node_y_dof-1][0],3);
        print('Displacement y: {} mm'.format(node_disp));
        node_disp=round(u_global[case][node_mz_dof-1][0],5);
        print('Rotation Rz: {} rads\n'.format(node_disp));

#Delete the unused variables
del(node,node_disp,node_mz_dof,node_mz_fix,node_reaction,node_x,node_x_dof,
    node_x_fix,node_y,node_y_dof,node_y_fix);

#Print the elements summary
for case in load_cases:
    print('LOAD CASE: {}- {}\n'.format(case,load_cases[case]['Case Name']))
    for ele in elements:
        node_i=elements[ele]['Node_1'];
        node_j=elements[ele]['Node_2'];
        section=elements[ele]['Section'];
        length=elements[ele]['Length'];
        px_1=round(elements[ele]['p_local'][case][0][0]*10**-3,1);
        py_1=round(elements[ele]['p_local'][case][1][0]*10**-3,1);
        mz_1=round(elements[ele]['p_local'][case][2][0]*10**-6,1);
        px_2=round(elements[ele]['p_local'][case][3][0]*10**-3,1);
        py_2=round(elements[ele]['p_local'][case][4][0]*10**-3,1);
        mz_2=round(elements[ele]['p_local'][case][5][0]*10**-6,1);
        ux_1=round(elements[ele]['u_local'][case][0][0],3);
        uy_1=round(elements[ele]['u_local'][case][1][0],3);
        rz_1=round(elements[ele]['u_local'][case][2][0],6);
        ux_2=round(elements[ele]['u_local'][case][3][0],3);
        uy_2=round(elements[ele]['u_local'][case][4][0],3);
        rz_2=round(elements[ele]['u_local'][case][5][0],6);
        print('Element {}:'.format(ele));
        print('Node 1: {} ; Node 2: {}'.format(node_i,node_j));
        print('Section: {}'.format(section));
        print('Length: {}'.format(length));
        print('Memeber End Displacements');
        print('Ux_1: {} mm ; Uy_1: {} mm ; Rz_1: {} rads'.format(ux_1,uy_1,rz_1));
        print('Ux_2: {} mm ; Uy_2: {} mm ; Rz_2: {} rads'.format(ux_2,uy_2,rz_2));
        print('Memeber End Forces');
        print('Px_1: {} kN ; Py_1: {} kN ; Mz_1: {} kN-m'.format(px_1,py_1,mz_1));
        print('Px_2: {} kN ; Py_2: {} kN ; Mz_2: {} kN-m\n'.format(px_2,py_2,mz_2));

#Delete unused variables
# del(case,case_nm,g_x,g_y,ele,node_i,node_j,section,length,px_1,py_1,mz_1,px_2,py_2,mz_2,
#     ux_1,uy_1,rz_1,ux_2,uy_2,rz_2);

#Create a structure dictionary where we can saved pyplot data
str_fig,str_ax=plot_structure(nodes,elements);

structure={
    'structure_plot':[str_fig,str_ax]
    };


# for i in range(1,8):
#     print_nodal_data(i,nodes,node_loads,p_global,u_global);
# del(i);

# for i in range(1,7):
#     print_element_data(i,elements,nodes,p_global,u_global);
# del(i);

# get_element_section_cuts(3,elements,cuts=5)
