# -*- coding: utf-8 -*-
from decimal import Decimal;
import math as m;
import numpy as np;
import sys;
import sympy;
import time;

#Define the global freedoms
global_freedoms={
    'z-translation':True,
    'x-rotation':True,
    'y-rotation':True
    };

#Node data = [x,y,Rx,RMx,RMy]
nodes={
   1:{'x':0.0,'y':0.0,'fix_z':True,'fix_mx':False,'fix_my':False},
   2:{'x':1000.0,'y':0.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   3:{'x':2000.0,'y':0.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   4:{'x':3000.0,'y':0.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   5:{'x':4000.0,'y':0.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   6:{'x':5000.0,'y':0.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   7:{'x':6000.0,'y':0.0,'fix_z':True,'fix_mx':False,'fix_my':False},
   8:{'x':0.0,'y':1000.0,'fix_z':True,'fix_mx':False,'fix_my':False},
   9:{'x':1000.0,'y':1000.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   10:{'x':2000.0,'y':1000.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   11:{'x':3000.0,'y':1000.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   12:{'x':4000.0,'y':1000.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   13:{'x':5000.0,'y':1000.0,'fix_z':False,'fix_mx':False,'fix_my':False},
   14:{'x':6000.0,'y':1000.0,'fix_z':True,'fix_mx':False,'fix_my':False}
   };

#Define the materials       
materials={
    'CONC_32':{'E':32000,'Poisson':0.25}
    };

#Define the section
sections={
    'S200G32':{'Thickness':200,'Material':'CONC_32'}
    };

#element data = [node1,node2,node3,node4,section]
elements={
    1:{'Node_i':1,'Node_j':2,'Node_k':9,'Node_l':8,'Section':'S200G32'},
    2:{'Node_i':2,'Node_j':3,'Node_k':10,'Node_l':9,'Section':'S200G32'},
    3:{'Node_i':3,'Node_j':4,'Node_k':11,'Node_l':10,'Section':'S200G32'},
    4:{'Node_i':4,'Node_j':5,'Node_k':12,'Node_l':11,'Section':'S200G32'},
    5:{'Node_i':5,'Node_j':6,'Node_k':13,'Node_l':12,'Section':'S200G32'},
    6:{'Node_i':6,'Node_j':7,'Node_k':14,'Node_l':13,'Section':'S200G32'},
    };

elements_data=[[1,2,9,8,'S200G40'],
               [2,3,10,9,'S200G40'],
               [3,4,11,10,'S200G40'],
               [4,5,12,11,'S200G40'],
               [5,6,13,12,'S200G40'],
               [6,7,14,13,'S200G40']];

#Define the loads
node_loads={1:[-2500.0,0,0],
            2:[-5000.0,0,0],
            3:[-5000.0,0,0],
            4:[-5000.0,0,0],
            5:[-5000.0,0,0],
            6:[-5000.0,0,0],
            7:[-2500.0,0,0],
            8:[-2500.0,0,0],
            9:[-5000.0,0,0],
            10:[-5000.0,0,0],
            11:[-5000.0,0,0],
            12:[-5000.0,0,0],
            13:[-5000.0,0,0],
            14:[-2500.0,0,0]
            };

#print the size of nodes storage
print('{} nodes, node storage size = {} kB'.format(len(nodes),
    round(sys.getsizeof(nodes)/1024,1)));

#print the size of materials storage
print('{} materials, materials storage size = {} kB'.format(len(materials),
    round(sys.getsizeof(materials)/1024,1)));

#print the size of sections storage
print('{} sections, materials storage size = {} kB'.format(len(sections),
    round(sys.getsizeof(sections)/1024,1)));

#print the size of elements storage
print('{} elements, elements storage size = {} kB'.format(len(elements),
    round(sys.getsizeof(elements)/1024,1)));

#print the size of node loads storage
print('{} node loads, node load storage size = {} kB'.format(len(node_loads),
    round(sys.getsizeof(node_loads)/1024,1)));

def B_matrix_plate_mzc(x_len,y_len,x_g,y_g):
    """
    B_matrix is formulated on plate mzc theory
    
    Currently only handles rectangular elements in the x-y plane. The x and y
    gauss points are the x and y centroid and should be 0,0 for the localized
    element

    Parameters
    ----------
    x_len : float
        length in the x-direction
    y_len : float
        length in the y-direction
    x_g : float
        x gauss point
    y_g : float
        y gauss point

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    d2N=np.zeros((4,4));
    d2NN=np.zeros((4,4));
    d2NNN=np.zeros((4,4));

    d2N[0,0]=3*(x_g-x_g*y_g)/(4*x_len**2);
    d2N[1,0]=3*(-x_g+x_g*y_g)/(4*x_len**2);
    d2N[2,0]=3*(-x_g-x_g*y_g)/(4*x_len**2);
    d2N[3,0]=3*(x_g+x_g*y_g)/(4*x_len**2);

    d2N[0,1]=3*(y_g-x_g*y_g)/(4*y_len**2);
    d2N[1,1]=3*(y_g+x_g*y_g)/(4*y_len**2);
    d2N[2,1]=3*(-y_g-x_g*y_g)/(4*y_len**2);
    d2N[3,1]=3*(-y_g+x_g*y_g)/(4*y_len**2);

    d2N[0,2]=2*(1/2-3*x_g**2/8-3*y_g**2/8)/(x_len*y_len)
    d2N[1,2]=2*(-1/2+3*x_g**2/ 8 + 3 * y_g** 2 / 8) / (x_len*y_len)
    d2N[2,2]=2*(1/2-3*x_g**2/ 8 - 3 * y_g** 2 / 8) / (x_len*y_len)
    d2N[3,2]=2*(-1/2+3*x_g**2/ 8 + 3 * y_g** 2 / 8) / (x_len*y_len)

    # = == == == == == == == == =

    d2NN[0,0]=((3*x_len*x_g-3*x_len*x_g*y_g-x_len+x_len*y_g)/4)/x_len**2;
    d2NN[1,0]=((3*x_len*x_g-3*x_len*x_g*y_g+x_len-x_len*y_g)/4)/x_len**2;
    d2NN[2,0]=((3*x_len*x_g+3*x_len*x_g*y_g+x_len+x_len*y_g)/4)/x_len**2;
    d2NN[3,0]=((3*x_len*x_g+3*x_len*x_g*y_g-x_len-x_len*y_g)/4)/x_len**2;

    d2NN[0,1]=0;
    d2NN[1,1]=0;
    d2NN[2,1]=0;
    d2NN[3,1]=0;

    d2NN[0,2]=2*(-3/8*x_len*x_g**2+x_len*x_g/4+x_len/8)/(x_len*y_len);
    d2NN[1,2]=2*(-3/8*x_len*x_g**2-x_len*x_g/4+x_len/8)/(x_len*y_len);
    d2NN[2,2]=2*(3/8*x_len*x_g**2+x_len*x_g/4-x_len/8)/(x_len*y_len);
    d2NN[3,2]=2*(3/8*x_len*x_g**2-x_len*x_g/4-x_len/8)/(x_len*y_len);

    # = == == == == == == == == =

    d2NNN[0,0]=0;
    d2NNN[1,0]=0;
    d2NNN[2,0]=0;
    d2NNN[3,0]=0;

    d2NNN[0,1]=((3*y_len*y_g-3*y_len*x_g*y_g-y_len+y_len*x_g)/4)/y_len**2;
    d2NNN[1,1]=((3*y_len*y_g+3*y_len*x_g*y_g-y_len-y_len*x_g)/4)/y_len**2;
    d2NNN[2,1]=((3*y_len*y_g+3*y_len*x_g*y_g+y_len+y_len*x_g)/4)/y_len**2;
    d2NNN[3,1]=((3*y_len*y_g-3*y_len*x_g*y_g+y_len-y_len*x_g)/4)/y_len**2;

    d2NNN[0,2]=2*(-3/8*y_len*y_g**2+y_len*y_g/4+y_len/8)/(x_len*y_len);
    d2NNN[1,2]=2*(3/8*y_len*y_g**2-y_len*y_g/4-y_len/8)/(x_len*y_len);
    d2NNN[2,2]=2*(3/8*y_len*y_g**2+y_len*y_g/4-y_len/8)/(x_len*y_len);
    d2NNN[3,2]=2*(-3/8*y_len*y_g**2-y_len*y_g/4+y_len/8)/(x_len*y_len);

    # = == == == == == == == == =

    bmat_1=np.array([[-d2N[0,0],-d2NN[0,0],-d2NNN[0,0],-d2N[1,0],
                      -d2NN[1,0],-d2NNN[1,0],-d2N[2,0],-d2NN[2,0],
                      -d2NNN[2,0],-d2N[3,0],-d2NN[3,0],-d2NNN[3,0]]]);

    bmat_2=np.array([[-d2N[0,1],-d2NN[0,1],-d2NNN[0,1],-d2N[1,1],
                      -d2NN[1,1],-d2NNN[1,1],-d2N[2,1],-d2NN[2,1],
                      -d2NNN[2,1],-d2N[3,1],-d2NN[3,1],-d2NNN[3,1]]]);

    bmat_3=np.array([[-d2N[0,2],-d2NN[0,2],-d2NNN[0,2],-d2N[1,2],
                      -d2NN[1,2],-d2NNN[1,2],-d2N[2,2],-d2NN[2,2],
                      -d2NNN[2,2],-d2N[3,2],-d2NN[3,2],-d2NNN[3,2]]]);

    bmat=np.concatenate((bmat_1,bmat_2,bmat_3));
    
    return bmat;

def B_matrix_plate(a,b,x,y):
    """
    B_matrix is strain-displacement matrix
    
    Currently only handles rectangular elements in the x-y plane. The x and y
    gauss points are the x and y centroid and should be 0,0 for the localized
    element

    Parameters
    ----------
    a : float
        half-length in the x-direction
    b : float
        half-length in the y-direction
    x : float
        x gauss point
    y : float
        y gauss point

    Returns
    -------
    TYPE
        numpy array.

    """
    b1_01=-3*(x-x*y)/(4*a**2);
    b1_02=-((3*a*x-3*a*x*y-a+a*y)/4)/a**2;
    b1_03=0;
    b1_04=-3*(-x+x*y)/(4*a**2);
    b1_05=-((3*a*x-3*a*x*y+a-a*y)/4)/a**2;
    b1_06=0;
    b1_07=-3*(-x-x*y)/(4*a**2);
    b1_08=-((3*a*x+3*a*x*y+a+a*y)/4)/a**2;
    b1_09=0;
    b1_10=-3*(x+x*y)/(4*a**2);
    b1_11=-((3*a*x+3*a*x*y-a-a*y)/4)/a**2;
    b1_12=0;
    b2_01=-3*(y-x*y)/(4*b**2);
    b2_02=0;
    b2_03=-((3*b*y-3*b*x*y-b+b*x)/4)/b**2;
    b2_04=-3*(y+x*y)/(4*b**2);
    b2_05=0;
    b2_06=-((3*b*y+3*b*x*y-b-b*x)/4)/b**2;
    b2_07=-3*(-y-x*y)/(4*b**2);
    b2_08=0;
    b2_09=-((3*b*y+3*b*x*y+b+b*x)/4)/b**2;
    b2_10=-3*(-y+x*y)/(4*b**2);
    b2_11=0;
    b2_12=-((3*b*y-3*b*x*y+b-b*x)/4)/b**2;
    b3_01=-2*(1/2-3*x**2/8-3*y**2/8)/(a*b)
    b3_02=-2*(-3/8*a*x**2+a*x/4+a/8)/(a*b);
    b3_03=-2*(-3/8*b*y**2+b*y/4+b/8)/(a*b);
    b3_04=-2*(-1/2+3*x**2/8+3*y**2/8)/(a*b)
    b3_05=-2*(-3/8*a*x**2-a*x/4+a/8)/(a*b);
    b3_06=-2*(3/8*b*y**2-b*y/4-b/8)/(a*b);
    b3_07=-2*(1/2-3*x**2/8-3*y**2/8)/(a*b);
    b3_08=-2*(3/8*a*x**2+a*x/4-a/8)/(a*b);
    b3_09=-2*(3/8*b*y**2+b*y/4-b/8)/(a*b);
    b3_10=-2*(-1/2+3*x**2/8+3*y**2/8)/(a*b);
    b3_11=-2*(3/8*a*x**2-a*x/4-a/8)/(a*b);
    b3_12=-2*(-3/8*b*y**2-b*y/4+b/8)/(a*b);
    bmat=np.array([
        [b1_01,b1_02,b1_03,b1_04,b1_05,b1_06,b1_07,b1_08,b1_09,b1_10,b1_11,b1_12],
        [b2_01,b2_02,b2_03,b2_04,b2_05,b2_06,b2_07,b2_08,b2_09,b2_10,b2_11,b2_12],
        [b3_01,b3_02,b3_03,b3_04,b3_05,b3_06,b3_07,b3_08,b3_09,b3_10,b3_11,b3_12]
        ]);
    return bmat;

#build the dof dictionary and add the dof to the nodes matrix
dof_db={};
dof=0;
for node in nodes:
    z_t_restraint=not(global_freedoms['z-translation']) or nodes[node]['fix_z'];
    x_r_restraint=not(global_freedoms['x-rotation']) or nodes[node]['fix_mx'];
    y_r_restraint=not(global_freedoms['y-rotation']) or nodes[node]['fix_my'];
    Fz_load=node_loads.get(node,[0.0,0.0,0.0])[0];
    Mx_load=node_loads.get(node,[0.0,0.0,0.0])[1];
    My_load=node_loads.get(node,[0.0,0.0,0.0])[2];
    dof+=1;
    dof_db[dof]={
        'node':node,
        'freedom':'z-translation',
        'restraint':z_t_restraint,
        'load':Fz_load};
    nodes[node]['Uz_dof']=dof;
    dof+=1;
    dof_db[dof]={
        'node':node,
        'freedom':'x-rotation',
        'restraint':x_r_restraint,
        'load':Mx_load};
    nodes[node]['Rx_dof']=dof;
    dof+=1;
    dof_db[dof]={
        'node':node,
        'freedom':'y-rotation',
        'restraint':y_r_restraint,
        'load':My_load};
    nodes[node]['Ry_dof']=dof;
print('Structure consists {} Nodes, {} Elements, {} Degrees of Freedom'.format(
    len(nodes),len(elements),len(dof_db)));
#delete unused variables
del(Fz_load,Mx_load,My_load);


#Create the gauss point list for each node point
gauss_x=[-1/3**0.5,1/3**0.5,1/3**0.5,-1/3**0.5];
gauss_y=[-1/3**0.5,-1/3**0.5,1/3**0.5,1/3**0.5];

#Formulate the plate stiffness matrices
elements_pre={};
print('Pre-processing Elements')
for ele in elements:
    n1=elements[ele]['Node_i'];
    x1,y1=nodes[n1]['x'],nodes[n1]['y'];
    n2=elements[ele]['Node_j'];
    x2,y2=nodes[n2]['x'],nodes[n2]['y'];
    n3=elements[ele]['Node_k'];
    x3,y3=nodes[n3]['x'],nodes[n3]['y'];
    n4=elements[ele]['Node_l'];
    x4,y4=nodes[n4]['x'],nodes[n4]['y'];
    a1=(x2-x1)/2;#Add /2
    a3=(x3-x4)/2;#Add /2
    a=(a1+a3)/2;
    b2=(y3-y2)/2;#Add /2
    b4=(y4-y1)/2;#Add /2
    b=(b2+b4)/2;
    #Initiate the elemental stiffness matrix
    k_local=np.zeros((12,12));
    #calculate the D-matrix
    section=elements[ele]['Section'];
    thk=sections[section]['Thickness'];
    mat=sections[section]['Material'];
    youngs=materials[mat]['E'];
    poisson=materials[mat]['Poisson'];
    dmat_quad=youngs*thk**3/(12*(1-poisson**2))*np.array(
        [[1,poisson,0],
         [poisson,1,0],
         [0,0,0.5*(1-poisson)]]);
    d_mat=youngs*thk**3/(12*(1-poisson**2))*np.array(
        [[1,poisson,0,0,0,0,0,0,0,0,0,0],
         [poisson,1,0,0,0,0,0,0,0,0,0,0],
         [0,0,0.5*(1-poisson),0,0,0,0,0,0,0,0,0],
         [0,0,0,1,poisson,0,0,0,0,0,0,0],
         [0,0,0,poisson,1,0,0,0,0,0,0,0],
         [0,0,0,0,0,0.5*(1-poisson),0,0,0,0,0,0],
         [0,0,0,0,0,0,1,poisson,0,0,0,0],
         [0,0,0,0,0,0,poisson,1,0,0,0,0],
         [0,0,0,0,0,0,0,0,0.5*(1-poisson),0,0,0],
         [0,0,0,0,0,0,0,0,0,1,poisson,0],
         [0,0,0,0,0,0,0,0,0,poisson,1,0],
         [0,0,0,0,0,0,0,0,0,0,0,0.5*(1-poisson)]])
    #Placeholders for the and B-matrix
    b_mat=None;
    for node_pt in range(4):
        x_g=gauss_x[node_pt];
        y_g=gauss_y[node_pt];
        #calculate the b_matrix
        bmat_b=B_matrix_plate(a,b,x_g,y_g);#Using mzc plate
        if(node_pt==0):
            b_mat=bmat_b;
        else:
            b_mat=np.concatenate((b_mat,bmat_b));
        #The @ symbol is used for matrix multiplication
        #k_local+=np.transpose(bmat_b)@dmat_quad@bmat_b*a*b;
    #The @ symbol is used for matrix multiplication
    k_local=np.transpose(b_mat)@d_mat@b_mat*a*b;#Change the area
    elements_pre[ele]={
        'x1':x1,'x2':x2,'x3':x3,'x4':x4,
        'y1':y1,'y2':y2,'y3':y3,'y4':y4,
        'a_edge1':a1,'a_edge3':a3,'a':a,
        'b_edge2':b2,'b_edge4':b4,'b':b,
        'd_mat_quad':dmat_quad,
        'd_mat':d_mat,
        'bmat_b':b_mat,
        'k_local':k_local};    
#delete unused variables
del(a,a1,a3,b,b2,b4,b_mat,bmat_b,d_mat,dmat_quad,dof,ele,k_local,mat,n1,n2,n3,
    n4,node,node_pt,poisson,section,thk,x1,x2,x3,x4,x_g,x_r_restraint,y1,y2,
    y3,y4,y_g,y_r_restraint,youngs,z_t_restraint);

#print the size of elements pre-processing
print('Elements pre-processing storage size = {} kB'.format(
    round(sys.getsizeof(elements_pre)/1024,1)));

#Initialize the global stiffness matrix
k_global=np.zeros((len(dof_db),len(dof_db)));
print('Building the global stiffness matrix');
for ele in elements:
    node_i=elements[ele]['Node_i'];
    dof_1=nodes[node_i]['Uz_dof'];
    dof_2=nodes[node_i]['Rx_dof'];
    dof_3=nodes[node_i]['Ry_dof'];
    node_j=elements[ele]['Node_j'];
    dof_4=nodes[node_j]['Uz_dof'];
    dof_5=nodes[node_j]['Rx_dof'];
    dof_6=nodes[node_j]['Ry_dof'];
    node_k=elements[ele]['Node_k'];
    dof_7=nodes[node_k]['Uz_dof'];
    dof_8=nodes[node_k]['Rx_dof'];
    dof_9=nodes[node_k]['Ry_dof'];
    node_l=elements[ele]['Node_l'];
    dof_10=nodes[node_l]['Uz_dof'];
    dof_11=nodes[node_l]['Rx_dof'];
    dof_12=nodes[node_l]['Ry_dof'];
    dof=[dof_1,dof_2,dof_3,dof_4,dof_5,dof_6,
         dof_7,dof_8,dof_9,dof_10,dof_11,dof_12];
    for i in range(len(dof)):
        for j in range(len(dof)):
            k_global[dof[i]-1][dof[j]-1]+=elements_pre[ele]['k_local'][i][j];

#Delete unused variables
del(ele,node_i,node_j,node_k,node_l,i,j,dof_1,dof_2,dof_3,dof_4,dof_5,dof_6,
    dof_7,dof_8,dof_9,dof_10,dof_11,dof_12,dof,
    );
#print the size of Global stiffness matrix
print('Global Stiffness Matrix size = {} kB'.format(
    round(sys.getsizeof(k_global)/1024,1)));

#Create the a list of unrestrained degrees of freedom
dof_free=[];
for node in nodes:
    if(nodes[node]['fix_z']==False):
        dof_free+=[nodes[node]['Uz_dof']];
    if(nodes[node]['fix_mx']==False):
        dof_free+=[nodes[node]['Rx_dof']];
    if(nodes[node]['fix_my']==False):
        dof_free+=[nodes[node]['Ry_dof']];
del(node);

#Initialize and build the reduced global stiffness matrix
k_global_red=np.zeros((len(dof_free),len(dof_free)));
for i in range(len(dof_free)):
    for j in range(len(dof_free)):
        k_global_red[i][j]+=k_global[dof_free[i]-1][dof_free[j]-1];
del(i,j);
#print the size of the reduced global stiffness matrix
print('Reduced Global Stiffness Matrix size = {} kB'.format(
    round(sys.getsizeof(k_global_red)/1024,1)));

#Initialize the global loads matrix
p_global=np.zeros((len(dof_db),1));
for node in node_loads:
    dof_1=nodes[node]['Uz_dof'];
    dof_2=nodes[node]['Rx_dof'];
    dof_3=nodes[node]['Ry_dof'];
    if(node_loads[node][0]!=0):
        p_global[dof_1-1]+=node_loads[node][0];
    if(node_loads[node][1]!=0):
        p_global[dof_2-1]+=node_loads[node][1];
    if(node_loads[node][2]!=0):
        p_global[dof_3-1]+=node_loads[node][2];
del(node,dof_1,dof_2,dof_3);
#print the size of Global loads matrix
print('Global load Matrix size = {} kB'.format(
    round(sys.getsizeof(p_global)/1024,1)));
#Initialize the reduce global loads matrix
p_global_red=np.zeros((len(dof_free),1));
for i in range(len(dof_free)):
    p_global_red[i][0]+=p_global[dof_free[i]-1][0];
del(i);
print('Reduced Global Loads Matrix size = {} kB'.format(
    round(sys.getsizeof(p_global_red)/1024,1)));

# for i in range(len(nodes_pre['dof_free'])):
#     ii=nodes_pre['dof_free'][i]-1;
#     p_global_red[i]=p_global[ii];
#     for j in range(len(nodes_pre['dof_free'])):
#         jj=nodes_pre['dof_free'][j]-1;
#         k_global_red[i][j]=k_global[ii][jj];
# del(i,j,ii,jj);
# print('Reduced Global Stiffness Matrix size = {} kB'.format(len(k_global_red),
#     round(sys.getsizeof(k_global_red)/1024,1)));

print('Solving Displacements');
time_0=time.time();
# #Calculate the local displacements
u_global_red=np.linalg.inv(k_global_red) @ p_global_red;

#Populate the global displacements
u_global=np.zeros((3*len(nodes),1));
for i in range(len(dof_free)):
    u_global[dof_free[i]-1][0]=u_global_red[i][0];

time_1=time.time();
print('Solved in {} seconds'.format(round(time_1-time_0,3)));
del(i,time_0,time_1)

p_global=k_global@u_global;

#Begin post Processing
print('\nRestrained Nodes Reactions\n');
for node in nodes:
    Rz,Mx,My=0,0,0;
    if(nodes[node]['fix_z'] or nodes[node]['fix_mx'] or nodes[node]['fix_my']):
        print('Node {}'.format(node));
        if(nodes[node]['fix_z']):
            dof_i=3*node-3;
            Rz=round(p_global[dof_i][0]*10**-3,2);
            print(' Rz = {} kN'.format(Rz));
        if(nodes[node]['fix_mx']):
            dof_i=3*node-2;
            Mx=round(p_global[dof_i][0]*10**-6,1);
            print(' Mx = {} kN-m'.format(Mx));
        if(nodes[node]['fix_my']):
            dof_i=3*node-1;
            My=round(p_global[dof_i][0]*10**-6,1);
            print(' My = {} kN-m'.format(My));
        print('');
del(node,dof_i,Rz,Mx,My)

for node in nodes:
    print('Node {} Displacements'.format(node));
    u_z=round(u_global[3*node-3][0],3);
    r_x=round(u_global[3*node-2][0],6);
    r_y=round(u_global[3*node-1][0],6);
    print(' Uz = {} mm'.format(u_z));
    print(' Rx = {} rads'.format(r_x));
    print(' Ry = {} rads'.format(r_y));
    print('');
del(node,u_z,r_x,r_y);

for ele in elements:
    k_local=elements_pre[ele]['k_local'];
    u_local=np.zeros((12,1));
    node_i=elements[ele]['Node_i'];
    node_j=elements[ele]['Node_j'];
    node_k=elements[ele]['Node_k'];
    node_l=elements[ele]['Node_l'];
    dof_inds=[nodes[node_i]['Uz_dof']-1,nodes[node_i]['Rx_dof']-1,nodes[node_i]['Ry_dof']-1,
              nodes[node_j]['Uz_dof']-1,nodes[node_j]['Rx_dof']-1,nodes[node_j]['Ry_dof']-1,
              nodes[node_k]['Uz_dof']-1,nodes[node_k]['Rx_dof']-1,nodes[node_k]['Ry_dof']-1,
              nodes[node_l]['Uz_dof']-1,nodes[node_l]['Rx_dof']-1,nodes[node_l]['Ry_dof']-1];
    for i in range(12):
        u_local[i]=u_global[dof_inds[i]]
    p_local=k_local @ u_local;
    pz_1=round(p_local[0][0]*10**-3,2);
    mx_1,my_1=round(p_local[1][0]*10**-6,1),round(p_local[2][0]*10**-6,1);
    pz_2=round(p_local[3][0]*10**-3,2);
    mx_2,my_2=round(p_local[4][0]*10**-6,1),round(p_local[5][0]*10**-6,1);
    pz_3=round(p_local[6][0]*10**-3,2);
    mx_3,my_3=round(p_local[7][0]*10**-6,1),round(p_local[8][0]*10**-6,1);
    pz_4=round(p_local[9][0]*10**-3,2);
    mx_4,my_4=round(p_local[10][0]*10**-6,1),round(p_local[11][0]*10**-6,1);
    print('Element {} Forces'.format(ele));
    print('Node i: Pz = {} kN ; Mx = {} kN-m ; My = {} kN-m ;'.format(
        pz_1,mx_1,my_1));
    print('Node j: Pz = {} kN ; Mx = {} kN-m ; My = {} kN-m ;'.format(
        pz_2,mx_2,my_2));
    print('Node k: Pz = {} kN ; Mx = {} kN-m ; My = {} kN-m ;'.format(
        pz_3,mx_3,my_3));
    print('Node l: Pz = {} kN ; Mx = {} kN-m ; My = {} kN-m ;'.format(
        pz_4,mx_4,my_4));

#Delete the unused variables
del(ele,i,mx_1,mx_2,mx_3,mx_4,my_1,my_2,my_3,my_4,node_i,node_j,node_k,
    node_l,pz_1,pz_2,pz_3,pz_4);

#Post-processing elemental stress
def B_matrix_plate_stress():
    a=3**0.5;
    bmat_stress=np.array([
        [   (a+1)**2,(1-a)*(a+1),   (1-a)**2,(1-a)*(a+1)],
        [(1-a)*(a+1),   (a+1)**2,(1-a)*(a+1),   (1-a)**2],
        [   (1-a)**2,(1-a)*(a+1),   (a+1)**2,(1-a)*(a+1)],
        [(1-a)*(a+1),   (1-a)**2,(1-a)*(a+1),   (a+1)**2]
        ])*1/4;
    return bmat_stress;

#Initiate an elements post-processing dictionary
elements_post={};
#Loop through the elements and store the elemental forces and
#displacements to the dictioary
for ele in elements:
    #Get the elements nodes
    node_i=elements[ele]['Node_i'];
    node_j=elements[ele]['Node_j'];
    node_k=elements[ele]['Node_k'];
    node_l=elements[ele]['Node_l'];
    #Get the elemental DOF
    ele_dof=[nodes[node_i]['Uz_dof'],nodes[node_i]['Rx_dof'],
              nodes[node_i]['Ry_dof'],nodes[node_j]['Uz_dof'],
              nodes[node_j]['Rx_dof'],nodes[node_j]['Ry_dof'],
              nodes[node_k]['Uz_dof'],nodes[node_k]['Rx_dof'],
              nodes[node_k]['Ry_dof'],nodes[node_l]['Uz_dof'],
              nodes[node_l]['Rx_dof'],nodes[node_l]['Ry_dof']];
    #Get the nodal displacements and forces
    u_local=np.zeros((12,1));
    for i in range(len(ele_dof)):
        u_local[i][0]=u_global[ele_dof[i]-1][0];
    p_local=k_local @ u_local;
    #store to dictionary
    elements_post[ele]={
        'p_local':p_local,
        'u_local':u_local
        };
    
    #Initiate the elemental stress matrix
    str_11=np.zeros((4,1));
    str_22=np.zeros((4,1));
    str_12=np.zeros((4,1));
    #Get the element properties D_mat, a and b
    d_mat=elements_pre[ele]['d_mat_quad'];
    a=elements_pre[ele]['a'];
    b=elements_pre[ele]['b'];
    #Loop through the elements Gauss points and calculate plate stresses
    for i in range(4):
        x_g=gauss_x[i];
        y_g=gauss_x[i];
        b_mat=B_matrix_plate(a,b,x_g,y_g)
        #b_mat=elements_pre[ele]['bmat_b'][3*i:3*(i+1)];
        plate_stress=d_mat @ b_mat @ u_local;
        #Add stresses to relevant matrices
        str_11[i]=plate_stress[0][0];
        str_22[i]=plate_stress[1][0];
        str_12[i]=plate_stress[2][0];
    b_mat_stress=B_matrix_plate_stress();
    str_1=b_mat_stress @ str_11;
    str_2=b_mat_stress @ str_22;
    str_3=b_mat_stress @ str_12;

def print_element_displacements(ele_no,elements,elements_post):
    #Get the elements nodes
    node_i=elements[ele_no]['Node_i'];
    node_j=elements[ele_no]['Node_j'];
    node_k=elements[ele_no]['Node_k'];
    node_l=elements[ele_no]['Node_l'];
    #Get the elements nodal displacemetns
    uzi,rxi,ryi=(round(elements_post[ele_no]['u_local'][0][0],3),
                 round(elements_post[ele_no]['u_local'][1][0],6),
                 round(elements_post[ele_no]['u_local'][2][0],6));
    uzj,rxj,ryj=(round(elements_post[ele_no]['u_local'][3][0],3),
                 round(elements_post[ele_no]['u_local'][4][0],6),
                 round(elements_post[ele_no]['u_local'][5][0],6));
    uzk,rxk,ryk=(round(elements_post[ele_no]['u_local'][6][0],3),
                 round(elements_post[ele_no]['u_local'][7][0],6),
                 round(elements_post[ele_no]['u_local'][8][0],6));
    uzl,rxl,ryl=(round(elements_post[ele_no]['u_local'][9][0],3),
                 round(elements_post[ele_no]['u_local'][10][0],6),
                 round(elements_post[ele_no]['u_local'][11][0],6));
    print('Element: {}  Nodes: {}, {}, {}, {}'.format(
        ele_no,node_i,node_j,node_k,node_l));
    print('Node {}:\tUz= {} mm\tRx= {} rads\tRy= {} rads'.format(
        node_i,uzi,rxi,ryi));
    print('Node {}:\tUz= {} mm\tRx= {} rads\tRy= {} rads'.format(
        node_j,uzj,rxj,ryj));
    print('Node {}:\tUz= {} mm\tRx= {} rads\tRy= {} rads'.format(
        node_k,uzk,rxk,ryk));
    print('Node {}:\tUz= {} mm\tRx= {} rads\tRy= {} rads'.format(
        node_l,uzl,rxl,ryl));    
    return None;

def print_element_edge_forces(ele_no,elements,elements_post):
    #Get the elements nodes
    node_i=elements[ele_no]['Node_i'];
    node_j=elements[ele_no]['Node_j'];
    node_k=elements[ele_no]['Node_k'];
    node_l=elements[ele_no]['Node_l'];
    #Get the elements edge forces
    vx_1=round((elements_post[ele_no]['p_local'][0][0]+
                elements_post[ele_no]['p_local'][9][0])*10**-3,2);
    vx_2=round((elements_post[ele_no]['p_local'][3][0]+
                elements_post[ele_no]['p_local'][6][0])*10**-3,2);
    vy_1=round((elements_post[ele_no]['p_local'][0][0]+
                elements_post[ele_no]['p_local'][3][0])*10**-3,2);
    vy_2=round((elements_post[ele_no]['p_local'][6][0]+
                elements_post[ele_no]['p_local'][9][0])*10**-3,2);    
    mx_1=round((elements_post[ele_no]['p_local'][1][0]+
                elements_post[ele_no]['p_local'][10][0])*10**-6,2);
    mx_2=round((elements_post[ele_no]['p_local'][4][0]+
                elements_post[ele_no]['p_local'][7][0])*10**-6,2);
    my_1=round((elements_post[ele_no]['p_local'][2][0]+
                elements_post[ele_no]['p_local'][5][0])*10**-6,2);
    my_2=round((elements_post[ele_no]['p_local'][8][0]+
                elements_post[ele_no]['p_local'][11][0])*10**-6,2);
    print('Element: {}  Nodes: {}, {}, {}, {}'.format(
        ele_no,node_i,node_j,node_k,node_l));
    print('Lft-Edge {}-{}\tVx= {} kN \tMx= {} kN-m'.format(
        node_i,node_l,vx_1,mx_1));
    print('Rgt-Edge {}-{}\tVx= {} kN \tMx= {} kN-m'.format(
        node_j,node_k,vx_2,mx_2));
    print('Btm-Edge {}-{}\tVy= {} kN \tMy= {} kN-m'.format(
        node_i,node_j,vy_1,my_1));
    print('Top-Edge {}-{}\tVy= {} kN \tMy= {} kN-m'.format(
        node_l,node_k,vy_2,my_2));    
    return None;