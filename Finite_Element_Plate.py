# -*- coding: utf-8 -*-
import numpy as np;
import time;

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

def run_analysis(nodes,elements,sections,materials,node_loads,load_cases):
    """
    Parameters
    ----------
    nodes : TYPE
        DESCRIPTION.
    elements : TYPE
        DESCRIPTION.
    sections : TYPE
        DESCRIPTION.
    materials : TYPE
        DESCRIPTION.
    node_loads : TYPE
        DESCRIPTION.
    load_cases : TYPE
        DESCRIPTION.

    Returns
    -------
    nodes : TYPE
        DESCRIPTION.
    elements : TYPE
        DESCRIPTION.
    k_global : TYPE
        DESCRIPTION.
    k_global_red : TYPE
        DESCRIPTION.
    p_global : TYPE
        DESCRIPTION.
    p_global_red : TYPE
        DESCRIPTION.
    p_global_post : TYPE
        DESCRIPTION.
    u_global : TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    """
    nodes=preprocess_nodes(nodes);
    elements=preprocess_elements(elements,nodes,sections,materials);
    k_global,k_global_red=preprocess_global_stiffness(nodes,elements);
    p_global,p_global_red=preprocess_global_loads(node_loads,nodes,load_cases);
    p_global_post,u_global,u_global_red=analyze(nodes,k_global_red,k_global,
                                                p_global_red);
    return (nodes,elements,
            k_global,k_global_red,
            p_global,p_global_red,p_global_post,
            u_global,u_global_red);

def preprocess_nodes(nodes):
    """
    Returns nodes
    """
    dof=0;
    for node in nodes:
        dof+=1;
        nodes[node]['Uz_dof']=dof;
        dof+=1;
        nodes[node]['Rx_dof']=dof;
        dof+=1;
        nodes[node]['Ry_dof']=dof;
    return nodes;

def preprocess_elements(elements,nodes,sections,materials):
    """
    Returns elements
    
    """
    #Create the gauss point list for each node point
    gauss_x=[-1/3**0.5,1/3**0.5,1/3**0.5,-1/3**0.5];
    gauss_y=[-1/3**0.5,-1/3**0.5,1/3**0.5,1/3**0.5];
    #Formulate the plate stiffness matrices
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
        elements[ele].update([('x1',x1),('x2',x2),('x3',x3),('x4',x4),
                              ('y1',y1),('y2',y2),('y3',y3),('y4',y4),
                              ('a_edge1',a1),('a_edge3',a3),('a',a),
                              ('b_edge2',b2),('b_edge4',b4),('b',b),
                              ('d_mat_quad',dmat_quad),('d_mat',d_mat),
                              ('bmat_b',b_mat),('k_local',k_local)]);
    return elements;
    
def preprocess_global_stiffness(nodes,elements):
    """
    Returns k_gloabl,k_global_red
    """
    k_global=np.zeros((3*len(nodes),3*len(nodes)));
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
                k_global[dof[i]-1][dof[j]-1]+=elements[ele]['k_local'][i][j];
    #Create the a list of unrestrained degrees of freedom
    dof_free=[];
    for node in nodes:
        if(nodes[node]['fix_z']==False):
            dof_free+=[nodes[node]['Uz_dof']];
        if(nodes[node]['fix_mx']==False):
            dof_free+=[nodes[node]['Rx_dof']];
        if(nodes[node]['fix_my']==False):
            dof_free+=[nodes[node]['Ry_dof']];
    #Initialize and build the reduced global stiffness matrix
    k_global_red=np.zeros((len(dof_free),len(dof_free)));
    for i in range(len(dof_free)):
        for j in range(len(dof_free)):
            k_global_red[i][j]+=k_global[dof_free[i]-1][dof_free[j]-1];
    return k_global,k_global_red;

def preprocess_global_loads(node_loads,nodes,load_cases):
    """
    Returns 
    """
    #Initialize the global loads matrix
    p_global={};
    for case in load_cases:
        p_global[case]=np.zeros((3*len(nodes),1));
    for node in node_loads:
        dof_1=nodes[node]['Uz_dof'];
        dof_2=nodes[node]['Rx_dof'];
        dof_3=nodes[node]['Ry_dof'];
        for case in node_loads[node]:
            fz,mx,my=node_loads[node][case];
            p_global[case][dof_1-1][0]=round(p_global[case][dof_1-1][0]+fz*10**3,1);
            p_global[case][dof_2-1][0]=round(p_global[case][dof_2-1][0]+mx*10**6,1);
            p_global[case][dof_3-1][0]=round(p_global[case][dof_3-1][0]+my*10**6,1);
    #Create the a list of unrestrained degrees of freedom
    dof_free=[];
    for node in nodes:
        if(nodes[node]['fix_z']==False):
            dof_free+=[nodes[node]['Uz_dof']];
        if(nodes[node]['fix_mx']==False):
            dof_free+=[nodes[node]['Rx_dof']];
        if(nodes[node]['fix_my']==False):
            dof_free+=[nodes[node]['Ry_dof']];        
    #Initialize the reduce global loads matrix
    p_global_red={};
    for case in load_cases:
        p_global_red[case]=np.zeros((len(dof_free),1));
        for i in range(len(dof_free)):
            p_global_red[case][i][0]+=p_global[case][dof_free[i]-1][0];
    return p_global,p_global_red;

def analyze(nodes,k_global_red,k_global,p_global_red):
    """
    Returns
    
    """
    
    print('Beginning Analysis');
    time_0=time.time();
    k_global_red_inv=np.linalg.inv(k_global_red);
    u_global_red={};
    for case in p_global_red:
        u_global_red[case]=k_global_red_inv @ p_global_red[case];
    #Create the a list of unrestrained degrees of freedom
    dof_free=[];
    for node in nodes:
        if(nodes[node]['fix_z']==False):
            dof_free+=[nodes[node]['Uz_dof']];
        if(nodes[node]['fix_mx']==False):
            dof_free+=[nodes[node]['Rx_dof']];
        if(nodes[node]['fix_my']==False):
            dof_free+=[nodes[node]['Ry_dof']];  
    #Populate the global displacements and loads
    u_global={};
    p_global_post={};
    for case in p_global_red:
        u_global[case]=np.zeros((3*len(nodes),1));
        for i in range(len(dof_free)):
            u_global[case][dof_free[i]-1][0]=u_global_red[case][i][0];
        p_global_post[case]=k_global @ u_global[case];
    time_1=time.time();
    print('Complete {} load case in {} secs'.format(
        len(p_global_red),round(time_1-time_0,3)));
    return p_global_post,u_global,u_global_red;

def get_reactions(nodes,p_global,p_global_post):
    """
    Returns node_reactions
    """
    node_reactions={};
    for node in nodes:
        if(nodes[node]['fix_z']==True or
           nodes[node]['fix_mx']==True or
           nodes[node]['fix_my']==True):
            node_reactions[node]={};
            for case in p_global:
                node_reactions[node][case]={};
            if(nodes[node]['fix_z']==True):
                dof=nodes[node]['Uz_dof'];
                for case in p_global:
                    p_applied=p_global[case][dof-1][0];
                    p_react=p_global_post[case][dof-1][0];
                    p_react=round((p_react-p_applied)*10**-3,3)
                    node_reactions[node][case]={'Rz':p_react};
            if(nodes[node]['fix_mx']==True):
                dof=nodes[node]['Rx_dof'];
                for case in p_global:
                    p_applied=p_global[case][dof-1][0];
                    p_react=p_global_post[case][dof-1][0];
                    p_react=round((p_react-p_applied)*10**-6,3)
                    node_reactions[node][case]={'RMx':p_react};
            if(nodes[node]['fix_my']==True):
                dof=nodes[node]['Ry_dof'];
                for case in p_global:
                    p_applied=p_global[case][dof-1][0];
                    p_react=p_global_post[case][dof-1][0];
                    p_react=round((p_react-p_applied)*10**-6,3)
                    node_reactions[node][case]={'RMy':p_react};  
    return node_reactions;

def get_nodal_displacement(node,nodes,u_global):
    """
    Returns node_disp
    """
    node_disp={};
    dof_uz=nodes[node]['Uz_dof'];
    dof_rx=nodes[node]['Rx_dof'];
    dof_ry=nodes[node]['Ry_dof'];
    for case in u_global:
        uz=round(u_global[case][dof_uz-1][0],3);
        rx=round(u_global[case][dof_rx-1][0],6);
        ry=round(u_global[case][dof_ry-1][0],6);
        node_disp[case]={
            'uz':uz,
            'rx':rx,
            'ry':ry,};
    return node_disp;

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

def get_displacements_combo_case(load_combo,load_combos,nodes,u_global):
    u_combo={};
    dof_id={};
    for node in nodes:
        u_combo[node]={
            'x':nodes[node]['x'],
            'y':nodes[node]['y'],
            'Uz':0,
            'Rx':0,
            'Ry':0,
            'dof_uz':nodes[node]['Uz_dof'],
            'dof_rx':nodes[node]['Rx_dof'],
            'dof_ry':nodes[node]['Ry_dof'],
            };
        dof_id[nodes[node]['Uz_dof']-1]=[node,0];
        dof_id[nodes[node]['Rx_dof']-1]=[node,1];
        dof_id[nodes[node]['Ry_dof']-1]=[node,2];
    for case_combo in load_combos[load_combo]['combinations'][0]:
        case=case_combo[0];
        case_factor=case_combo[1];
        for i in range(len(u_global[case])):
            dof_node,dof=dof_id[i];
            dof=['Uz','Rx','Ry'][dof];
            u_combo[dof_node][dof]+=case_factor*u_global[case][i][0];
    for node in u_combo:
        u_combo[node]['Uz']=round(u_combo[node]['Uz'],3);
        u_combo[node]['Rx']=round(u_combo[node]['Rx'],6);
        u_combo[node]['Ry']=round(u_combo[node]['Ry'],6);
    return u_combo;

def print_load_and_reaction_summary(load_cases,nodes,p_global,p_global_post):
    """
    Parameters
    ----------
    p_global : dict
        the global loads.

    Returns
    -------
    None.

    """
    #Initial the applied loads dictionary
    cases_applied_loads={};
    cases_reactions={};
    for case in load_cases:
        cases_applied_loads[case]={
            'Fz':0.0
            };
        cases_reactions[case]={
            'Rz':0.0
            };
    #Loop thru p_global and add to cases applied
    for case in p_global:
        for i in range(len(p_global[case])):
            if(i%3==0):
                cases_applied_loads[case]['Fz']+=p_global[case][i][0];
    #Get the node reactions
    node_reactions=get_reactions(nodes,p_global,p_global_post);
    for node in node_reactions:
        for case in node_reactions[node]:
            for reaction in node_reactions[node][case]:
                if(reaction=='Rz'):
                    cases_reactions[case]['Rz']+=node_reactions[node][case][reaction];
    #Print the cases
    print('\nLOAD CASE SUMMARY\n')
    for case in load_cases:
        print('Load Case : {}\nCase No : {}\nCase Type : {}'.format(
            case,load_cases[case]['number'],
            load_cases[case]['case type']));
        print('Applied Load : {} kN    Sum Reactions : {} kN\n'.format(
            round(cases_applied_loads[case]['Fz']*10**-3,3),
            round(cases_reactions[case]['Rz'],3)));
    return None;