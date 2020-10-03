# -*- coding: utf-8 -*-
import numpy as np;
import matplotlib.pyplot as plt;
from matplotlib import cm;
from scipy.interpolate import griddata;

def plot_mesh_3d(nodes,elements):
    
    fig=plt.figure();
    ax_3d=fig.add_subplot(111,projection='3d');
    #get the points from nodes
    x_pts=[];
    y_pts=[];
    for node in nodes:
        x_pts+=[nodes[node]['x']];
        y_pts+=[nodes[node]['y']];
    #Create the plot;
    ax_3d.scatter(x_pts,y_pts,0,
               zdir='z',
               s=2,
               c='Black');
    return fig,ax_3d;

def plot_mesh_2d(nodes,elements):
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
    x_pts,y_pts=[],[];
    for node in nodes:
        x_1,y_1=nodes[node]['x'],nodes[node]['y'];
        x_pts+=[x_1];
        y_pts+=[y_1];
    #Get the x and y range
    x_min,x_max=min(x_pts),max(x_pts);
    x_range=x_max-x_min;
    y_min,y_max=min(y_pts),max(y_pts);
    y_range=y_max-y_min;   
    x_size=2**2;
    y_size=int(x_size*y_range/x_range);
    fig,ax=plt.subplots(figsize=(x_size,y_size)); #Create the plot
    #Loop thru elements dictionary and add frames to plot
    for ele in elements:
        node_i=elements[ele]['Node_i'];
        node_j=elements[ele]['Node_j'];
        node_k=elements[ele]['Node_k'];
        node_l=elements[ele]['Node_l'];
        x_i,y_i=nodes[node_i]['x'],nodes[node_i]['y'];
        x_j,y_j=nodes[node_j]['x'],nodes[node_j]['y'];
        x_k,y_k=nodes[node_k]['x'],nodes[node_k]['y'];
        x_l,y_l=nodes[node_l]['x'],nodes[node_l]['y'];
        ax.fill([x_i,x_j,x_k,x_l,x_i],
                [y_i,y_j,y_k,y_l,y_i],
                edgecolor='black',
                facecolor='blue',
                zorder=0)
    #Loop thru nodes dictionary and add nodes to the plot
    ax.scatter(x_pts,y_pts,
               c='red',
               s=6); 
    #Display the restrained points
    for node in nodes:
        if(nodes[node]['fix_z']==True):
            x=nodes[node]['x'];
            y=nodes[node]['y'];
            ax.scatter(x,y,
               c='red',
               marker='^',
               s=150); 
    #ax.imshow(fig,aspect='auto');
    ax.set_xlabel('x-Axis');
    ax.set_ylabel('y-Axis');
    ax.set_title('Structure');
    #plt.grid(b=True);
    return fig,ax;

def save_plot(fig,filename):
    fig.savefig(filename);
    return None;

def plot_displacement_contours(u_case):
    x,y,z=[],[],[];
    for node in u_case:
        x+=[u_case[node]['x']];
        y+=[u_case[node]['y']];
        z+=[u_case[node]['Uz']];
    x,y,z=np.array(x),np.array(y),np.array(z);
    x_min,x_max=x.min(),x.max();
    y_min,y_max=y.min(),y.max();
    nx,ny=10,10;
    xi=np.linspace(x_min,x_max,nx);
    yi=np.linspace(y_min,y_max,ny);
    X,Y=np.meshgrid(xi,yi);
    Z=griddata((x,y),z,(X,Y),method='nearest');
    fig,ax=plt.subplots(figsize=(8,int(8*(y_max-y_min)/(x_max-x_min))));
    ax.contourf(X,Y,Z,20,cmap=cm.RdYlBu);
    ax.set_xlabel('x-Axis');
    ax.set_ylabel('y-Axis');
    ax.set_title('Uz Displacement');
    plt.grid(b=True);    
    return fig,ax;
