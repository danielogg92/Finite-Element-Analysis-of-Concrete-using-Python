# -*- coding: utf-8 -*-
import sys;

#Add parent directory
if('..' not in sys.path):
    sys.path.insert(0,'..')

import Finite_Element_Plate as fep;
import Graphics.Graphics as fep_graphics;
from Generate_Floor_Mesh import read_text_file;

#Test read_text_file
file_name='Slab_Example_01.txt'
(nodes,elements,sections,materials,global_freedoms,load_cases,
 load_combos,criteria,mesh_input,area_loads,node_loads)=read_text_file(file_name);

#Test run analysis
analysis_output=fep.run_analysis(nodes,elements,sections,materials,
                                 node_loads,load_cases);
#Seperate output of analysis output
nodes=analysis_output[0];
elements=analysis_output[1];
k_global=analysis_output[2];
k_global_red=analysis_output[3];
p_global=analysis_output[4];
p_global_red=analysis_output[5];
p_global_post=analysis_output[6];
u_global=analysis_output[7];
u_global_red=analysis_output[8];

#Test print_load_and_reaction_summary
fep.print_load_and_reaction_summary(load_cases,nodes,p_global,p_global_post);

#Test plot_mesh_2d
fig,ax=fep_graphics.plot_mesh_2d(nodes,elements);
