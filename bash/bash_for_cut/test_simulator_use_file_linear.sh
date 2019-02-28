#! /bin/bash

#parameter for solvers

kind=use_data_from_file

##############simulation###############
# simulation  dynamic or static 
simulation=dynamic 

# if dynamic ,the num of frame
frame=800

# if dynamic ,the dt is needed
dt=0.01  # belong to object

#gravity
gravity=10

#density
density=20 # belong to object
#############common para################
# whether it uses line_search
line_search=0  #belong to object

# if using line_search, the weight
weight_line_search=1e-4  #belong to object

##########################################################
project_dir=../..
script_dir=${project_dir}/scripts
exe_dir=${project_dir}/build/bin
exe=${exe_dir}/numerical_coarsening
base_output_dir=${project_dir}/example
##########################################################
# constitutive model
constitutive_model=linear #belong to object

fine_coarsen=fine

object_name=cut

out_dir=${base_output_dir}/${object_name}/linear/

#input_object decides the fine_coarsen is fine or coarsen 
input_object=${project_dir}/input_data/${object_name}/${object_name}.vtk

#input_mat 1.decides the ori or new and 2.decides the constitutive_model is with stiffness_tensor or not and 3.decides the input_mat_kind is 1 or 0
input_mat=  # here we do not use file input

input_constraint=${project_dir}/input_data/${object_name}/${object_name}.csv

input_mat_kind=1

force_function=gravity

mkdir -p ${out_dir}
source ${script_dir}/simulator.sh
