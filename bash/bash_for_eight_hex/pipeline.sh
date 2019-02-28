#! /bin/bash
##########################################################
project_dir=${HOME}/project/mytest/numerical_coarsening_test
script_dir=${project_dir}/scripts
exe_dir=${project_dir}/build/bin
exe=${exe_dir}/numerical_coarsening
base_dir=${project_dir}/data
cd ${base_dir}/input_data
num_obj=`ls -l |grep "d"|wc -l`
echo num_obj: ${num_obj}

mkdir ${base_dir}/output_data/

mass=(0 1 200) #(0 200 ) #(0 1 30 200) #for each object 
magnitude=( 10 150 ) #for each force_function
magnitude_gravity=(0 10 10) #for each object's gravity's magnitude
magnitude_rotation=(0 150 10) #for each object's rotation's magnitude
force_function_array=( gravity rotation )

object_name_array=( "wtyatzoo" "beam_new" "mybeam18") #( "wtyatzoo" "beam_new" "cube_new")


:<<EOF

for((i=1;i<=1;i+=1))
do

done
EOF

    cd ${base_dir}/input_data
    object_name=heart #${object_name_array[${i}]}
    #object_name= #hand #`ls | xargs -n1 | sed -n ${i}'p'`    
    file_vtk_fine=${base_dir}/input_data/${object_name}/${object_name}.sub1.vtk
    file_out_fine=${base_dir}/output_data/${object_name}/${object_name}_material
    y_homo=5000
    y1=50000
    y2=1000
    p=0.49
    N=38 #18
    C=-0.2 #-0.4
    mkdir -p ${file_out_fine}
    source ${script_dir}/mate_factory.sh


#spectrum_analyser
:<<EOF
for((i=1;i<=1;i+=1))
do
    cd ${base_dir}/input_data
    object_name=mybeam18 #`ls | xargs -n1 | sed -n ${i}'p'`

    file_vtk_fine=${base_dir}/input_data/${object_name}/${object_name}.sub1.vtk

    for((j=2;j<=2;j+=1)) #材料 
    do
	mat_in_fine=${base_dir}/output_data/${object_name}/${object_name}_material/${object_name}_${j}.txt
	mat_out_fine=${base_dir}/input_data/${object_name}/${object_name}.sub1.mat
	source ${script_dir}/exchanger.sh

	#####################################################
	object_name=mybeam18
        out_dir=${base_dir}/output_data/${object_name}/${object_name}_${j}
	input_fine_object=${base_dir}/input_data/${object_name}/${object_name}.sub1.vtk
	input_fine_mat=${base_dir}/input_data/${object_name}/${object_name}.sub1.mat
	input_coarsen_object=${base_dir}/input_data/${object_name}/${object_name}.sub0.vtk
	input_coarsen_mat_new=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}.sub0.mat
	level=2
	mkdir -p ${out_dir}
	source ${script_dir}/spectrum_analyser.sh
	######################################################################
    done    
done
EOF


for((i=2;i<=1;i+=1))
do
    cd ${base_dir}/input_data
    object_name=${object_name_array[${i}]} #`ls | xargs -n1 | sed -n ${i}'p'`
    file_vtk_fine=${base_dir}/input_data/${object_name}/${object_name}.sub1.vtk
    for((j=5;j<=5;j+=1)) #材料 
    do
	mat_in_fine=${base_dir}/output_data/${object_name}/${object_name}_material/${object_name}_${j}.txt
	mat_out_fine=${base_dir}/input_data/${object_name}/${object_name}.sub1.mat
	source ${script_dir}/exchanger.sh

	#####################################################
	kind=use_data_from_file
	scale=400
	level=2
	input_dir=${base_dir}/input_data/${object_name}
	out_dir=${base_dir}/output_data/${object_name}/${object_name}_${j}
	mkdir -p ${out_dir}
	#source ${script_dir}/harmonic_solver.sh

	#####################################################

	input_fine_object=${base_dir}/input_data/${object_name}/${object_name}.sub1.vtk
	input_fine_mat=${base_dir}/input_data/${object_name}/${object_name}.sub1.mat
	input_coarsen_object=${base_dir}/input_data/${object_name}/${object_name}.sub0.vtk
	input_harmonic_def=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}_harmonic_def
	level=2
	mkdir -p ${out_dir}
	#source ${script_dir}/coarsener.sh

	#####################################################
	for((k=0;k<1;k++)) #force pattern
	do
	    force_function=${force_function_array[${k}]}		
	    for((z=9;z<10;z+=1)) #force magnitude
	    do
		if [ ${force_function} = "rotation" ] ; then
		    magnitude_here=${magnitude_rotation[${i}]}
		else
		    magnitude_here=${magnitude_gravity[${i}]}
		fi
		gravity=$( expr "${magnitude_here}*(1.0+19.0/9.0*$z)+0.5" | bc -l | cut -d. -f1 )
	#	gravity=$( expr "scale=10; $z*${magnitude[${k}]}" | bc ) #10位小数
		###################################################
		kind=use_data_from_file
		simulation=static  
		frame=2000
		dt=0.01  # belong to object
		line_search=1  #belong to object
		weight_line_search=1e-6 #belong to object
		density=${mass[${i}]}
		####################################################
		constitutive_model=co_rotated_linear_with_stiffness_tensor #belong to object
		fine_coarsen=coarsen
		input_object=${base_dir}/input_data/${object_name}/${object_name}.sub0.vtk	     
		input_constraint=${base_dir}/input_data/${object_name}/${object_name}.sub0.csv
		input_mat_kind=0
		level=2
		
		out_dir=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}_${fine_coarsen}_${constitutive_model}_new_equa8
		input_mat=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}.sub0.1.mat

		mkdir -p ${out_dir}	        
	#	source ${script_dir}/simulator.sh
	
		######################################################
		:<<EOF
input_fine_object=${base_dir}/input_data/${object_name}/${object_name}.sub1.vtk
		input_coarsen_deform_object=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}_${fine_coarsen}_${constitutive_model}_new_equa5/${object_name}_${force_function}_${simulation}_${gravity}_0.vtk
		input_coarsen_rest_object=${base_dir}/input_data/${object_name}/${object_name}.sub0.vtk
		input_harmonic_def=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}_harmonic_def
		level=2
		mkdir -p ${out_dir}
		source ${script_dir}/CtoFer.sh
EOF

		######################################################
		
	        :<<EOF
out_dir=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}_${fine_coarsen}_${constitutive_model}_new_equa5
		input_mat=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}.sub0.mat
		
		mkdir -p ${out_dir}	        
		source ${script_dir}/simulator.sh
EOF

		######################################################

		:<<EOF
input_fine_object=${base_dir}/input_data/${object_name}/${object_name}.sub1.vtk
		input_coarsen_deform_object=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}_${fine_coarsen}_${constitutive_model}_new_equa8/${object_name}_${force_function}_${simulation}_${gravity}_0.vtk
		input_coarsen_rest_object=${base_dir}/input_data/${object_name}/${object_name}.sub0.vtk
		input_harmonic_def=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}_harmonic_def
		level=2
		mkdir -p ${out_dir}
		source ${script_dir}/CtoFer.sh
EOF
		
		######################################################
		
		constitutive_model=co_rotated_linear #belong to object
		fine_coarsen=fine
		out_dir=${base_dir}/output_data/${object_name}/${object_name}_${j}/${object_name}_${fine_coarsen}_${constitutive_model}_ori
		input_object=${base_dir}/input_data/${object_name}/${object_name}.sub1.vtk
		input_mat=${base_dir}/input_data/${object_name}/${object_name}.sub1.mat
		input_constraint=${base_dir}/input_data/${object_name}/${object_name}.sub1.csv
		input_mat_kind=1
		mkdir -p ${out_dir}
		source ${script_dir}/simulator.sh

	    done   
	done
    done    
done

:<<EOF
EOF
