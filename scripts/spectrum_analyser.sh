#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} spectrum_analyser EXE=${exe} OBJECT_NAME=${object_name} INPUT_FINE_OBJECT=${input_fine_object} INPUT_FINE_MAT=${input_fine_mat} OUT_DIR=${out_dir} INPUT_COARSEN_OBJECT=${input_coarsen_object} INPUT_COARSEN_MAT_NEW=${input_coarsen_mat_new} LEVEL=${level}
