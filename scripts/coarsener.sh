#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} coarsener EXE=${exe} OBJECT_NAME=${object_name} INPUT_FINE_OBJECT=${input_fine_object} INPUT_FINE_MAT=${input_fine_mat} OUT_DIR=${out_dir} INPUT_COARSEN_OBJECT=${input_coarsen_object} INPUT_HARMONIC_DEF=${input_harmonic_def} LEVEL=${level}
