#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} CtoFer EXE=${exe} OBJECT_NAME=${object_name} INPUT_FINE_OBJECT=${input_fine_object} OUT_DIR=${out_dir} INPUT_COARSEN_REST_OBJECT=${input_coarsen_rest_object} INPUT_COARSEN_DEFORM_OBJECT=${input_coarsen_deform_object} INPUT_HARMONIC_DEF=${input_harmonic_def} LEVEL=${level}
