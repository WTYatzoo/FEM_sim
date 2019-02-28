#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} simulator EXE=${exe} KIND=${kind} SIMULATION=${simulation} FRAME=${frame} DT=${dt} POISSONRATIO=${PoissonRatio} YOUNGMODULUS=${YoungModulus} LINE_SEARCH=${line_search} WEIGHT_LINE_SEARCH=${weight_line_search} CONSTITUTIVE_MODEL=${constitutive_model} GRAVITY=${gravity} DENSITY=${density} OUT_DIR=${out_dir} FINE_COARSEN=${fine_coarsen} OBJECT_NAME=${object_name} INPUT_OBJECT=${input_object} INPUT_MAT=${input_mat} INPUT_CONSTRAINT=${input_constraint} INPUT_MAT_KIND=${input_mat_kind} LEVEL=${level} FORCE_FUNCTION=${force_function}

