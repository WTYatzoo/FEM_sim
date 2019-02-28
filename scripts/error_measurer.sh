#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} error_measurer EXE=${exe} FINE_DEFORM_MESH=${fine_deform_mesh} COARSEN_DEFORM_MESH=${coarsen_deform_mesh} MAX_ERROR_TXT=${max_error_txt} AVERAGE_ERROR_TXT=${average_error_txt}
