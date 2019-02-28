#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} exchanger EXE=${exe} FILE_VTK_FINE=${file_vtk_fine} MAT_IN_FINE=${mat_in_fine} MAT_OUT_FINE=${mat_out_fine}
