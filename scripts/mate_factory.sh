#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} mate_factory EXE=${exe} FILE_VTK_FINE=${file_vtk_fine} FILE_OUT_FINE=${file_out_fine} Y_HOMO=${y_homo} Y1=${y1} Y2=${y2} P=${p} OBJECT_NAME=${object_name} N=${N} C=${C}
