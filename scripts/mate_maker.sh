#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} mate_maker EXE=${exe} FILE_VTK_FINE=${file_vtk_fine} MAT_OUT_FINE=${mat_out_fine} VTK_OUT_FINE=${vtk_out_fine} Y1=${y1} Y2=${y2} P=${p} AXIS=${axis} MAT_OUT_FINE_TXT=${mat_out_fine_txt}
