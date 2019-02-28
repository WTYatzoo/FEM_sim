#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} object_creator EXE=${exe} L_SIZE_FINE=${l_size_fine} W_SIZE_FINE=${w_size_fine} H_SIZE_FINE=${h_size_fine} DMETRIC_FINE=${dmetric_fine} OBJECT_NAME=${object_name} OUT_DIR=${out_dir} POISSONRATIO=${PoissonRatio} YOUNGMODULUS1=${YoungModulus1} YOUNGMODULUS2=${YoungModulus2} LEVEL=${level}

