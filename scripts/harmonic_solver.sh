#! /bin/bash
MAKEFILE=${script_dir}/makefile
make -f ${MAKEFILE} harmonic_solver EXE=${exe} KIND=${kind} SCALE=${scale} POISSONRATIO=${PoissonRatio} YOUNGMODULUS=${YoungModulus} OUT_DIR=${out_dir} INPUT_DIR=${input_dir} OBJECT_NAME=${object_name} LEVEL=${level}

