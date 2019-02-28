$(message "GET INTO makefile")
simulator:
	$(EXE) prog=simulator kind=${KIND} simulation=${SIMULATION} frame=${FRAME} dt=${DT} poissonratio=${POISSONRATIO} youngmodulus=${YOUNGMODULUS} line_search=${LINE_SEARCH} weight_line_search=${WEIGHT_LINE_SEARCH} constitutive_model=${CONSTITUTIVE_MODEL} gravity=${GRAVITY} density=${DENSITY} out_dir=${OUT_DIR} fine_coarsen=${FINE_COARSEN} object_name=${OBJECT_NAME} input_object=${INPUT_OBJECT} input_mat=${INPUT_MAT} input_constraint=${INPUT_CONSTRAINT} input_mat_kind=${INPUT_MAT_KIND} level=${LEVEL} force_function=${FORCE_FUNCTION}
harmonic_solver:
	$(EXE) prog=harmonic_solver kind=${KIND} scale=${SCALE} poissonratio=${POISSONRATIO} youngmodulus=${YOUNGMODULUS} out_dir=${OUT_DIR} input_dir=${INPUT_DIR} object_name=${OBJECT_NAME} level=${LEVEL}
object_creator:
	$(EXE) prog=object_creator l_size_fine=${L_SIZE_FINE} w_size_fine=${W_SIZE_FINE} h_size_fine=${H_SIZE_FINE} dmetric_fine=${DMETRIC_FINE} object_name=${OBJECT_NAME} out_dir=${OUT_DIR} poissonratio=${POISSONRATIO} youngmodulus1=${YOUNGMODULUS1} youngmodulus2=${YOUNGMODULUS2} level=${LEVEL}
coarsener:
	$(EXE) prog=coarsener object_name=${OBJECT_NAME} input_fine_object=${INPUT_FINE_OBJECT} input_fine_mat=${INPUT_FINE_MAT} out_dir=${OUT_DIR} input_coarsen_object=${INPUT_COARSEN_OBJECT} input_harmonic_def=${INPUT_HARMONIC_DEF} level=${LEVEL}
exchanger:
	${EXE} prog=exchanger file_vtk_fine=${FILE_VTK_FINE} mat_in_fine=${MAT_IN_FINE} mat_out_fine=${MAT_OUT_FINE}
spectrum_analyser:
	${EXE} prog=spectrum_analyser object_name=${OBJECT_NAME} input_fine_object=${INPUT_FINE_OBJECT} input_fine_mat=${INPUT_FINE_MAT} out_dir=${OUT_DIR} input_coarsen_object=${INPUT_COARSEN_OBJECT} input_coarsen_mat_new=${INPUT_COARSEN_MAT_NEW} level=${LEVEL}
CtoFer:
	$(EXE) prog=CtoFer object_name=${OBJECT_NAME} input_fine_object=${INPUT_FINE_OBJECT} out_dir=${OUT_DIR} input_coarsen_rest_object=${INPUT_COARSEN_REST_OBJECT} input_coarsen_deform_object=${INPUT_COARSEN_DEFORM_OBJECT} input_harmonic_def=${INPUT_HARMONIC_DEF} level=${LEVEL}
mate_maker:
	${EXE} prog=mate_maker file_vtk_fine=${FILE_VTK_FINE} mat_out_fine=${MAT_OUT_FINE} vtk_out_fine=${VTK_OUT_FINE} y1=${Y1} y2=${Y2} p=${P} axis=${AXIS} mat_out_fine_txt=${MAT_OUT_FINE_TXT}
mate_factory:
	${EXE} prog=mate_factory file_vtk_fine=${FILE_VTK_FINE} file_out_fine=${FILE_OUT_FINE} y_homo=${Y_HOMO} y1=${Y1} y2=${Y2} p=${P} object_name=${OBJECT_NAME} N=${N} C=${C}
error_measurer:
	${EXE} prog=error_measurer fine_deform_mesh=${FINE_DEFORM_MESH} coarsen_deform_mesh=${COARSEN_DEFORM_MESH} max_error_txt=${MAX_ERROR_TXT} average_error_txt=${AVERAGE_ERROR_TXT}
