ml_agg_MIS.o : ../Obj/ml_agg_MIS.c ../Obj/ml_aggregate.h \
        ../Obj/ml_comm.h ../Obj/ml_comminfoop.h ../Obj/ml_defs.h ../Obj/ml_ggraph.h ../Obj/ml_gridfunc.h ../Obj/ml_lapack.h ../Obj/ml_mat_formats.h ../Obj/ml_memory.h ../Obj/ml_operator.h
	$(CC) -c $(CFLAGS) ../Obj/ml_agg_MIS.c -o $@

ml_agg_coupled.o : ../Obj/ml_agg_coupled.c ../Obj/ml_aggregate.h \
        ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_lapack.h ../Obj/ml_memory.h ../Obj/ml_operator.h
	$(CC) -c $(CFLAGS) ../Obj/ml_agg_coupled.c -o $@

ml_agg_dd.o : ../Obj/ml_agg_dd.c ../Obj/ml_aggregate.h \
        ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_lapack.h ../Obj/ml_memory.h ../Obj/ml_operator.h
	$(CC) -c $(CFLAGS) ../Obj/ml_agg_dd.c -o $@

ml_agg_genP.o : ../Obj/ml_agg_genP.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_op_utils.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/ml_agg_genP.c -o $@

ml_agg_uncoupled.o : ../Obj/ml_agg_uncoupled.c ../Obj/ml_aggregate.h \
        ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_lapack.h ../Obj/ml_memory.h ../Obj/ml_operator.h
	$(CC) -c $(CFLAGS) ../Obj/ml_agg_uncoupled.c -o $@

ml_aggregate.o : ../Obj/ml_aggregate.c ../Obj/ml_aggregate.h \
        ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_lapack.h ../Obj/ml_memory.h ../Obj/ml_operator.h
	$(CC) -c $(CFLAGS) ../Obj/ml_aggregate.c -o $@

ml_check.o : ../Obj/ml_check.c ../Obj/ml_1level.h ../Obj/ml_bdrypts.h \
        ../Obj/ml_check.h ../Obj/ml_comm.h ../Obj/ml_comminfoagx.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_grid.h ../Obj/ml_gridfunc.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_memory.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_check.c -o $@

ml_ggraph.o : ../Obj/ml_ggraph.c ../Obj/ml_comm.h ../Obj/ml_comminfoop.h \
        ../Obj/ml_defs.h ../Obj/ml_ggraph.h ../Obj/ml_gridfunc.h ../Obj/ml_mat_formats.h ../Obj/ml_memory.h ../Obj/ml_operator.h
	$(CC) -c $(CFLAGS) ../Obj/ml_ggraph.c -o $@

ml_comm.o : ../Obj/ml_comm.c ../Obj/ml_comm.h ../Obj/ml_defs.h \
        ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_comm.c -o $@

ml_comminfoagx.o : ../Obj/ml_comminfoagx.c ../Obj/ml_comminfoagx.h \
        ../Obj/ml_defs.h ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_comminfoagx.c -o $@

ml_comminfoop.o : ../Obj/ml_comminfoop.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/ml_comminfoop.c -o $@

ml_exch_row.o : ../Obj/ml_exch_row.c ../Obj/ml_1level.h ../Obj/ml_bdrypts.h \
        ../Obj/ml_comm.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_grid.h ../Obj/ml_gridfunc.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_exch_row.c -o $@

ml_bdrypts.o : ../Obj/ml_bdrypts.c ../Obj/ml_bdrypts.h ../Obj/ml_defs.h \
        ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_bdrypts.c -o $@

ml_elementagx.o : ../Obj/ml_elementagx.c ../Obj/ml_elementagx.h \
        ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_elementagx.c -o $@

ml_get_basis.o : ../Obj/ml_get_basis.c ../Obj/ml_defs.h ../Obj/ml_gridfunc.h \
        ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_get_basis.c -o $@

ml_grid.o : ../Obj/ml_grid.c ../Obj/ml_defs.h ../Obj/ml_grid.h \
        ../Obj/ml_gridfunc.h ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_grid.c -o $@

ml_gridagx.o : ../Obj/ml_gridagx.c ../Obj/ml_defs.h ../Obj/ml_elementagx.h \
        ../Obj/ml_gridagx.h ../Obj/ml_intlist.h ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_gridagx.c -o $@

ml_gridfunc.o : ../Obj/ml_gridfunc.c ../Obj/ml_defs.h ../Obj/ml_gridfunc.h \
        ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_gridfunc.c -o $@

ml_mapper.o : ../Obj/ml_mapper.c ../Obj/ml_defs.h ../Obj/ml_mapper.h \
        ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_mapper.c -o $@

ml_pde.o : ../Obj/ml_pde.c ../Obj/ml_1level.h ../Obj/ml_bdrypts.h \
        ../Obj/ml_comm.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_grid.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_memory.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/ml_pde.c -o $@

ml_setup.o : ../Obj/ml_setup.c ../Obj/ml_1level.h ../Obj/ml_bdrypts.h \
        ../Obj/ml_comm.h ../Obj/ml_comminfoagx.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_memory.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_setup.c -o $@

ml_bicgstabl.o : ../Obj/ml_bicgstabl.c ../Obj/ml_bicgstabl.h \
        ../Obj/ml_krylov.h
	$(CC) -c $(CFLAGS) ../Obj/ml_bicgstabl.c -o $@

ml_cg.o : ../Obj/ml_cg.c ../Obj/ml_cg.h ../Obj/ml_krylov.h
	$(CC) -c $(CFLAGS) ../Obj/ml_cg.c -o $@

ml_gmres.o : ../Obj/ml_gmres.c ../Obj/ml_gmres.h ../Obj/ml_krylov.h
	$(CC) -c $(CFLAGS) ../Obj/ml_gmres.c -o $@

ml_krylov.o : ../Obj/ml_krylov.c ../Obj/ml_comm.h ../Obj/ml_krylov.h \
        ../Obj/ml_operator.h
	$(CC) -c $(CFLAGS) ../Obj/ml_krylov.c -o $@

driver.o : ../Obj/driver.c ../Obj/ml_agg_genP.h ../Obj/ml_aggregate.h \
        ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/driver.c -o $@

ml_seg_precond.o : ../Obj/ml_seg_precond.c ../Obj/ml_1level.h \
        ../Obj/ml_agg_genP.h ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bdrypts.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_rap.h ../Obj/ml_seg_precond.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/ml_seg_precond.c -o $@

ml_struct.o : ../Obj/ml_struct.c ../Obj/ml_1level.h ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bdrypts.h ../Obj/ml_comm.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_grid.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_krylov.h ../Obj/ml_lapack.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_memory.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_struct.c -o $@

mli_solver.o : ../Obj/mli_solver.c ../Obj/ml_1level.h ../Obj/ml_bdrypts.h \
        ../Obj/ml_comm.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_grid.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/mli_solver.c -o $@

ml_mat_formats.o : ../Obj/ml_mat_formats.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/ml_mat_formats.c -o $@

ml_matmat_mult.o : ../Obj/ml_matmat_mult.c ../Obj/ml_1level.h \
        ../Obj/ml_bdrypts.h ../Obj/ml_comm.h ../Obj/ml_comminfoagx.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_grid.h ../Obj/ml_gridfunc.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_matmat_mult.c -o $@

ml_op_utils.o : ../Obj/ml_op_utils.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/ml_op_utils.c -o $@

ml_operator.o : ../Obj/ml_operator.c ../Obj/ml_1level.h ../Obj/ml_bdrypts.h \
        ../Obj/ml_defs.h ../Obj/ml_memory.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_operator.c -o $@

ml_operatoragx.o : ../Obj/ml_operatoragx.c ../Obj/ml_comm.h \
        ../Obj/ml_comminfoagx.h ../Obj/ml_memory.h ../Obj/ml_operatoragx.h ../Obj/ml_struct.h
	$(CC) -c $(CFLAGS) ../Obj/ml_operatoragx.c -o $@

ml_rap.o : ../Obj/ml_rap.c ../Obj/ml_1level.h ../Obj/ml_bdrypts.h \
        ../Obj/ml_comm.h ../Obj/ml_comminfoagx.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_grid.h ../Obj/ml_gridfunc.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_rap.c -o $@

ml_rap_utils.o : ../Obj/ml_rap_utils.c ../Obj/ml_1level.h \
        ../Obj/ml_bdrypts.h ../Obj/ml_comm.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_grid.h ../Obj/ml_gridfunc.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_rap_utils.c -o $@

ml_csolve.o : ../Obj/ml_csolve.c ../Obj/ml_1level.h ../Obj/ml_csolve.h \
        ../Obj/ml_defs.h ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_csolve.c -o $@

ml_smoother.o : ../Obj/ml_smoother.c ../Obj/ml_1level.h ../Obj/ml_aztec_utils.h \
        ../Obj/ml_defs.h ../Obj/ml_include.h ../Obj/ml_lapack.h ../Obj/ml_memory.h ../Obj/ml_smoother.h
	$(CC) -c $(CFLAGS) ../Obj/ml_smoother.c -o $@

ml_solver.o : ../Obj/ml_solver.c ../Obj/ml_agg_genP.h ../Obj/ml_aggregate.h \
        ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/ml_solver.c -o $@

ml_superlu.o : ../Obj/ml_superlu.c ../Obj/ml_1level.h ../Obj/ml_bdrypts.h \
        ../Obj/ml_comm.h ../Obj/ml_comminfoop.h ../Obj/ml_csolve.h ../Obj/ml_defs.h ../Obj/ml_grid.h ../Obj/ml_gridfunc.h ../Obj/ml_krylov.h ../Obj/ml_mapper.h ../Obj/ml_mat_formats.h ../Obj/ml_memory.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_rap.h ../Obj/ml_smoother.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_utils.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_superlu.c -o $@

ml_xxt.o : ../Obj/ml_xxt.c ../Obj/ml_agg_genP.h ../Obj/ml_aggregate.h \
        ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/ml_xxt.c -o $@

ml_aztec_utils.o : ../Obj/ml_aztec_utils.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../Obj/ml_aztec_utils.c -o $@

ml_intlist.o : ../Obj/ml_intlist.c ../Obj/ml_defs.h ../Obj/ml_intlist.h \
        ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_intlist.c -o $@

ml_lapack.o : ../Obj/ml_lapack.c ../Obj/ml_aztec_lapack.h \
        ../Obj/ml_defs.h ../Obj/ml_superlu_lapack.h ../Obj/ml_vendor_lapack.h
	$(CC) -c $(CFLAGS) ../Obj/ml_lapack.c -o $@

ml_memory.o : ../Obj/ml_memory.c ../Obj/ml_comm.h ../Obj/ml_defs.h \
        ../Obj/ml_memory.h
	$(CC) -c $(CFLAGS) ../Obj/ml_memory.c -o $@

ml_rbm.o : ../Obj/ml_rbm.c ../Obj/ml_rbm.h
	$(CC) -c $(CFLAGS) ../Obj/ml_rbm.c -o $@

ml_utils.o : ../Obj/ml_utils.c ../Obj/ml_comm.h ../Obj/ml_defs.h \
        ../Obj/ml_memory.h ../Obj/ml_utils.h
	$(CC) -c $(CFLAGS) ../Obj/ml_utils.c -o $@

ml_vec.o : ../Obj/ml_vec.c ../Obj/ml_comm.h ../Obj/ml_defs.h \
        ../Obj/ml_memory.h ../Obj/ml_vec.h
	$(CC) -c $(CFLAGS) ../Obj/ml_vec.c -o $@

ml_recirc : ../Obj/ml_recirc
	$(CC) -c $(CFLAGS) ../Obj/ml_recirc -o $@

ml_read_elas.o : ../examples/ml_read_elas.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_read_elas.c -o $@

convertSund2AZdatafile.o : ../examples/convertSund2AZdatafile.c
	$(CC) -c $(CFLAGS) ../examples/convertSund2AZdatafile.c -o $@

ml_ex1d.o : ../examples/ml_ex1d.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_ex1d.c -o $@

ml_example1d.o : ../examples/ml_example1d.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example1d.c -o $@

ml_example1dGS.o : ../examples/ml_example1dGS.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example1dGS.c -o $@

ml_example2d.o : ../examples/ml_example2d.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example2d.c -o $@

ml_example3d.o : ../examples/ml_example3d.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_example3d.c -o $@

ml_readex.o : ../examples/ml_readex.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_readex.c -o $@

ml_star2d.o : ../examples/ml_star2d.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_star2d.c -o $@

mlguide.o : ../examples/mlguide.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/mlguide.c -o $@

mlguide_par.o : ../examples/mlguide_par.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/mlguide_par.c -o $@

new_readex.o : ../examples/new_readex.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/new_readex.c -o $@

oldml_readex.o : ../examples/oldml_readex.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/oldml_readex.c -o $@

seg_readex.o : ../examples/seg_readex.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/seg_readex.c -o $@

ml_recirc.o : ../examples/ml_recirc.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_recirc.c -o $@

ml_read_salsa.o : ../examples/ml_read_salsa.c ../Obj/ml_agg_genP.h \
        ../Obj/ml_aggregate.h ../Obj/ml_aztec_utils.h ../Obj/ml_bicgstabl.h ../Obj/ml_cg.h ../Obj/ml_comm.h ../Obj/ml_defs.h ../Obj/ml_elementagx.h ../Obj/ml_gmres.h ../Obj/ml_grid.h ../Obj/ml_gridagx.h ../Obj/ml_gridfunc.h ../Obj/ml_include.h ../Obj/ml_intlist.h ../Obj/ml_krylov.h ../Obj/ml_operator.h ../Obj/ml_operatoragx.h ../Obj/ml_pde.h ../Obj/ml_solver.h ../Obj/ml_struct.h ../Obj/ml_vec.h ../Obj/mli_solver.h
	$(CC) -c $(CFLAGS) ../examples/ml_read_salsa.c -o $@

