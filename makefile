include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
CFLAGS += -pedantic -std=c99

expcircle: expcircle.o
	-${CLINKER} -o expcircle expcircle.o ${PETSC_LIB}
	${RM} expcircle.o

ecjac: ecjac.o
	-${CLINKER} -o ecjac ecjac.o ${PETSC_LIB}
	${RM} ecjac.o

reaction: reaction.o
	-${CLINKER} -o reaction reaction.o ${PETSC_LIB}
	${RM} reaction.o
robin_constraints: robin_constraints.o
	-${CLINKER} -o robin_constraints robin_constraints.o ${PETSC_LIB}
	${RM} robin_constraints.o
copy_robin: copy_robin.o
	-${CLINKER} -o copy_robin copy_robin.o ${PETSC_LIB}
	${RM} copy_robin.o

robin_constraints_N2: robin_constraints_N2.o
	-${CLINKER} -o robin_constraints_N2 robin_constraints_N2.o ${PETSC_LIB}
	${RM} robin_constraints_N2.o
puncture: puncture.o
	-${CLINKER} -o puncture puncture.o ${PETSC_LIB}
	${RM} puncture.o
other_robin_constraints: other_robin_constraints.o
	-${CLINKER} -o other_robin_constraints other_robin_constraints.o ${PETSC_LIB}
	${RM} other_robin_constraints.o




neumann_constraints_kerr_Psi: neumann_constraints_kerr_Psi.o
	-${CLINKER} -o neumann_constraints_kerr_Psi neumann_constraints_kerr_Psi.o ${PETSC_LIB}
	${RM} neumann_constraints_kerr_Psi.o

parameters: parameters.o
	-${CLINKER} -o parameters parameters.o ${PETSC_LIB}
	${RM} parameters.o
kerr: kerr.o
	-${CLINKER} -o kerr kerr.o ${PETSC_LIB}
	${RM} kerr.o

2Kerr: 2Kerr.o
	-${CLINKER} -o 2Kerr 2Kerr.o ${PETSC_LIB}
	${RM} 2Kerr.o

3Kerr: 3Kerr.o
	-${CLINKER} -o 3Kerr 3Kerr.o ${PETSC_LIB}
	${RM} 3Kerr.o





# testing
runexpcircle_1:
	-@../testit.sh expcircle "-snes_fd" 1 1

runexpcircle_2:
	-@../testit.sh expcircle "-snes_mf" 1 2

runecjac_1:
	-@../testit.sh ecjac "" 1 1

runecjac_2:
	-@../testit.sh ecjac "-snes_mf_operator" 1 2

runreaction_1:
	-@../testit.sh reaction "-snes_converged_reason -da_refine 1" 1 1

runreaction_2:
	-@../testit.sh reaction "-da_grid_x 7 -mat_is_symmetric" 1 2

runreaction_3:
	-@../testit.sh reaction "-snes_fd_color -da_refine 1 -snes_view" 4 3

runreaction_4:
	-@../testit.sh reaction "-snes_converged_reason -snes_linesearch_type cp" 1 4

runreaction_5:
	-@../testit.sh reaction "-snes_test_jacobian" 1 5

runreaction_6:
	-@../testit.sh reaction "-snes_converged_reason -rct_noRinJ" 1 6

test_expcircle: runexpcircle_1 runexpcircle_2

test_ecjac: runecjac_1 runecjac_2

test_reaction: runreaction_1 runreaction_2 runreaction_3 runreaction_4 runreaction_5 runreaction_6

test: test_expcircle test_ecjac test_reaction

# etc

.PHONY: distclean runexpcircle_1 runexpcircle_2 runecjac_1 runecjac_2 runreaction_1 runreaction_2 runreaction_3 runreaction_4 runreaction_5 runreaction_6 test test_expcircle test_ecjac test_reaction

distclean:
	@rm -f *~ expcircle ecjac reaction *tmp

