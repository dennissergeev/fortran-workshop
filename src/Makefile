EXE = heat_eq.x
FORTRAN = nagfor

run: ${EXE}
	./${EXE}

${EXE}: fd1d_heat_explicit.o libfmw.a
	${FORTRAN} $^ -L. -lfmw -o $@

fd1d_heat_explicit.o : fd1d_heat_explicit.f90 libfmw.a 
	${FORTRAN} -c -I. $< -o $@

libfmw.a: CFL_mod.o IO_mod.o RHS_mod.o Solver_mod.o Types_mod.o
	ar rcs $@ $^ 

CFL_mod.o : CFL_mod.f90 Types_mod.o 
	${FORTRAN} -c -I. $< -o $@
IO_mod.o : IO_mod.f90 Types_mod.o 
	${FORTRAN} -c -I. $< -o $@
RHS_mod.o : RHS_mod.f90 Types_mod.o 
	${FORTRAN} -c -I. $< -o $@
Solver_mod.o : Solver_mod.f90 RHS_mod.o Types_mod.o 
	${FORTRAN} -c -I. $< -o $@
Types_mod.o : Types_mod.f90 
	${FORTRAN} -c $< -o $@

.PHONY: clean
clean:
	-rm -f *.o
	-rm -f *.mod
	-rm -f ${EXE}
