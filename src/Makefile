EXE = heat_eq.x
FORTRAN = nagfor

NC_DIR = /usr/local/netcdf-4.6.1
NC_INC = -I${NC_DIR}/include
NC_LIB = -L${NC_DIR}/lib -lnetcdff -lnetcdf

run: ${EXE}
	./${EXE}

${EXE}: fd1d_heat_explicit.o libfmw.a
	${FORTRAN} $^ -L. -lfmw ${NC_LIB} -o $@

fd1d_heat_explicit.o : fd1d_heat_explicit.f90 libfmw.a 
	${FORTRAN} -c ${NC_INC} -I. $< -o $@

libfmw.a: CFL_mod.o IO_mod.o RHS_mod.o Solver_mod.o Types_mod.o
	ar rcs $@ $^ 

CFL_mod.o : CFL_mod.f90 Types_mod.o 
	${FORTRAN} -c -I. $< -o $@
IO_mod.o : IO_mod.f90 Types_mod.o 
	${FORTRAN} -c ${NC_INC} -I. $< -o $@
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
