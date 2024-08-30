# Makefile for compiling the Serial Spectral Code Version 3.0
# This uses the Intel Fortran Compiler, with FFTW 
# ---------------------------------------------------------------------
compile =mpiifort -g 
link = -L/usr/local/lib -llapack -lblas 
#---------------------------------------------------------------------
OBJ_Ham =  main.o initiate.o hybrd_solve.o minpack.o fcn.o ham.o sc_eq.o open_files.o write_data.o open_final.o write_final.o mfe.o plaquette.o gap.o conv.o
MOD_FILE_1=zeeman_mod.o
test : 
	echo $(OBJ_FILE) 
# -------------------------------------------------------------
self_consistent.exe:  $(MOD_FILE_1) $(OBJ_Ham)
	$(compile)   $(MOD_FILE_1) $(OBJ_Ham) -o self_consistent.exe $(link)
#
zeeman_mod.o: zeeman_mod.f90
	$(compile) -c zeeman_mod.f90
#
main.o: main.f90 $(MOD_FILE_1) 
	$(compile) -c main.f90 
#
initiate.o: initiate.f90
	$(compile) -c initiate.f90
#
hybrd_solve.o: hybrd_solve.f90
	$(compile) -c hybrd_solve.f90
#
minpack.o: minpack.f90
	$(compile) -c minpack.f90
#
fcn.o: fcn.f90
	$(compile) -c fcn.f90
#
ham.o: ham.f90
	$(compile) -c ham.f90
#
sc_eq.o: sc_eq.f90
	$(compile) -c sc_eq.f90
#
open_files.o: open_files.f90
	$(compile) -c open_files.f90
#
write_data.o: write_data.f90
	$(compile) -c write_data.f90
#
open_final.o: open_final.f90
	$(compile) -c open_final.f90
#
write_final.o: write_final.f90
	$(compile) -c write_final.f90
#
mfe.o: mfe.f90
	$(compile) -c mfe.f90
#
plaquette.o: plaquette.f90
	$(compile) -c plaquette.f90
#
gap.o: gap.f90
	$(compile) -c gap.f90
#
conv.o: conv.f90
	$(compile) -c conv.f90
# ------------------------------------------------------------------------ 
clean:
	rm -f *.o *.mod core *.cmp.f90 *.kmo work.* *.d
total_clean:
	rm -f *.o *.mod core *.exe 
## --------------------------------------------------------------------##  
