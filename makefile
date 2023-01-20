FC = gfortran

LAPACK_FILES = /usr/local/src/lapack-3.4.1/liblapack.a /usr/local/src/lapack-3.4.1/librefblas.a

exe: resinst.exe

fonctions.o: fonctions.f90 

	$(FC) -c fonctions.f90 

resinst.o: fonctions.f90 resinst.f90

	$(FC) -c resinst.f90

resinst.exe: fonctions.o resinst.o

	$(FC) -o resinst.exe fonctions.o resinst.o $(LAPACK_FILES)

clean :

	rm *.o *.exe
