# Fortran 90 resonant-tunneling diodes simulator  
  
Solves the Schrödinger-Poisson model for the prediction of the steady-state and time dependent electron distribution in a resonant-tunneling diode.  
Compares the efficiencies of a direct solver (reference) and a solver based on a projection method (fast solver) on resonant states.  
The projection method is based on an accurate algorithm to compute quantum resonances.  
  
- Makefile:  
Compile using the makefile with the command "make".
In the makefile, the variable "FC" should be set to the fortran compiler name and the "LAPACK_FILES" variable should contain the paths to the lapack and blas libraries   
  
- Input:  
The initial steady state solution is read from input files. The input files are in the compressed folder valeurs.zip.  
The input files must be in a folder named "valeurs" located in the working folder     
The time dependent solver is initialized at the initial steady state solution. At time t=0, a shift is applied to the external potential bias (forcing potential). The solver computes the time evolution induced by this potential shift
  
- Ouput:  
The output files are written in the folder "valeurs" 
  - The steps performed by the simulator are written in an output file called "resultats":
  - The electronic density and potential computed by the Schrödinger-Poisson solvers are written in text files
