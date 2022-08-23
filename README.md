# ShellSet - Parallel Neotectonic Modelling

ShellSet is an MPI parallelization of the combination of OrbData, Shells and OrbScore, to which we have added two input options: a model list and a grid search algorithm.  
  
Within this repository you will find four directories: src, containing all Fortran program files; INPUT, location for storage of input files (containing some examples); EXAMPLES, contains two example ShellSet uses (grid search and list input options); Docs, contains the ShellSet user guide. Also present in the repository root is the Makefile used to compile ShellSet. The Makefile contains three build targets: ShellSet, optimal, debug, and two clean options: clean, cleanall. All Makefile options are outlined in the user guide, as well as instructions for altering the Makefile MACROs for the build options.


#### Outline of each file in src:  
DMods.f90          -> Updated double precision routines by Peter Bird  
MOD_Data.f90       -> Routines used by OrbData5  
MOD_Score.f90      -> Routines used by OrbScore2  
MOD_SharedSubs.f90 -> Globally shared variables and data blocks  
MOD_ShellSet.f90   -> Routines used within the ShellSet update  
MOD_Shells.f90     -> Routines used by Shells  
MOD_VarCheck.f90   -> Routines for pre-run variable checks - separated to simplify user personalisation  
OrbData5.f90       -> OrbData5 program unit - mostly original with few modifications  
OrbScore2.f90      -> OrbScore2 program unit - mostly original with few modifications  
SHELLS_v5.0.f90    -> Shells program unit - mostly original with few modifications  
ShellSetMain.f90   -> The main ShellSet program  
lapack.f90         -> Intel interface for Lapack  
mkl_service.f90    -> Intel interface for MKL  
  
  
ShellSet was developed by Jon May (INGV), Peter Bird (UCLA) & Michele Carafa (INGV)
