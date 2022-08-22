# Makefile for ShellSet

# ----- Macros -----
# Compilers
FC  = ifort
MPI = mpiifort

# Compiler Flags
INCL  = -I/opt/intel/oneapi/mkl/2021.3.0/include/intel64/ilp64 -I"/opt/intel/oneapi/mkl/2021.3.0/include"
DBFLAGS  = -check all -traceback -debug
OPTFLAGS = -xHost -ipo -unroll-aggressive

# Linker Flags
MKLFLAGS = -qmkl=parallel

# Paths for prerequisites
vpath %.f90 src
vpath %.mod lib
vpath %.o lib
vpath %.DB.o lib_DB
vpath %.OPT.o lib_OPT

# Object Files
OBJS    = ShellSetMain.o OrbScore2.o MOD_Score.o SHELLS_v5.0.o MOD_Shells.o OrbData5.o MOD_Data.o \
       DMods.o MOD_VarCheck.o MOD_ShellSet.o MOD_SharedVars.o mkl_service.o lapack.o
	   
OBJSDB  = ShellSetMain.DB.o OrbScore2.DB.o MOD_Score.DB.o SHELLS_v5.0.DB.o MOD_Shells.DB.o OrbData5.DB.o MOD_Data.DB.o \
       DMods.DB.o MOD_VarCheck.DB.o MOD_ShellSet.DB.o MOD_SharedVars.DB.o mkl_service.DB.o lapack.DB.o

OBJSOPT = ShellSetMain.OPT.o OrbScore2.OPT.o MOD_Score.OPT.o SHELLS_v5.0.OPT.o MOD_Shells.OPT.o OrbData5.OPT.o MOD_Data.OPT.o \
       DMods.OPT.o MOD_VarCheck.OPT.o MOD_ShellSet.OPT.o MOD_SharedVars.OPT.o mkl_service.OPT.o lapack.OPT.o

# Executable Files
PROG = ShellSet.exe
DB   = ShellSet_DB.exe
OPT  = ShellSet_OPT.exe


# ----- Targets -----
# Program
ShellSet: $(PROG)
ShellSet: CMPLFLAGS = -I"lib"

# Debug
debug: $(DB)
debug: CMPLFLAGS = -I"lib_DB" $(DBFLAGS)

# Optimised
optimal: $(OPT)
optimal: CMPLFLAGS = -I"lib_OPT" $(OPTFLAGS)

# Cleaning
clean:
	rm -rf src/*.bak INPUT/*.bak iteration* Error/ ModelError* *.mod *.o *.bak
	clear
cleanall:
	rm -rf src/*.bak INPUT/*.bak iteration* Error/ ModelError* *.mod *.o *.bak lib/  lib_DB/  lib_OPT/ FatalError* *.exe
	clear


# ----- Building & Linking -----

# ----- Create Program -----
%.o: %.f90
	@echo "compiling $<"
	@$(FC) -diag-disable 10145 $(INCL) $(CMPLFLAGS) -c $^ ||:
	
$(PROG): $(OBJS)
	@echo "compilation command:"
	@echo $(FC) -diag-disable 10145 $(INCL) $(CMPLFLAGS)
	@[ -d lib ] || mkdir -p lib
	@echo "creating $(PROG):"
	$(MPI) $(MKLFLAGS) $(INCL) -o $@ $^
	@mv *.o *.mod lib 2>/dev/null ||:

# ----- Create Debug -----
%.DB.o: %.f90
	@echo "compiling $<"
	@$(FC) -diag-disable 10145 $(INCL) $(CMPLFLAGS) -c $^ ||:
	@mv $*.o $*.DB.o	||:
	
$(DB): $(OBJSDB)
	@echo "compilation command:"
	@echo $(FC) -diag-disable 10145 $(INCL) $(CMPLFLAGS)
	@[ -d lib_DB ] || mkdir -p lib_DB
	@echo "creating $(DB):"
	$(MPI) -diag-disable 10182 $(MKLFLAGS) $(INCL) $(DBFLAGS) -o $@ $^
	@mv *.DB.o *.mod lib_DB 2>/dev/null ||:

# ----- Create Optimal -----
%.OPT.o: %.f90
	@echo "compiling $<"
	@$(FC) -diag-disable 10145 $(INCL) $(CMPLFLAGS) -c $^ ||:
	@mv $*.o $*.OPT.o ||:
	
$(OPT): $(OBJSOPT)
	@echo "compilation command:"
	@echo $(FC) -diag-disable 10145 $(INCL) $(CMPLFLAGS)
	@[ -d lib_OPT ] || mkdir -p lib_OPT
	@echo "creating $(OPT):"
	$(MPI) $(MKLFLAGS) $(INCL) $(OPTFLAGS) -o $@ $^
	@mv *.OPT.o *.mod lib_OPT 2>/dev/null ||:


# ----- Dependency chains -----
ShellSetMain.o   : ShellSetMain.f90 MOD_SharedVars.o MOD_ShellSet.o mkl_service.o SHELLS_v5.0.o OrbScore2.o OrbData5.o MOD_VarCheck.o
OrbScore2.o      : OrbScore2.f90    MOD_SharedVars.o MOD_ShellSet.o DMods.o lapack.o MOD_Score.o
MOD_Shells.o     : MOD_Shells.f90   MOD_SharedVars.o MOD_ShellSet.o DMods.o lapack.o
MOD_Score.o      : MOD_Score.f90    MOD_SharedVars.o MOD_ShellSet.o DMods.o
SHELLS_v5.0.o    : SHELLS_v5.0.f90  MOD_SharedVars.o MOD_ShellSet.o MOD_Shells.o
OrbData5.o       : OrbData5.f90     MOD_SharedVars.o MOD_ShellSet.o MOD_Data.o
MOD_Data.o       : MOD_Data.f90     MOD_SharedVars.o MOD_ShellSet.o
MOD_VarCheck.o   : MOD_VarCheck.f90 MOD_SharedVars.o MOD_ShellSet.o
DMods.o          : DMods.f90        MOD_SharedVars.o MOD_ShellSet.o
MOD_ShellSet.o   : MOD_ShellSet.f90 MOD_SharedVars.o
MOD_SharedVars.o : MOD_SharedVars.f90
mkl_service.o    : mkl_service.f90
lapack.o         : lapack.f90


ShellSetMain.DB.o   : ShellSetMain.f90 MOD_SharedVars.DB.o MOD_ShellSet.DB.o mkl_service.DB.o SHELLS_v5.0.DB.o OrbScore2.DB.o OrbData5.DB.o MOD_VarCheck.DB.o
OrbScore2.DB.o      : OrbScore2.f90    MOD_SharedVars.DB.o MOD_ShellSet.DB.o DMods.DB.o lapack.DB.o MOD_Score.DB.o
MOD_Shells.DB.o     : MOD_Shells.f90   MOD_SharedVars.DB.o MOD_ShellSet.DB.o DMods.DB.o lapack.DB.o
MOD_Score.DB.o      : MOD_Score.f90    MOD_SharedVars.DB.o MOD_ShellSet.DB.o DMods.DB.o
SHELLS_v5.0.DB.o    : SHELLS_v5.0.f90  MOD_SharedVars.DB.o MOD_ShellSet.DB.o MOD_Shells.DB.o
OrbData5.DB.o       : OrbData5.f90     MOD_SharedVars.DB.o MOD_ShellSet.DB.o MOD_Data.DB.o
MOD_Data.DB.o       : MOD_Data.f90     MOD_SharedVars.DB.o MOD_ShellSet.DB.o
MOD_VarCheck.DB.o   : MOD_VarCheck.f90 MOD_SharedVars.DB.o MOD_ShellSet.DB.o
DMods.DB.o          : DMods.f90        MOD_SharedVars.DB.o MOD_ShellSet.DB.o
MOD_ShellSet.DB.o   : MOD_ShellSet.f90 MOD_SharedVars.DB.o
MOD_SharedVars.DB.o : MOD_SharedVars.f90
mkl_service.DB.o    : mkl_service.f90
lapack.DB.o         : lapack.f90


ShellSetMain.OPT.o   : ShellSetMain.f90 MOD_SharedVars.OPT.o MOD_ShellSet.OPT.o mkl_service.OPT.o SHELLS_v5.0.OPT.o OrbScore2.OPT.o OrbData5.OPT.o MOD_VarCheck.OPT.o
OrbScore2.OPT.o      : OrbScore2.f90    MOD_SharedVars.OPT.o MOD_ShellSet.OPT.o DMods.OPT.o lapack.OPT.o MOD_Score.OPT.o
MOD_Shells.OPT.o     : MOD_Shells.f90   MOD_SharedVars.OPT.o MOD_ShellSet.OPT.o DMods.OPT.o lapack.OPT.o
MOD_Score.OPT.o      : MOD_Score.f90    MOD_SharedVars.OPT.o MOD_ShellSet.OPT.o DMods.OPT.o
SHELLS_v5.0.OPT.o    : SHELLS_v5.0.f90  MOD_SharedVars.OPT.o MOD_ShellSet.OPT.o MOD_Shells.OPT.o
OrbData5.OPT.o       : OrbData5.f90     MOD_SharedVars.OPT.o MOD_ShellSet.OPT.o MOD_Data.OPT.o
MOD_Data.OPT.o       : MOD_Data.f90     MOD_SharedVars.OPT.o MOD_ShellSet.OPT.o
MOD_VarCheck.OPT.o   : MOD_VarCheck.f90 MOD_SharedVars.OPT.o MOD_ShellSet.OPT.o
DMods.OPT.o          : DMods.f90        MOD_SharedVars.OPT.o MOD_ShellSet.OPT.o
MOD_ShellSet.OPT.o   : MOD_ShellSet.f90 MOD_SharedVars.OPT.o
MOD_SharedVars.OPT.o : MOD_SharedVars.f90
mkl_service.OPT.o    : mkl_service.f90
lapack.OPT.o         : lapack.f90
