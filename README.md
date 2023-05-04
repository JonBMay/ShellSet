# ShellSet - Parallel Dynamic Neotectonic Modelling

ShellSet was developed by Dr Jon B. May (INGV), Prof. Peter Bird (UCLA) & Dr Michele M. C. Carafa (INGV). It is an MPI parallelization of the combination of OrbData5, Shells_v5.0 and OrbScore2, to which we have added two input options: a model list and a grid search algorithm.

Within this repository you will find four directories:
1. src -> Fortran program files and Python routines
2. INPUT -> Location for storage of input files (contains examples)
3. EXAMPLES -> Three ShellSet examples
4. Docs -> ShellSet user guide

Also present in the repository is a GNU licence file and a Makefile. The Makefile contains three build targets: ShellSet, optimal, debug, and two clean options: clean, cleanall. All Makefile options are outlined in the user guide, as well as instructions for altering the Makefile MACROs for the build options.

Each ShellSet version will be released with a personal DOI and users of ShellSet should cite the program version using its assigned DOI.

#### Contact:
Please contact jonbryan.may@ingv.it with any issues, questions or suggestions.


## Program requirements:
ShellSet has been designed to work in a Linux environment, either as the main or guest OS (e.g., WSL2). The following is a list of ShellSet requirements, each of which is freely available from Intel in their oneAPI toolkits, specifically the Base and HPC toolkits:
1. Intel Fortran compiler
2. Intel Math Kernel Library
3. Message Passing Interface (MPI) environment


## Examples:
There are 3 examples included within the ShellSet package.
1. List input showing 9 models distributed in a 3x3 grid.
2. Grid search over 3 levels (41 distinct models) updating results from the publication available [here](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2007JB005460).
3. Equal to second example with an additional fault slip rate (FSR) data set used to score models.


## Potential issues with Python routines
Both the plotting routine and GUI have only been tested in Python version 3, specifically version 3.8.10.
Some users have reported issues with GUI visualisation when using WSL2. These issues have been fixed by installing VcXsrv, following instructions from [here](https://techcommunity.microsoft.com/t5/modern-work-app-consult-blog/running-wsl-gui-apps-on-windows-10/ba-p/1493242).
Some users may need to manually install tkinter (as well as other Python libraries), the terminal error message will generally give information on how to do this if required but note that users should download the version which corresponds to their Python version (while also noting that the provided Python routines are only tested in v3.x).


## Outline of each file in src:
DMods.f90 -> Updated double precision routines by Peter Bird  
MOD_Data.f90 -> Routines used by OrbData5  
MOD_Score.f90 -> Routines used by OrbScore2  
MOD_SharedSubs.f90 -> Globally shared variables and data blocks  
MOD_ShellSet.f90 -> Routines used within the ShellSet update  
MOD_Shells.f90 -> Routines used by SHELLS_v5.0 
MOD_VarCheck.f90 -> Routines for pre-run variable checks - separated to simplify user personalisation  
OrbData5.f90 -> OrbData5 program unit - mostly original with few modifications  
OrbScore2.f90 -> OrbScore2 program unit - mostly original with few modifications  
SHELLS_v5.0.f90 -> Shells program unit - mostly original with few modifications  
ShellSetMain.f90 -> The main ShellSet program  
lapack.f90 -> Intel interface for Lapack  
mkl_service.f90 -> Intel interface for MKL  
ShellSetGUI.py -> Graphical User Interface for ShellSet  
ShellSetScatter.py -> Scatter plotter for ShellSet output (1D, 2D & 3D)  


## Outline INPUT files:
age_1p5.grd -> Seafloor ages  
aggregated_offset_rates.dig -> Long-term fault heave and/or throw rates  
Baum887.dig -> [Baumgardner, J.R.](https://digital.library.unt.edu/ark:/67531/metadc1108970/) Figure 7A-F  
CRUST2.grd -> Crustal thickness  
delta_ts.grd -> Travel-time anomaly (s) for vertical S-waves traveling through upper mantle  
Earth5R.feg -> Finite element grid file  
Earth5R-type4A.bcs -> Boundary conditions (no plate-interior VBCs)  
Earth5R-type4AplusA.bcs -> Boundary conditions  
ETOPO20.grd -> 20-minute gridded global relief data  
Fouch_2004_SKS_splitting-selected.dat -> Upper-mantle anisotropy data as fast-polarization azimuths from SKS splitting with SKS splitting times in s  
GCMT_shallow_m5p7_1977-2017.eqc -> Seismic catalog  
GPS2006_selected_subset.gps -> Geodetic benchmark positions and relative horizontal velocities  
GridInput.in -> Information related to the grid search option  
HOC79ii.dig -> [Hager and O'Connell (1979)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JB084iB03p01031?casa_token=Oc-Qr2482YoAAAAA:9VsVcCBfkEIdokY1WOJURg2VK2BqSRHg-HbmRKCWDxzeW94KnzWajcN-jMI1Yps6H_-aL9QbPAzqFw) Model II.   
iEarth5-049.in -> Parameter input file  
InputFiles.in -> Names of input files required for ShellSet  
ListInput.in -> List of models to be run  
magnetic_PB2002.dat -> Seafloor spreading (full) rates at mid-ocean ridges  
PB2002_boundaries.dig -> Digitized plate boundaries, with identification of the polarity of any subduction zones  
PB2002_plates.dig -> Digitized plate outlines  
robust_interpolated_stress_for_OrbScore2.dat -> Most-compressive horizontal principal stress azimuths  
UpVar.in -> List of variables to be updated between the 1st & 2nd iterations of Shells  
