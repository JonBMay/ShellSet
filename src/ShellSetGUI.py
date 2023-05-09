import tkinter as tk
from tkinter import *
from tkinter import messagebox
from tkinter import ttk
import os

global ParamFile, InFile, InDir, Suppress, Default

# -------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------ User Edit --------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------

ParamFile = 'iEarth5-049.in'     # Name of existing Parameter input file (to update existing values)
InFile    = 'InputFiles.in'      # Name of existing Input file list file (to update existing file names)
# -------------------------------------------------------------------------------------------------------------------------
InDir     = 'INPUT'              # Name of directory where all ShellSet input files are stored
Suppress  =  False               # Prevent warning messages appearing after selecting to update Parameter input file
Default   = 'Standard'           # Set default ShellSet version - Standard, Optimal, Debug

# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------------



# Param file update
def ParamInput():

  global data, ParamFile, Suppress
  data = []


  def PlanetsInfo():
    string = ('This interface to the input-parameter file performs reasonableness checks '
              'on variable values, assuming that they are in SI units, and that they refer '
              'to the Earth.  If you are modeling another planet, different ranges may '
              'be appropriate (e.g., for gravitational acceleration on the surface). '
              'In such cases, you may wish to skip this GUI interface, and build your '
              'input-parameter file directly with any simple ASCII text editor.')

    messagebox.showinfo('Which planet are you modeling?',string)


  def NoParamFile():
    string = ('No Parameter input file "'+ParamFile+'" found in the input directory "'+InDir+'", '
              'therefore these variables will be prefilled with example values. '
              'If you expect the file to exist then check for spelling errors. '
              'If the file name or input directory are incorrect then edit them at the '
              'beginning of the Python file "ShellSetGUI.py" and restart the GUI.')

    messagebox.showinfo('Parameter input file not found',string)


  def UnitsInfo():
    string = ('All files input to OrbData and to Shells must use the SI system of units. '
              'Therefore, use m for distance, s for time, m/s for velocity, kg/m^3 for '
              'density, W/m^2 for heat-flow, and so on. '
              'However, note that the scoring datasets that are read by OrbScore use '
              'conventional units of mm/a for velocities (fault offset rates, seafloor spreading rates, '
              'geodetic velocities) and also use degrees to express azimuths measured '
              'clockwise from North (principal stress directions, seismic anisotropy).')

    messagebox.showinfo('System of Units',string)


  def Save():

    global data, ParamFile

    diff = tk.BooleanVar(Param)
    fail = tk.BooleanVar(Param)
    diff = False
    fail = False

    PlateNames = ['AF','AM','AN','AP','AR','AS','AT','AU','BH','BR','BS','BU','CA','CL','CO','CR','EA',
                  'EU','FT','GP','IN','JF','JZ','KE','MA','MN','MO','MS','NA','NB','ND','NH','NI','NZ',
                  'OK','ON','PA','PM','PS','RI','SA','SB','SC','SL','SO','SS','SU','SW','TI','TO','WL','YA']
# Check input are correct type
    try :
      float(FFric.get())
    except :
      messagebox.showerror('Error','Fault friction coefficient is expected to be a float')
      fail = True
    try :
      float(CFric.get())
    except :
      messagebox.showerror('Error','Continuum friction coefficient is expected to be a float')
      fail = True
    try :
      float(Biot.get())
    except :
      messagebox.showerror('Error','Biot coefficient is expected to be a float')
      fail = True
    try :
      float(Byerly.get())
    except :
      messagebox.showerror('Error','Byerly entry is expected to be a float')
      fail = True
    try :
      float(ACreepC.get())
    except :
      messagebox.showerror('Error','ACreep first entry is expected to be a float')
      fail = True
    try :
      float(ACreepM.get())
    except :
      messagebox.showerror('Error','ACreep second entry is expected to be a float')
      fail = True
    try :
      float(BCreepC.get())
    except :
      messagebox.showerror('Error','BCreep first entry is expected to be a float')
      fail = True
    try :
      float(BCreepM.get())
    except :
      messagebox.showerror('Error','BCreep second entry is expected to be a float')
      fail = True
    try :
      float(CCreepC.get())
    except :
      messagebox.showerror('Error','CCreep first entry is expected to be a float')
      fail = True
    try :
      float(CCreepM.get())
    except :
      messagebox.showerror('Error','CCreep second entry is expected to be a float')
      fail = True
    try :
      float(DCreepC.get())
    except :
      messagebox.showerror('Error','DCreep first entry is expected to be a float')
      fail = True
    try :
      float(DCreepM.get())
    except :
      messagebox.showerror('Error','DCreep second entry is expected to be a float')
      fail = True
    try :
      float(ECreep.get())
    except :
      messagebox.showerror('Error','ECreep is expected to be a float')
      fail = True
    try :
      float(AdiabatC.get())
    except :
      messagebox.showerror('Error','Intercept of upper mantle is expected to be a float')
      fail = True
    try :
      float(AdiabatM.get())
    except :
      messagebox.showerror('Error','Slope of upper mantle is expected to be a float')
      fail = True
    try :
      float(Zbasth.get())
    except :
      messagebox.showerror('Error','Base of Asthenosphere is expected to be a float')
      fail = True
    if AFRef.get() not in PlateNames:
      string = ('Plate reference frame is not in the approved list, choose from (case sensitive):\n'
                'AF, AM, AN, AP, AR, AS, AT, AU, BH, BR, BS, BU, CA, CL, CO, CR, EA,'
                'EU, FT, GP, IN, JF, JZ, KE, MA, MN, MO, MS, NA, NB, ND, NH, NI, NZ,'
                'OK, ON, PA, PM, PS, RI, SA, SB, SC, SL, SO, SS, SU, SW, TI, TO, WL, YA')
      messagebox.showerror('Error',string)
      fail = True
    if not Iconve1.get().isdigit():
      messagebox.showerror('Error','Iconve first entry is expected to be an integer')
      fail = True
    try :
      float(Iconve2.get())
    except :
      messagebox.showerror('Error','Iconve second entry is expected to be a float')
      fail = True
    try :
      float(TRHmax.get())
    except :
      messagebox.showerror('Error','TRHmax is expected to be a float')
      fail = True
    try :
      float(TAUmax.get())
    except :
      messagebox.showerror('Error','Taumax is expected to be a float')
      fail = True
    try :
      float(RhoH2O.get())
    except :
      messagebox.showerror('Error','RhoH2O is expected to be a float')
      fail = True
    try :
      float(RhoBarC.get())
    except :
      messagebox.showerror('Error','RhoBar first entry is expected to be a float')
      fail = True
    try :
      float(RhoBarM.get())
    except :
      messagebox.showerror('Error','RhoBar second entry is expected to be a float')
      fail = True
    try :
      float(RhoAst.get())
    except :
      messagebox.showerror('Error','RhoAst is expected to be a float')
      fail = True
    try :
      float(Grav.get())
    except :
      messagebox.showerror('Error','Gravity is expected to be a float')
      fail = True
    try :
      float(OneKM.get())
    except :
      messagebox.showerror('Error','One KM is expected to be a float')
      fail = True
    try :
      float(Radius.get())
    except :
      messagebox.showerror('Error','Radius is expected to be a float')
      fail = True
    try :
      float(VolExpC.get())
    except :
      messagebox.showerror('Error','Volumetric expansion first entry is expected to be a float')
      fail = True
    try :
      float(VolExpM.get())
    except :
      messagebox.showerror('Error','Volumetric expansion second entry is expected to be a float')
      fail = True
    try :
      float(ThrConC.get())
    except :
      messagebox.showerror('Error','Thermal conductivity first entry is expected to be a float')
      fail = True
    try :
      float(ThrConM.get())
    except :
      messagebox.showerror('Error','Thermal conductivity second entry is expected to be a float')
      fail = True
    try :
      float(RadHeatC.get())
    except :
      messagebox.showerror('Error','Radioactive heat production first entry is expected to be a float')
      fail = True
    try :
      float(RadHeatM.get())
    except :
      messagebox.showerror('Error','Radioactive heat production second entry is expected to be a float')
      fail = True
    try :
      float(Temp.get())
    except :
      messagebox.showerror('Error','Surface temperature is expected to be a float')
      fail = True
    try :
      float(UppManC.get())
    except :
      messagebox.showerror('Error','Upper temperature limits first entry is expected to be a float')
      fail = True
    try :
      float(UppManM.get())
    except :
      messagebox.showerror('Error','Upper temperature limits second entry is expected to be a float')
      fail = True
    if not MaxIter.get().isdigit():
      messagebox.showerror('Error','Maximum iterations is expected to be an integer')
      fail = True
    try :
      float(AccConv.get())
    except :
      messagebox.showerror('Error','Acceptable convergence level is expected to be a float')
      fail = True
    try :
      float(RefStress.get())
    except :
      messagebox.showerror('Error','Reference level of shear stress is expected to be a float')
      fail = True
    try :
      float(AccVel.get())
    except :
      messagebox.showerror('Error','Acceptable level of velocity errors is expected to be a float')
      fail = True
    if Out.get() not in ['T','F']:
      messagebox.showerror('Error','Output node velocities is expected to be either T or F (case sensitive)')
      fail = True

# Parameter value checks
    if not fail:
      if float(FFric.get()) < 0.0  or  float(FFric.get()) > 1.0 :
        messagebox.showerror('Error','Fault friction coefficient should in range [0,1] \n(all inclusive)')
        fail = True
      if float(CFric.get()) < 0.0  or  float(CFric.get()) > 1.0:
        messagebox.showerror('Error','Continuum friction coefficient should in range [0,1] \n(all inclusive)')
        fail = True
      if float(Byerly.get()) < 0.0  or  float(Byerly.get()) >= 1.0:
        messagebox.showerror('Error','Byerly variable should in range [0,1) \n(0 inclusive, 1 not)')
        fail = True

# Check for differences & save
    if not fail:
      if os.path.exists(InDir+'/'+ParamFile):
        with open(InDir+'/'+ParamFile, 'r') as file:
          data = file.readlines()

        if FFric.get() != data[1].split()[0]:
          data[1] = FFric.get()+' \t '+' '.join(data[1].split()[1:])+'\n'
          diff = True
        if CFric.get() != data[2].split()[0]:
          data[2] = CFric.get()+' \t '+' '.join(data[2].split()[1:])+'\n'
          diff = True
        if Biot.get() != data[3].split()[0]:
          data[3] = Biot.get()+' \t '+' '.join(data[3].split()[1:])+'\n'
          diff = True
        if Byerly.get() != data[4].split()[0]:
          data[4] = Byerly.get()+' \t '+' '.join(data[4].split()[1:])+'\n'
          diff = True
        if ACreepC.get() != data[5].split(',')[0] or ACreepM.get() != data[5].split(',')[1].split()[0]:
          data[5] = ACreepC.get()+','+ACreepM.get()+' \t '+' '.join(data[5].split()[1:])+'\n'
          diff = True
        if BCreepC.get() != data[6].split(',')[0] or BCreepM.get() != data[6].split(',')[1].split()[0]:
          data[6] = BCreepC.get()+','+BCreepM.get()+' \t '+' '.join(data[6].split()[1:])+'\n'
          diff = True
        if CCreepC.get() != data[7].split(',')[0] or CCreepM.get() != data[7].split(',')[1].split()[0]:
          data[7] = CCreepC.get()+','+CCreepM.get()+' \t '+' '.join(data[7].split()[1:])+'\n'
          diff = True
        if DCreepC.get() != data[8].split(',')[0] or DCreepM.get() != data[8].split(',')[1].split()[0]:
          data[8] = DCreepC.get()+','+DCreepM.get()+' \t '+' '.join(data[8].split()[1:])+'\n'
          diff = True
        if ECreep.get() != data[9].split()[0]:
          data[9] = ECreep.get()+' \t '+' '.join(data[9].split()[1:])+'\n'
          diff = True
        if AdiabatC.get() != data[10].split(',')[0] or AdiabatM.get() != data[10].split(',')[1].split()[0]:
          data[10] = AdiabatC.get()+','+AdiabatM.get()+' \t '+' '.join(data[10].split()[1:])+'\n'
          diff = True
        if Zbasth.get() != data[11].split()[0]:
          data[11] = Zbasth.get()+' \t '+' '.join(data[11].split()[1:])+'\n'
          diff = True
        if AFRef.get() != data[12].split()[0]:
          data[12] = AFRef.get()+' \t '+' '.join(data[12].split()[1:])+'\n'
          diff = True
        if Iconve1.get() != data[13].split(',')[0] or Iconve2.get() != data[13].split(',')[1].split()[0]:
          data[13] = Iconve1.get()+','+Iconve2.get()+' \t '+' '.join(data[13].split()[1:])+'\n'
          diff = True
        if TRHmax.get() != data[14].split()[0]:
          data[14] = TRHmax.get()+' \t '+' '.join(data[14].split()[1:])+'\n'
          diff = True
        if TAUmax.get() != data[15].split()[0]:
          data[15] = TAUmax.get()+' \t '+' '.join(data[15].split()[1:])+'\n'
          diff = True
        if RhoH2O.get() != data[16].split()[0]:
          data[16] = RhoH2O.get()+' \t '+' '.join(data[16].split()[1:])+'\n'
          diff = True
        if RhoBarC.get() != data[17].split(',')[0] or RhoBarM.get() != data[17].split(',')[1].split()[0]:
          data[17] = RhoBarC.get()+','+RhoBarM.get()+' \t '+' '.join(data[17].split()[1:])+'\n'
          diff = True
        if RhoAst.get() != data[18].split()[0]:
          data[18] = RhoAst.get()+' \t '+' '.join(data[18].split()[1:])+'\n'
          diff = True
        if Grav.get() != data[19].split()[0]:
          data[19] = Grav.get()+' \t '+' '.join(data[19].split()[1:])+'\n'
          diff = True
        if OneKM.get() != data[20].split()[0]:
          data[20] = OneKM.get()+' \t '+' '.join(data[20].split()[1:])+'\n'
          diff = True
        if Radius.get() != data[21].split()[0]:
          data[21] = Radius.get()+' \t '+' '.join(data[21].split()[1:])+'\n'
          diff = True
        if VolExpC.get() != data[22].split(',')[0] or VolExpM.get() != data[22].split(',')[1].split()[0]:
          data[22] = VolExpC.get()+','+VolExpM.get()+' \t '+' '.join(data[22].split()[1:])+'\n'
          diff = True
        if ThrConC.get() != data[23].split(',')[0] or ThrConM.get() != data[23].split(',')[1].split()[0]:
          data[23] = ThrConC.get()+','+ThrConM.get()+' \t '+' '.join(data[23].split()[1:])+'\n'
          diff = True
        if RadHeatC.get() != data[24].split(',')[0] or RadHeatM.get() != data[24].split(',')[1].split()[0]:
          data[24] = RadHeatC.get()+','+RadHeatM.get()+' \t '+' '.join(data[24].split()[1:])+'\n'
          diff = True
        if Temp.get() != data[25].split()[0]:
          data[25] = Temp.get()+' \t '+' '.join(data[25].split()[1:])+'\n'
          diff = True
        if UppManC.get() != data[26].split(',')[0] or UppManM.get() != data[26].split(',')[1].split()[0]:
          data[26] = UppManC.get()+','+UppManM.get()+' \t '+' '.join(data[26].split()[1:])+'\n'
          diff = True
        if MaxIter.get() != data[27].split()[0]:
          data[27] = MaxIter.get()+' \t '+' '.join(data[27].split()[1:])+'\n'
          diff = True
        if AccConv.get() != data[28].split()[0]:
          data[28] = AccConv.get()+' \t '+' '.join(data[28].split()[1:])+'\n'
          diff = True
        if RefStress.get() != data[29].split()[0]:
          data[29] = RefStress.get()+' \t '+' '.join(data[29].split()[1:])+'\n'
          diff = True
        if AccVel.get() != data[30].split()[0]:
          data[30] = AccVel.get()+' \t '+' '.join(data[30].split()[1:])+'\n'
          diff = True
        if Out.get() != data[31].split()[0]:
          data[31] = Out.get()+' \t '+' '.join(data[31].split()[1:])+'\n'
          diff = True

      else:
        diff = True
        data = []

        data.append(ParamFile+':\n')
        data.append(FFric.get()+' \t \t \t'+' fFric = fault friction coefficient (during slip)\n')
        data.append(CFric.get()+' \t \t \t'+' cFric = continuum friction coefficient\n')
        data.append(Biot.get()+' \t \t \t'+' Biot  = coefficient of pore pressure in effective-normal-stress equation\n')
        data.append(Byerly.get()+' \t \t \t'+' Byerly= fractional friction reduction on master fault [Bird & Kong, 1994]\n')
        data.append(ACreepC.get()+','+ACreepM.get()+' \t \t \t'+' aCreep = shear-stress coefficient of creep law, crust/mantle\n')
        data.append(BCreepC.get()+','+BCreepM.get()+' \t \t \t'+' bCreep = creep activation energy/n/R, in K, crust/mantle\n')
        data.append(CCreepC.get()+','+CCreepM.get()+' \t \t \t'+' cCreep = derivative of exponential numerator w.r.t. depth, in K/m, crust/mantle\n')
        data.append(DCreepC.get()+','+DCreepM.get()+' \t \t \t'+' dCreep = maximum shear stress at any temperature/strain-rate, Pa, crust/mantle\n')
        data.append(ECreep.get()+' \t \t \t'+' eCreep = exponent on strain-rate in creep-strength law, = 1/n, same for crust & mantle\n')
        data.append(AdiabatC.get()+','+AdiabatM.get()+' \t \t \t'+' tAdiab, gradie = intercept (in K) and slope (in K/m) of upper mantle adiabat\n')
        data.append(Zbasth.get()+' \t \t \t'+' zBAsth = depth (in m) of base of upper mantle (end of olivine=rich layer)\n')
        data.append(AFRef.get()+' \t \t \t'+' pltRef = plate held fixed in boundary conditions (or reference frame in global model)\n')
        data.append(Iconve1.get()+','+Iconve2.get()+' \t \t \t'+' iConve, vTimes = convection under lithosphere (codes 0:6 given below, vTimes needed for iConve > 0; also see trHMax below)\n')
        data.append(TRHmax.get()+' \t \t \t'+' trHMax = upper limit on basal tractions from mantle convection (in Pa; may be 0.0 for free-slip, regardless of iConve)\n')
        data.append(TAUmax.get()+' \t \t \t'+'  tauMax = limit on down-dip integral of tractions in subduction megathrusts, oceanic/continental (or one value for both)\n')
        data.append(RhoH2O.get()+' \t \t \t'+' rhoH2O = density of (salty?) water, in kg/m**3, at atmospheric pressure and surface T (e.g., 4 C)\n')
        data.append(RhoBarC.get()+','+RhoBarM.get()+' \t \t \t'+' rhoBar = mean densities of crust/mantle rocks, in kg/m**3, at 0K and low P\n')
        data.append(RhoAst.get()+' \t \t \t'+' rhoAst = density of asthenosphere, in kg/m**3, at its typical adiabatic temperature, but low P\n')
        data.append(Grav.get()+' \t \t \t'+' gMean  = gravitational acceleration at surface of planet, m/s**2 (e.g., 9.8 for Earth)\n')
        data.append(OneKM.get()+' \t \t \t'+' oneKm  = length of 1 kilometer, expressed in current length units (e.g., 1000. meters if using SI units)\n')
        data.append(Radius.get()+' \t \t \t'+' radius = mean radius of the planet, in m (if using SI units) (e.g., 6371000. for Earth)\n')
        data.append(VolExpC.get()+','+VolExpM.get()+' \t \t \t'+' alphaT = volumetric thermal expansion coefficients, in /C or /K, crust/mantle\n')
        data.append(ThrConC.get()+','+ThrConM.get()+' \t \t \t'+' conduc = thermal conductivity, crust/mantle (for SI units, in W/m/C)\n')
        data.append(RadHeatC.get()+','+RadHeatM.get()+' \t \t \t'+' radio = volumetric radioactive heat production (for SI units, in W/m**3)\n')
        data.append(Temp.get()+' \t \t \t'+' tSurf = surface temperature of planet, in K\n')
        data.append(UppManC.get()+','+UppManM.get()+' \t \t \t'+' temLim = temperature limits (due to melting) in crust/mantle-lithosphere, in Kelvin(!)\n')
        data.append(MaxIter.get()+' \t \t \t'+' maxItr = maximum number of iterations of the velocity solution (e.g., 80?)\n')
        data.append(AccConv.get()+' \t \t \t'+' okToQt = acceptable level of fractional change in RMS velocity which stops iteration\n')
        data.append(RefStress.get()+' \t \t \t'+' refStre = reference level of shear stress in lithosphere, for initiating linearization of rheology, in Pa\n')
        data.append(AccVel.get()+' \t \t \t'+' okDelV = acceptable level of velocity errors due to necessary limit on highest viscosity, in m/s\n')
        data.append(Out.get()+' \t \t \t'+' everyP = switch: Shall velocities of nodes be output in every iteration (for convergence studies)?\n')
        data.append('------------------------------------------------------------------------------------------\n')
        data.append('MEMO about possible values of iConve:\n')
        data.append('0 = lower mantle is static (with respect to AF)\n')
        data.append('1 = Hager and O\'Connell (1979) minimum-flow kinematic Model II\n')
        data.append('2 = Baumgardner (1988) dynamic model of Figure 7A-F, * 10. because Ra was too low\n')
        data.append('3 = PB2002 (Bird, 2003) (i.e., lower mantle flow is identical to surface velocities)\n')
        data.append('4 = PB2002 flow which drags continents; but NO drag on oceanic lithosphere\n')
        data.append('5 = drag on base of subduction forearc only (for local models)\n')
        data.append('6 = sense & traction from traction pole vector for each plate (use vTimes = 1.0)\n')
        data.append('\n')
        data.append('It is important within the parameter file that lines with multiple inputs seperated by a comma must not have a space before\n')
        data.append('or after the comma (e.g. 5.E8,5.E8  not  5.E8, 5.E8).')

      if diff:
        if os.path.exists(InDir+'/'+ParamFile):
          os.system('cp ' + InDir+'/'+ParamFile +' '+ InDir+'/'+'OLD_'+ParamFile)
        with open(InDir+'/'+ParamFile, 'w') as file:
          file.writelines(data)

      Param.quit()
      Param.destroy()


  if not Suppress:
    UnitsInfo()
    PlanetsInfo()
# Make Frame
  Param = tk.Toplevel()
  Param.title('Parameter file update')

# Labels
  tk.Label(Param, text='FAULT FRICTION COEFFICIENT', font=('Roman',12)).grid(row=0,sticky='w')
  tk.Label(Param, text='CONTINUUM FRICTION COEFFICIENT', font=('Roman',12)).grid(row=1,sticky='w')
  tk.Label(Param, text='BIOT COEFFICIENT (EFFICACY OF PORE PRESSURE)', font=('Roman',12)).grid(row=2,sticky='w')
  tk.Label(Param, text='BYERLY (0.-.99); FRACTIONAL STRENGTH REDUCTION OF MASTER FAULT', font=('Roman',12)).grid(row=3,sticky='w')
  tk.Label(Param, text='ACREEP (SHEAR STRESS COEFFICIENT OF CREEP LAW)', font=('Roman',12)).grid(row=4,sticky='w')
  tk.Label(Param, text='BCREEP (ACTIVATION ENERGY/N/GAS-CONSTANT) (IN KELVIN)', font=('Roman',12)).grid(row=5,sticky='w')
  tk.Label(Param, text='CCREEP (DERIVITIVE OF BCREEP WITH RESPECT TO DEPTH; CRUST/MANTLE)', font=('Roman',12)).grid(row=6,sticky='w')
  tk.Label(Param, text='DCREEP (MAXIMUM SHEAR STRESS; CRUST/MANTLE)', font=('Roman',12)).grid(row=7,sticky='w')
  tk.Label(Param, text='ECREEP (EXPONENT ON STRAIN-RATE IN CREEP-STRESS LAWS)=(1/N)', font=('Roman',12)).grid(row=8,sticky='w')
  tk.Label(Param, text='INTERCEPT AND SLOPE OF UPPER MANTLE ADIABAT (K, K/M)', font=('Roman',12)).grid(row=9,sticky='w')
  tk.Label(Param, text='DEPTH OF BASE OF ASTHENOSPHERE', font=('Roman',12)).grid(row=10,sticky='w')
  tk.Label(Param, text='AFRICAN PLATE REFERENCE FRAME', font=('Roman',12)).grid(row=11,sticky='w')
  tk.Label(Param, text='ICONVE:0=NONE;1=HAGER;2=BAUMGARDNER;3=PB2002;4=CONTINENTAL PB2002;5=forearc;6=inferred', font=('Roman',12)).grid(row=12,sticky='w')
  tk.Label(Param, text='TRHMAX (LIMIT ON BASAL TRACTION)', font=('Roman',12)).grid(row=13,sticky='w')
  tk.Label(Param, text='TAUMAX (DOWN-DIP INTEGRAL OF SUBDUCTION ZONE TRACTION; OCEAN\LAND)', font=('Roman',12)).grid(row=14,sticky='w')
  tk.Label(Param, text='RHOH2O (DENSITY OF WATER, AT P=0 AND T=SURFACE TEMPERATURE)', font=('Roman',12)).grid(row=15,sticky='w')
  tk.Label(Param, text='RHOBAR (MEAN DENSITY AT P=0 AND T=0; CRUST/MANTLE)', font=('Roman',12)).grid(row=16,sticky='w')
  tk.Label(Param, text='RHOAST (DENSITY OF ASTHENOSPHERE, AT P=0 AND AMBIENT T)', font=('Roman',12)).grid(row=17,sticky='w')
  tk.Label(Param, text='GRAVITATIONAL ACCELERATION', font=('Roman',12)).grid(row=18,sticky='w')
  tk.Label(Param, text='ONE KILOMETER, EXPRESSED IN CURRENT LENGTH UNITS', font=('Roman',12)).grid(row=19,sticky='w')
  tk.Label(Param, text='RADIUS OF THE PLANET', font=('Roman',12)).grid(row=20,sticky='w')
  tk.Label(Param, text='VOLUMETRIC THERMAL EXPANSION COEFFICIENT; CRUST/MANTLE', font=('Roman',12)).grid(row=21,sticky='w')
  tk.Label(Param, text='THERMAL CONDUCTIVITY; CRUST/MANTLE', font=('Roman',12)).grid(row=22,sticky='w')
  tk.Label(Param, text='RADIOACTIVE HEAT PRODUCTION, ON VOLUME (NOT MASS) BASIS', font=('Roman',12)).grid(row=23,sticky='w')
  tk.Label(Param, text='SURFACE TEMPERATURE, IN KELVIN', font=('Roman',12)).grid(row=24,sticky='w')
  tk.Label(Param, text='UPPER TEMPERATURE LIMITS, IN KELVIN; CRUST/MANTLE-LITHOSPHERE', font=('Roman',12)).grid(row=25,sticky='w')
  tk.Label(Param, text='MAXIMUM NUMBER OF ITERATIONS', font=('Roman',12)).grid(row=26,sticky='w')
  tk.Label(Param, text='ACCEPTABLE CONVERGENCE LEVEL (FRACTIONAL VELOCITY CHANGE)', font=('Roman',12)).grid(row=27,sticky='w')
  tk.Label(Param, text='REFERENCE LEVEL OF SHEAR STRESS', font=('Roman',12)).grid(row=28,sticky='w')
  tk.Label(Param, text='ACCEPTABLE LEVEL OF VELOCITY ERRORS (1 MM/A = 3.17E-11 M/S)', font=('Roman',12)).grid(row=29,sticky='w')
  tk.Label(Param, text='OUTPUT NODE VELOCITIES EVERY ITERATION? (FOR CONVERGENCE STUDIES)', font=('Roman',12)).grid(row=30,sticky='w')

# Entries
  FFric = tk.Entry(Param,font=('Roman',12))
  CFric = tk.Entry(Param,font=('Roman',12))
  Biot = tk.Entry(Param,font=('Roman',12))
  Byerly = tk.Entry(Param,font=('Roman',12))
  ACreepC = tk.Entry(Param,font=('Roman',12))
  ACreepM = tk.Entry(Param,font=('Roman',12))
  BCreepC = tk.Entry(Param,font=('Roman',12))
  BCreepM = tk.Entry(Param,font=('Roman',12))
  CCreepC = tk.Entry(Param,font=('Roman',12))
  CCreepM = tk.Entry(Param,font=('Roman',12))
  DCreepC = tk.Entry(Param,font=('Roman',12))
  DCreepM = tk.Entry(Param,font=('Roman',12))
  ECreep = tk.Entry(Param,font=('Roman',12))
  AdiabatC = tk.Entry(Param,font=('Roman',12))
  AdiabatM = tk.Entry(Param,font=('Roman',12))
  Zbasth = tk.Entry(Param,font=('Roman',12))
  AFRef = tk.Entry(Param,font=('Roman',12))
  Iconve1 = tk.Entry(Param,font=('Roman',12))
  Iconve2 = tk.Entry(Param,font=('Roman',12))
  TRHmax = tk.Entry(Param,font=('Roman',12))
  TAUmax = tk.Entry(Param,font=('Roman',12))
  RhoH2O = tk.Entry(Param,font=('Roman',12))
  RhoBarC = tk.Entry(Param,font=('Roman',12))
  RhoBarM = tk.Entry(Param,font=('Roman',12))
  RhoAst = tk.Entry(Param,font=('Roman',12))
  Grav = tk.Entry(Param,font=('Roman',12))
  OneKM = tk.Entry(Param,font=('Roman',12))
  Radius = tk.Entry(Param,font=('Roman',12))
  VolExpC = tk.Entry(Param,font=('Roman',12))
  VolExpM = tk.Entry(Param,font=('Roman',12))
  ThrConC = tk.Entry(Param,font=('Roman',12))
  ThrConM = tk.Entry(Param,font=('Roman',12))
  RadHeatC = tk.Entry(Param,font=('Roman',12))
  RadHeatM = tk.Entry(Param,font=('Roman',12))
  Temp = tk.Entry(Param,font=('Roman',12))
  UppManC = tk.Entry(Param,font=('Roman',12))
  UppManM = tk.Entry(Param,font=('Roman',12))
  MaxIter = tk.Entry(Param,font=('Roman',12))
  AccConv = tk.Entry(Param,font=('Roman',12))
  RefStress = tk.Entry(Param,font=('Roman',12))
  AccVel = tk.Entry(Param,font=('Roman',12))
  Out = tk.Entry(Param,font=('Roman',12))


# Entry Defaults
  if os.path.exists(InDir+'/'+ParamFile):
    with open(InDir+'/'+ParamFile, 'r') as file:
      data = file.readlines()

    FFric.insert(0, data[1].split()[0])
    CFric.insert(0, data[2].split()[0])
    Biot.insert(0, data[3].split()[0])
    Byerly.insert(0, data[4].split()[0])
    ACreepC.insert(0, data[5].split(',')[0])
    ACreepM.insert(0, data[5].split(',')[1].split()[0])
    BCreepC.insert(0, data[6].split(',')[0])
    BCreepM.insert(0, data[6].split(',')[1].split()[0])
    CCreepC.insert(0, data[7].split(',')[0])
    CCreepM.insert(0, data[7].split(',')[1].split()[0])
    DCreepC.insert(0, data[8].split(',')[0])
    DCreepM.insert(0, data[8].split(',')[1].split()[0])
    ECreep.insert(0, data[9].split()[0])
    AdiabatC.insert(0, data[10].split(',')[0])
    AdiabatM.insert(0, data[10].split(',')[1].split()[0])
    Zbasth.insert(0, data[11].split()[0])
    AFRef.insert(0, data[12].split()[0])
    Iconve1.insert(0, data[13].split(',')[0])
    Iconve2.insert(0, data[13].split(',')[1].split()[0])
    TRHmax.insert(0, data[14].split()[0])
    TAUmax.insert(0, data[15].split()[0])
    RhoH2O.insert(0, data[16].split()[0])
    RhoBarC.insert(0, data[17].split(',')[0])
    RhoBarM.insert(0, data[17].split(',')[1].split()[0])
    RhoAst.insert(0, data[18].split()[0])
    Grav.insert(0, data[19].split()[0])
    OneKM.insert(0, data[20].split()[0])
    Radius.insert(0, data[21].split()[0])
    VolExpC.insert(0, data[22].split(',')[0])
    VolExpM.insert(0, data[22].split(',')[1].split()[0])
    ThrConC.insert(0, data[23].split(',')[0])
    ThrConM.insert(0, data[23].split(',')[1].split()[0])
    RadHeatC.insert(0, data[24].split(',')[0])
    RadHeatM.insert(0, data[24].split(',')[1].split()[0])
    Temp.insert(0, data[25].split()[0])
    UppManC.insert(0, data[26].split(',')[0])
    UppManM.insert(0, data[26].split(',')[1].split()[0])
    MaxIter.insert(0, data[27].split()[0])
    AccConv.insert(0, data[28].split()[0])
    RefStress.insert(0, data[29].split()[0])
    AccVel.insert(0, data[30].split()[0])
    Out.insert(0, data[31].split()[0])

  else:
    NoParamFile()

    FFric.insert(0, '0.10')
    CFric.insert(0, '0.85')
    Biot.insert(0, '1.00')
    Byerly.insert(0, '0.00')
    ACreepC.insert(0, '2.3E9')
    ACreepM.insert(0, '9.5E4')
    BCreepC.insert(0, '4000.')
    BCreepM.insert(0, '18314.')
    CCreepC.insert(0, '0.')
    CCreepM.insert(0, '0.0171')
    DCreepC.insert(0, '5.E8')
    DCreepM.insert(0, '5.E8')
    ECreep.insert(0, '0.333333')
    AdiabatC.insert(0, '1412.')
    AdiabatM.insert(0, '6.1E-4')
    Zbasth.insert(0, '400.E3')
    AFRef.insert(0, 'AF')
    Iconve1.insert(0, '0')
    Iconve2.insert(0, '1.00')
    TRHmax.insert(0, '0.0')
    TAUmax.insert(0, '2.0E+12')
    RhoH2O.insert(0, '1032.')
    RhoBarC.insert(0, '2889.')
    RhoBarM.insert(0, '3332.')
    RhoAst.insert(0, '3125.')
    Grav.insert(0, '9.8')
    OneKM.insert(0, '1000.')
    Radius.insert(0, '6371000.')
    VolExpC.insert(0, '2.4E-5')
    VolExpM.insert(0, '3.94E-5')
    ThrConC.insert(0, '2.7')
    ThrConM.insert(0, '3.20')
    RadHeatC.insert(0, '3.5E-7')
    RadHeatM.insert(0, '3.2E-8')
    Temp.insert(0, '273.')
    UppManC.insert(0, '1223.')
    UppManM.insert(0, '1673.')
    MaxIter.insert(0, '50')
    AccConv.insert(0, '0.0005')
    RefStress.insert(0, '50.E6')
    AccVel.insert(0, '1.00E-11')
    Out.insert(0, 'F')

# Entry locations
  FFric.grid(row=0, column=2)
  CFric.grid(row=1, column=2)
  Biot.grid(row=2, column=2)
  Byerly.grid(row=3, column=2)
  ACreepC.grid(row=4, column=2)
  ACreepM.grid(row=4, column=3)
  BCreepC.grid(row=5, column=2)
  BCreepM.grid(row=5, column=3)
  CCreepC.grid(row=6, column=2)
  CCreepM.grid(row=6, column=3)
  DCreepC.grid(row=7, column=2)
  DCreepM.grid(row=7, column=3)
  ECreep.grid(row=8, column=2)
  AdiabatC.grid(row=9, column=2)
  AdiabatM.grid(row=9, column=3)
  Zbasth.grid(row=10, column=2)
  AFRef.grid(row=11, column=2)
  Iconve1.grid(row=12, column=2)
  Iconve2.grid(row=12, column=3)
  TRHmax.grid(row=13, column=2)
  TAUmax.grid(row=14, column=2)
  RhoH2O.grid(row=15, column=2)
  RhoBarC.grid(row=16, column=2)
  RhoBarM.grid(row=16, column=3)
  RhoAst.grid(row=17, column=2)
  Grav.grid(row=18, column=2)
  OneKM.grid(row=19, column=2)
  Radius.grid(row=20, column=2)
  VolExpC.grid(row=21, column=2)
  VolExpM.grid(row=21, column=3)
  ThrConC.grid(row=22, column=2)
  ThrConM.grid(row=22, column=3)
  RadHeatC.grid(row=23, column=2)
  RadHeatM.grid(row=23, column=3)
  Temp.grid(row=24, column=2)
  UppManC.grid(row=25, column=2)
  UppManM.grid(row=25, column=3)
  MaxIter.grid(row=26, column=2)
  AccConv.grid(row=27, column=2)
  RefStress.grid(row=28, column=2)
  AccVel.grid(row=29, column=2)
  Out.grid(row=30, column=2)

# Save button
  btnStart=tk.Button(Param,text='Save',height=1, width=4,relief=RAISED, justify=CENTER,command=Save).grid(row=29, column=3)

  Param.mainloop()


# Update input files list
def InputFiles():


  def InputFiles2():

    global InDir, InFile, data
    data = []


    def NoInputFile():
      string = ('No input files list file "'+InFile+'" found in the input directory "'+InDir+'", '
                'therefore these filenames will be prefilled with example names of files provided '
                'with the original ShellSet package. '
                'If you expect the file to exist then check for spelling errors. '
                'If the file name or input directory are incorrect then edit them at the '
                'beginning of the Python file "ShellSetGUI.py" and restart the GUI.')

      messagebox.showinfo('Input file list not found',string)


    def Save():
      global InDir, InFile,data
      diff = tk.BooleanVar(DataIn)
      diff = False
      fail = tk.BooleanVar(DataIn)
      fail = False

# Check entered files exist
# Global
      if not os.path.exists(InDir+'/'+ParamIn.get()) or len(ParamIn.get()) == 0: # Required
        messagebox.showerror('Error','Parameter file does not exist or empty')
        fail = True

      if not os.path.exists(InDir+'/'+GridIn.get()) or len(GridIn.get()) == 0: # Required
        messagebox.showerror('Error','Grid file does not exist or empty')
        fail = True

      if not os.path.exists(InDir+'/'+LRIn.get()) and LRIn.get() != 'X': # Optional
        messagebox.showerror('Error','Lithospheric Rheologies file does not exist')
        fail = True

# OrbData
      if not os.path.exists(InDir+'/'+ElevIn.get()) or len(ElevIn.get()) == 0:
        messagebox.showerror('Error','OrbData Elevation/Topography/Bathymetry file does not exist or empty')
        fail = True

      if not os.path.exists(InDir+'/'+SeaIn.get()) or len(SeaIn.get()) == 0:
        messagebox.showerror('Error','OrbData Seafloor age file does not exist or empty')
        fail = True

      if not os.path.exists(InDir+'/'+CrustIn.get()) or len(CrustIn.get()) == 0:
        messagebox.showerror('Error','OrbData Crustal thickness file does not exist or empty')
        fail = True

      if not os.path.exists(InDir+'/'+SWaveIn.get()) or len(SWaveIn.get()) == 0:
        messagebox.showerror('Error','OrbData S-wave anomaly file does not exist or empty')
        fail = True

# Shells
      if not os.path.exists(InDir+'/'+BcsInS.get()) or len(BcsInS.get()) == 0:
        messagebox.showerror('Error','Shells Boundary file does not exist or empty')
        fail = True

      if not os.path.exists(InDir+'/'+PltOutInS.get()) or len(PltOutInS.get()) == 0:
        messagebox.showerror('Error','Shells Plate outline file does not exist or empty')
        fail = True

      if not os.path.exists(InDir+'/'+PltBoundInS.get()) or len(PltBoundInS.get()) == 0:
        messagebox.showerror('Error','Shells Plate-Pair boundary file does not exist or empty')
        fail = True
# Optional
      if not os.path.exists(InDir+'/'+VelSolInS.get()) and VelSolInS.get() != 'X':
        messagebox.showerror('Error','Shells Approx. Velocity file does not exist')
        fail = True

      if not os.path.exists(InDir+'/'+MtlFlwInS.get()) and MtlFlwInS.get() != 'X':
        messagebox.showerror('Error','Shells Mantle flow file does not exist')
        fail = True

      if not os.path.exists(InDir+'/'+TrqInS.get()) and TrqInS.get() != 'X':
        messagebox.showerror('Error','Shells Torque & force file does not exist')
        fail = True

# Shells final
      if not os.path.exists(InDir+'/'+BcsInF.get()) or len(BcsInF.get()) == 0:
        messagebox.showerror('Error','Shells Final Boundary file does not exist or empty')
        fail = True

      if not os.path.exists(InDir+'/'+PltOutInF.get()) or len(PltOutInF.get()) == 0:
        messagebox.showerror('Error','Shells Final Plate outline file does not exist or empty')
        fail = True

      if not os.path.exists(InDir+'/'+PltBoundInF.get()) or len(PltBoundInF.get()) == 0:
        messagebox.showerror('Error','Shells Final Plate-Pair boundary file does not exist or empty')
        fail = True
# Optional
      if not os.path.exists(InDir+'/'+VelSolInF.get()) and VelSolInF.get() != 'X':
        messagebox.showerror('Error','Shells Final Approx. Velocity file does not exist')
        fail = True

      if not os.path.exists(InDir+'/'+MtlFlwInF.get()) and MtlFlwInF.get() != 'X':
        messagebox.showerror('Error','Shells Final Mantle flow file does not exist')
        fail = True

      if not os.path.exists(InDir+'/'+TrqInF.get()) and TrqInF.get() != 'X':
        messagebox.showerror('Error','Shells Final Torque & force file does not exist')
        fail = True

# OrbScore
      if GVLog.get():
        if not os.path.exists(InDir+'/'+GVIn.get()) and GVIn.get() != 'X':
          messagebox.showerror('Error','File used for Geodetic Velocity scoring does not exist')
          fail = True

      if SDLog.get():
        if not os.path.exists(InDir+'/'+SDIn.get()) and SDIn.get() != 'X':
          messagebox.showerror('Error','File used for Stress Direction scoring does not exist')
          fail = True

      if SSRLog.get():
        if not os.path.exists(InDir+'/'+SSRIn.get()) and SSRIn.get() != 'X':
          messagebox.showerror('Error','File used for Sea-floor Spreading Rate scoring does not exist')
          fail = True

      if FSRLog.get():
        if not os.path.exists(InDir+'/'+FSRIn.get()) and FSRIn.get() != 'X':
          messagebox.showerror('Error','File used for Fault Slip Rate scoring does not exist')
          fail = True

      if SCLog.get():
        if not os.path.exists(InDir+'/'+SCIn.get()) and SCIn.get() != 'X':
          messagebox.showerror('Error','File used for Seismic Catalogue scoring does not exist')
          fail = True

      if SALog.get():
        if not os.path.exists(InDir+'/'+SAIn.get()) or not os.path.exists(InDir+'/'+SCIn.get()):
          if SAIn.get() != 'X' or SCIn.get() != 'X':
            messagebox.showerror('Error','One or more files used for Seismic Anisotropy scoring does not exist')
            fail = True


      if not fail:
        if os.path.exists(InDir+'/'+InFile):
# Global
          if ParamIn.get() != data[1].split()[0]:
            data[1] = ParamIn.get()+' \t '+' '.join(data[1].split()[1:])+'\n'
            diff = True
          if GridIn.get() != data[2].split()[0]:
            data[2] = GridIn.get()+' \t '+' '.join(data[2].split()[1:])+'\n'
            diff = True
          if LRIn.get() != data[3].split()[0]:
            data[3] = LRIn.get()+' \t '+' '.join(data[3].split()[1:])+'\n'
            diff = True
# OrbData
          if ElevIn.get() != data[7].split()[0]:
            data[7] = ElevIn.get()+' \t '+' '.join(data[7].split()[1:])+'\n'
            diff = True
          if SeaIn.get() != data[8].split()[0]:
            data[8] = SeaIn.get()+' \t '+' '.join(data[8].split()[1:])+'\n'
            diff = True
          if CrustIn.get() != data[9].split()[0]:
            data[9] = CrustIn.get()+' \t '+' '.join(data[9].split()[1:])+'\n'
            diff = True
          if SWaveIn.get() != data[10].split()[0]:
            data[10] = SWaveIn.get()+' \t '+' '.join(data[10].split()[1:])+'\n'
            diff = True
# Shells
          if BcsInS.get() != data[14].split()[0]:
            data[14] = BcsInS.get()+' \t '+' '.join(data[14].split()[1:])+'\n'
            diff = True
          if PltOutInS.get() != data[15].split()[0]:
            data[15] = PltOutInS.get()+' \t '+' '.join(data[15].split()[1:])+'\n'
            diff = True
          if PltBoundInS.get() != data[16].split()[0]:
            data[16] = PltBoundInS.get()+' \t '+' '.join(data[16].split()[1:])+'\n'
            diff = True
# Optional
          if VelSolInS.get() != data[18].split()[0]:
            data[18] = VelSolInS.get()+' \t '+' '.join(data[18].split()[1:])+'\n'
            diff = True
          if MtlFlwInS.get() != data[19].split()[0]:
            data[19] = MtlFlwInS.get()+' \t '+' '.join(data[19].split()[1:])+'\n'
            diff = True
          if TrqInS.get() != data[20].split()[0]:
            data[20] = TrqInS.get()+' \t '+' '.join(data[20].split()[1:])+'\n'
            diff = True
# Shells final
          if BcsInF.get() != data[24].split()[0]:
            data[24] = BcsInF.get()+' \t '+' '.join(data[24].split()[1:])+'\n'
            diff = True
          if PltOutInF.get() != data[25].split()[0]:
            data[25] = PltOutInF.get()+' \t '+' '.join(data[25].split()[1:])+'\n'
            diff = True
          if PltBoundInF.get() != data[26].split()[0]:
            data[26] = PltBoundInF.get()+' \t '+' '.join(data[26].split()[1:])+'\n'
            diff = True
# Optional
          if VelSolInF.get() != data[28].split()[0]:
            data[28] = VelSolInF.get()+' \t '+' '.join(data[28].split()[1:])+'\n'
            diff = True
          if MtlFlwInF.get() != data[29].split()[0]:
            data[29] = MtlFlwInF.get()+' \t '+' '.join(data[29].split()[1:])+'\n'
            diff = True
          if TrqInF.get() != data[30].split()[0]:
            data[30] = TrqInF.get()+' \t '+' '.join(data[30].split()[1:])+'\n'
            diff = True
# OrbScore
          if GVLog.get():
            if GVIn.get() != data[34].split()[0]:
              data[34] = GVIn.get()+' \t '+' '.join(data[34].split()[1:])+'\n'
              diff = True
          if SDLog.get():
            if SDIn.get() != data[35].split()[0]:
              data[35] = SDIn.get()+' \t '+' '.join(data[35].split()[1:])+'\n'
              diff = True
          if FSRLog.get():
            if FSRIn.get() != data[36].split()[0]:
              data[36] = FSRIn.get()+' \t '+' '.join(data[36].split()[1:])+'\n'
              diff = True
          if SSRLog.get():
            if SSRIn.get() != data[37].split()[0]:
              data[37] = SSRIn.get()+' \t '+' '.join(data[37].split()[1:])+'\n'
              diff = True
          if SCLog.get():
            if SCIn.get() != data[38].split()[0]:
              data[38] = SCIn.get()+' \t '+' '.join(data[38].split()[1:])+'\n'
              diff = True
          if SALog.get():
            if SAIn.get() != data[39].split()[0]:
              data[39] = SAIn.get()+' \t '+' '.join(data[39].split()[1:])+'\n'
              diff = True
            if SCIn.get() != data[38].split()[0]:
              data[38] = SCIn.get()+' \t '+' '.join(data[38].split()[1:])+'\n'
              diff = True

            if SAIn.get() == 'X'  or SCIn.get() == 'X':
              data[39] = 'X'+' \t '+' '.join(data[39].split()[1:])+'\n'

        else:
          diff = True
# Global
          data.append('All:\n')
          data.append(ParamIn.get()+' \t \t \t'+' Parameter input file\n')
          data.append(GridIn.get()+' \t \t \t'+' Grid input file\n')
          data.append(LRIn.get()+' \t \t \t'+' non-default Lithospheric Rheologies\n')
# OrbData
          data.append('------------------------------------------------------------------------------------------\n')
          data.append('\n')
          data.append('OrbData:\n')
          data.append(ElevIn.get()+' \t \t \t'+' elevation/topography/bathymetry\n')
          data.append(SeaIn.get()+' \t \t \t'+' age of seafloor\n')
          data.append(CrustIn.get()+' \t \t \t'+' thickness of crust\n')
          data.append(SWaveIn.get()+' \t \t \t'+' S-wave travel-time anomaly in upper mantle\n')
# Shells
          data.append('------------------------------------------------------------------------------------------\n')
          data.append('\n')
          data.append('Shells:\n')
          data.append(BcsInS.get()+' \t \t \t'+' Boundary condition file\n')
          data.append(PltOutInS.get()+' \t \t \t'+' Outlines of Plates\n')
          data.append(PltBoundInS.get()+' \t \t \t'+' Plate-Pair boundaries\n')
          data.append('------------------------------\n')
# Optional
          data.append(VelSolInS.get()+' \t \t \t'+' Approx Vel solution\n')
          data.append(MtlFlwInS.get()+' \t \t \t'+' Mantle Flow -> HOC79ii.dig (iConve=1), Baum887.dig (iConve=2), PB2002_plates.dig (iConve=3 or 4)\n')
          data.append(TrqInS.get()+' \t \t \t'+' Torque & Force balance -> Only if iConve=6\n')
# Shells final
          data.append('------------------------------------------------------------------------------------------\n')
          data.append('\n')
          data.append('Shells (Final run):\n')
          data.append(BcsInF.get()+' \t \t \t'+' Boundary condition file\n')
          data.append(PltOutInF.get()+' \t \t \t'+' Outlines of Plates\n')
          data.append(PltBoundInF.get()+' \t \t \t'+' Plate-Pair boundaries\n')
          data.append('------------------------------\n')
# Optional
          data.append(VelSolInF.get()+' \t \t \t'+' Approx Vel solution\n')
          data.append(MtlFlwInF.get()+' \t \t \t'+' Mantle Flow -> HOC79ii.dig (iConve=1), Baum887.dig (iConve=2), PB2002_plates.dig (iConve=3 or 4)\n')
          data.append(TrqInF.get()+' \t \t \t'+' Torque & Force balance -> Only if iConve=6\n')
# OrbScore
          data.append('------------------------------------------------------------------------------------------\n')
          data.append('\n')
          data.append('OrbScore:\n')
          data.append(GVIn.get()+' \t \t \t'+' (.gps | Geodetic Velocity)\n')
          data.append(SDIn.get()+' \t \t \t'+' (.dat | Stress Direction)\n')
          data.append(FSRIn.get()+' \t \t \t'+' (.txt | Fault Slip Rate)\n')
          data.append(SSRIn.get()+' \t \t \t'+' (.dat | Seafloor Spreading Rates)\n')
          data.append(SCIn.get()+' \t \t \t'+' (.eqc | Smoothed Seismic Correlation)\n')
          data.append(SAIn.get()+' \t \t \t'+' (.dat | Seismic Anisotropy | Only if SC)\n')
          data.append('------------------------------------------------------------------------------------------')

        if diff:
          if os.path.exists(InDir+'/'+InFile):
            os.system('cp ' +InDir+'/'+InFile+' ' +InDir+'/'+'OLD_'+InFile)
          with open(InDir+'/'+InFile, 'w') as file:
            file.writelines(data)

        DataIn.quit()
        DataIn.destroy()

# Make Frame
    DataIn = tk.Toplevel()
    DataIn.title('Input file list')

# Labels
    tk.Label(DataIn, text='Required by all program parts:', font=('Roman',15)).grid(row=0,sticky='w')
    tk.Label(DataIn, text='Parameter file', font=('Roman',12)).grid(row=1,sticky='w')
    tk.Label(DataIn, text='Grid file', font=('Roman',12)).grid(row=2,sticky='w')
    tk.Label(DataIn, text='Non-default Lithospheric Rheologies (opt)', font=('Roman',12)).grid(row=3,sticky='w')

    tk.Label(DataIn, text='\nOrbData:', font=('Roman',15)).grid(row=4,sticky='w')
    tk.Label(DataIn, text='Elevation/Topography/Bathymetry', font=('Roman',12)).grid(row=5,sticky='w')
    tk.Label(DataIn, text='Seafloor age', font=('Roman',12)).grid(row=6,sticky='w')
    tk.Label(DataIn, text='Crustal thickness', font=('Roman',12)).grid(row=7,sticky='w')
    tk.Label(DataIn, text='S-wave travel-time anomaly in upper mantle', font=('Roman',12)).grid(row=8,sticky='w')

    tk.Label(DataIn, text='\nShells:', font=('Roman',15)).grid(row=9,sticky='w')
    tk.Label(DataIn, text='Boundary condition', font=('Roman',12)).grid(row=10,sticky='w')
    tk.Label(DataIn, text='Outlines of Plates', font=('Roman',12)).grid(row=11,sticky='w')
    tk.Label(DataIn, text='Plate-Pair boundaries', font=('Roman',12)).grid(row=12,sticky='w')
    tk.Label(DataIn, text='Approx Vel solution (opt)', font=('Roman',12)).grid(row=13,sticky='w')
    tk.Label(DataIn, text='Mantle Flow (opt)', font=('Roman',12)).grid(row=14,sticky='w')
    tk.Label(DataIn, text='Torque & Force balance (opt)', font=('Roman',12)).grid(row=15,sticky='w')

    tk.Label(DataIn, text='\nShells (Final OR Non-Iterating run):', font=('Roman', 16)).grid(row=16,sticky='w')
    tk.Label(DataIn, text='Boundary condition', font=('Roman',12)).grid(row=17,sticky='w')
    tk.Label(DataIn, text='Outlines of Plates', font=('Roman',12)).grid(row=18,sticky='w')
    tk.Label(DataIn, text='Plate-Pair boundaries', font=('Roman',12)).grid(row=19,sticky='w')
    tk.Label(DataIn, text='Approx Vel solution (opt)', font=('Roman',12)).grid(row=20,sticky='w')
    tk.Label(DataIn, text='Mantle Flow (opt)', font=('Roman',12)).grid(row=21,sticky='w')
    tk.Label(DataIn, text='Torque & Force balance (opt)', font=('Roman',12)).grid(row=22,sticky='w')

    if GVLog.get() or SSRLog.get() or SDLog.get() or FSRLog.get() or SCLog.get() or SALog.get():
      tk.Label(DataIn, text='\nOrbScore:', font=('Roman',15)).grid(row=23,sticky='w')
      if GVLog.get():
        tk.Label(DataIn, text='Geodetic Velocities', font=('Roman',12)).grid(row=24,sticky='w')
      if SSRLog.get():
        tk.Label(DataIn, text='Sea-floor Spreading Rates', font=('Roman',12)).grid(row=27,sticky='w')
      if SDLog.get():
        tk.Label(DataIn, text='Stress Directions', font=('Roman',12)).grid(row=25,sticky='w')
      if FSRLog.get():
        tk.Label(DataIn, text='Fault Slip Rates', font=('Roman',12)).grid(row=26,sticky='w')
      if SCLog.get():
        tk.Label(DataIn, text='Seismic Catalogue', font=('Roman',12)).grid(row=28,sticky='w')
      if SALog.get():
        tk.Label(DataIn, text='Seismic Catalogue', font=('Roman',12)).grid(row=28,sticky='w')
        tk.Label(DataIn, text='Upper-mantle anisotropy', font=('Roman',12)).grid(row=29,sticky='w')

# Entries
    ParamIn = tk.Entry(DataIn, font=('Roman',12), width=40)
    GridIn  = tk.Entry(DataIn, font=('Roman',12), width=40)
    LRIn    = tk.Entry(DataIn, font=('Roman',12), width=40)

    ElevIn  = tk.Entry(DataIn, font=('Roman',12), width=40)
    SeaIn   = tk.Entry(DataIn, font=('Roman',12), width=40)
    CrustIn = tk.Entry(DataIn, font=('Roman',12), width=40)
    SWaveIn = tk.Entry(DataIn, font=('Roman',12), width=40)

    BcsInS      = tk.Entry(DataIn, font=('Roman',12), width=40)
    PltOutInS   = tk.Entry(DataIn, font=('Roman',12), width=40)
    PltBoundInS = tk.Entry(DataIn, font=('Roman',12), width=40)
    VelSolInS   = tk.Entry(DataIn, font=('Roman',12), width=40)
    MtlFlwInS   = tk.Entry(DataIn, font=('Roman',12), width=40)
    TrqInS      = tk.Entry(DataIn, font=('Roman',12), width=40)

    BcsInF      = tk.Entry(DataIn, font=('Roman',12), width=40)
    PltOutInF   = tk.Entry(DataIn, font=('Roman',12), width=40)
    PltBoundInF = tk.Entry(DataIn, font=('Roman',12), width=40)
    VelSolInF   = tk.Entry(DataIn, font=('Roman',12), width=40)
    MtlFlwInF   = tk.Entry(DataIn, font=('Roman',12), width=40)
    TrqInF      = tk.Entry(DataIn, font=('Roman',12), width=40)

    GVIn  = tk.Entry(DataIn, font=('Roman',12), width=40)
    SDIn  = tk.Entry(DataIn, font=('Roman',12), width=40)
    FSRIn = tk.Entry(DataIn, font=('Roman',12), width=40)
    SSRIn = tk.Entry(DataIn, font=('Roman',12), width=40)
    SCIn  = tk.Entry(DataIn, font=('Roman',12), width=40)
    SAIn  = tk.Entry(DataIn, font=('Roman',12), width=40)

# Entry Defaults
    if os.path.exists(InDir+'/'+InFile):
      with open(InDir+'/'+InFile, 'r') as file:
        data = file.readlines()
        ParamIn.insert(0, data[1].split()[0])
        GridIn.insert(0, data[2].split()[0])
        LRIn.insert(0, data[3].split()[0])

        ElevIn.insert(0, data[7].split()[0])
        SeaIn.insert(0, data[8].split()[0])
        CrustIn.insert(0, data[9].split()[0])
        SWaveIn.insert(0, data[10].split()[0])

        BcsInS.insert(0, data[14].split()[0])
        PltOutInS.insert(0, data[15].split()[0])
        PltBoundInS.insert(0, data[16].split()[0])
        VelSolInS.insert(0, data[18].split()[0])
        MtlFlwInS.insert(0, data[19].split()[0])
        TrqInS.insert(0, data[20].split()[0])

        BcsInF.insert(0, data[24].split()[0])
        PltOutInF.insert(0, data[25].split()[0])
        PltBoundInF.insert(0, data[26].split()[0])
        VelSolInF.insert(0, data[28].split()[0])
        MtlFlwInF.insert(0, data[29].split()[0])
        TrqInF.insert(0, data[30].split()[0])

        GVIn.insert(0, data[34].split()[0])
        SDIn.insert(0, data[35].split()[0])
        FSRIn.insert(0, data[36].split()[0])
        SSRIn.insert(0, data[37].split()[0])
        SCIn.insert(0, data[38].split()[0])
        SAIn.insert(0, data[39].split()[0])

    else:
      NoInputFile()

      ParamIn.insert(0, 'iEarth5-049.in')
      GridIn.insert(0, 'Earth5R.feg')
      LRIn.insert(0, 'X')

      ElevIn.insert(0, 'ETOPO20.grd')
      SeaIn.insert(0, 'age_1p5.grd')
      CrustIn.insert(0, 'CRUST2.grd')
      SWaveIn.insert(0, 'delta_ts.grd')

      BcsInS.insert(0, 'Earth5R-type4AplusA.bcs')
      PltOutInS.insert(0, 'PB2002_plates.dig')
      PltBoundInS.insert(0, 'PB2002_boundaries.dig')
      VelSolInS.insert(0, 'X')
      MtlFlwInS.insert(0, 'X')
      TrqInS.insert(0, 'X')

      BcsInF.insert(0, 'Earth5R-type4A.bcs')
      PltOutInF.insert(0, 'PB2002_plates.dig')
      PltBoundInF.insert(0, 'PB2002_boundaries.dig')
      VelSolInF.insert(0, 'X')
      MtlFlwInF.insert(0, 'X')
      TrqInF.insert(0, 'X')

      GVIn.insert(0, 'GPS2014_selected_subset.gps')
      SDIn.insert(0, 'robust_interpolated_stress_2023.dat')
      FSRIn.insert(0, 'aggregated_offset_rates.dig')
      SSRIn.insert(0, 'magnetic_PB2002.dat')
      SCIn.insert(0, 'GCMT_shallow_m5p7_1977-2017.eqc')
      SAIn.insert(0, 'Becker_SKS_2-degree_20220807-selected.dat')


# Entry locations
    ParamIn.grid(row=1, column=2)
    GridIn.grid(row=2, column=2)
    LRIn.grid(row=3, column=2)

    ElevIn.grid(row=5, column=2)
    SeaIn.grid(row=6, column=2)
    CrustIn.grid(row=7, column=2)
    SWaveIn.grid(row=8, column=2)

    BcsInS.grid(row=10, column=2)
    PltOutInS.grid(row=11, column=2)
    PltBoundInS.grid(row=12, column=2)
    VelSolInS.grid(row=13, column=2)
    MtlFlwInS.grid(row=14, column=2)
    TrqInS.grid(row=15, column=2)

    BcsInF.grid(row=17, column=2)
    PltOutInF.grid(row=18, column=2)
    PltBoundInF.grid(row=19, column=2)
    VelSolInF.grid(row=20, column=2)
    MtlFlwInF.grid(row=21, column=2)
    TrqInF.grid(row=22, column=2)


    if GVLog.get() or SSRLog.get() or SDLog.get() or FSRLog.get() or SCLog.get() or SALog.get():
      if GVLog.get():
        GVIn.grid(row=24, column=2)
      if SSRLog.get():
        SSRIn.grid(row=27, column=2)
      if SDLog.get():
        SDIn.grid(row=25, column=2)
      if FSRLog.get():
        FSRIn.grid(row=26, column=2)
      if SCLog.get():
        SCIn.grid(row=28, column=2)
      if SALog.get():
        SCIn.grid(row=28, column=2)
        SAIn.grid(row=29, column=2)

# Save button
    btnStart=tk.Button(DataIn,text='Save',height=1, width=4,relief=RAISED, justify=CENTER,command=Save, font=('Roman',12)).grid(row=32, column=2)

# Close OrbScore option frame
    ScoreOpt.quit()
    ScoreOpt.destroy()

    DataIn.mainloop()


# Make Frame
  ScoreOpt = tk.Toplevel()
  ScoreOpt.title('OrbScore options')

  GVLog  = tk.BooleanVar(ScoreOpt)
  SSRLog = tk.BooleanVar(ScoreOpt)
  SDLog  = tk.BooleanVar(ScoreOpt)
  FSRLog = tk.BooleanVar(ScoreOpt)
  SCLog  = tk.BooleanVar(ScoreOpt)
  SALog  = tk.BooleanVar(ScoreOpt)

# Labels
  tk.Label(ScoreOpt, text='Make changes to OrbScore options?', font=('Roman',15)).grid(row=1)
  tk.Label(ScoreOpt, text='Geodetic velocities', font=('Roman',12)).grid(row=2,sticky='w')
  tk.Label(ScoreOpt, text='Sea-floor Spreading Rates', font=('Roman',12)).grid(row=3,sticky='w')
  tk.Label(ScoreOpt, text='Stress Directions', font=('Roman',12)).grid(row=4,sticky='w')
  tk.Label(ScoreOpt, text='Fault Slip Rates', font=('Roman',12)).grid(row=5,sticky='w')
  tk.Label(ScoreOpt, text='Seismic Catalogue', font=('Roman',12)).grid(row=6,sticky='w')
  tk.Label(ScoreOpt, text='Seismic Anisotropy', font=('Roman',12)).grid(row=7,sticky='w')

# Check boxes
  GV  = tk.Checkbutton(ScoreOpt, text='Yes/No',variable=GVLog, onvalue=True, offvalue=False, font=('Roman',12))
  SSR = tk.Checkbutton(ScoreOpt, text='Yes/No',variable=SSRLog, onvalue=True, offvalue=False, font=('Roman',12))
  SD  = tk.Checkbutton(ScoreOpt, text='Yes/No',variable=SDLog, onvalue=True, offvalue=False, font=('Roman',12))
  FSR = tk.Checkbutton(ScoreOpt, text='Yes/No',variable=FSRLog, onvalue=True, offvalue=False, font=('Roman',12))
  SC  = tk.Checkbutton(ScoreOpt, text='Yes/No',variable=SCLog, onvalue=True, offvalue=False, font=('Roman',12))
  SA  = tk.Checkbutton(ScoreOpt, text='Yes/No',variable=SALog, onvalue=True, offvalue=False, font=('Roman',12))

# Entry locations
  GV.grid(row=2, column=1)
  SSR.grid(row=3, column=1)
  SD.grid(row=4, column=1)
  FSR.grid(row=5, column=1)
  SC.grid(row=6, column=1)
  SA.grid(row=7, column=1)

# Enter button
  btnStart=tk.Button(ScoreOpt,text='Enter',height=1, width=5,relief=RAISED, justify=CENTER,command=InputFiles2, font=('Roman',12)).grid(row=9, column=1)

  ScoreOpt.mainloop()


# List input file update
def ListInput():

  global InDir,VarCount,data,info
  data = []
  info = []
  VarCount = 0

# Input information about each variable
  def ListInput2():
    global InDir, ModTarget, VarCount, VarCountOld, Names_Used, VarTarget
    Names_Used = []
    ModTarget = NumMods.get()
    VarTarget = VarNum.get()

# Update saved variable count
    def Count():
      global VarCount,VarCountOld
      VarCountOld = VarCount
      VarCount += 1

# Save button work
    def Save():

# Clear the Entry Widget Content
      def clear_text():
        VarName.delete(0, END)
        VarValues.delete(1.0, END)


      global InDir, Names_Used
      fail = False

      Names = ['fFric','cFric','Biot','Byerly','aCreep_C','aCreep_M','bCreep_C','bCreep_M', # Approved variable list
               'cCreep_C','cCreep_M','dCreep_C','dCreep_M','eCreep','tAdiab','gradie','zBAsth','pltRef','iConve',
               'trHMax','tauMax','tauMax_S','tauMax_L','rhoH2O','rhoBar_C','rhoBar_M','rhoAst','gMean','oneKm',
               'radius','alphaT_C','alphaT_M','conduc_C','conduc_M','radio_C','radio_M','tSurf','temLim_C','temLim_M']

      if len(VarName.get()) != 0: # Check Variable name is allowed
        if VarName.get() not in Names:
          string = ('Variable name is not a Shells variable, choose from (case sensitive):\n'
                    'fFric, cFric, Biot, Byerly, aCreep_C, aCreep_M, bCreep_C, bCreep_M, cCreep_C, cCreep_M,'
                    'dCreep_C, dCreep_M, eCreep, tAdiab, gradie, zBAsth, pltRef, iConve, trHMax, tauMax,'
                    'tauMax_S, tauMax_L, rhoH2O, rhoBar_C, rhoBar_M, rhoAst, gMean, oneKm, radius, alphaT_C,'
                    'alphaT_M, conduc_C, conduc_M, radio_C, radio_M, tSurf, temLim_C, temLim_M')
          messagebox.showerror('Error',string)
          fail = True
# Check Variable name not used
        if VarName.get() in Names_Used:
          messagebox.showerror('Error','Variable has already been used')
          fail = True
      else:
        messagebox.showerror('Error','Nothing entered for variable name')
        fail = True
# Check enough models entered
      VarVals = [i.strip()+'\n' for i in VarValues.get(1.0,END).splitlines()]
      if len(VarVals) != int(ModTarget):
        messagebox.showerror('Error','Not enough models entered')
        fail = True

# Reset variable counter & go again
      if fail:
        global VarCount,VarCountOld
        VarCount = VarCountOld
      else:
# Updatee Names lists
        Names_Used.append(VarName.get())
        Names.remove(VarName.get())
# Save variable info to data
        global data,info
        if VarCount == 1: # only on first call
          info.append(VarNum.get() +','+ NumMods.get() +'\n')
          info.append(VarName.get() +'\n')
          data.append('---------- \n')
          data += VarVals
          data.append('---------- \n')
# Close old frame to avoid confusion
          ListRoot.quit()
          ListRoot.destroy()
# Clear input
          clear_text()
          List.title('List input variable ' + str(VarCount+1))
        else:
          info.append(VarName.get()+' \n')
          data += VarVals
          data.append('---------- \n')
# Save to file & close
          if VarCount == int(VarTarget):
            with open(InDir+'/'+'ListInput.in', 'w') as file:
              file.writelines(info)
              file.writelines(data)
            List.quit()
            List.destroy()
          else:
# Clear input
            clear_text()
            List.title('List input variable ' + str(VarCount+1))


# Check previously entered variables are correct
    fail = False
    if not VarNum.get().isnumeric():
      messagebox.showerror('Error','Number of variables must be an integer')
      fail = True
    if not NumMods.get().isnumeric():
      messagebox.showerror('Error','Number of models must be an integer')
      fail = True


# Ask for variable information
    if not fail:
# Initial frame
      List = tk.Toplevel()
      List.title('List input variable ' + str(VarCount+1))

# Labels
      tk.Label(List, text='Variable Name ', font=('Roman',12)).grid(row=0,sticky='w')
      ttk.Label(List, text='Variable Values :', font=('Roman',12)).grid(row=3,sticky='w',padx=10,pady=25)

# Entries
      VarName = tk.Entry(List)
      VarValues = tk.Text(List, width=20, height=ModTarget)

# Entry locations
      VarName.grid(row=0, column=2)
      VarValues.grid(column=2, row=3)

# Save button
      btnStart=tk.Button(List,text='Save',height=1, width=4,relief=RAISED, justify=CENTER,command=lambda: [Count(), Save()], font=('Roman',12)).grid(row=4, column=2)

      List.mainloop()


# Make Frame
  ListRoot = tk.Toplevel()
  ListRoot.title('Listed models input file')

# Labels
  tk.Label(ListRoot, text='Number of variables', font=('Roman',12)).grid(row=0,sticky='w')
  tk.Label(ListRoot, text='Number of models', font=('Roman',12)).grid(row=1,sticky='w')

# Entries
  VarNum = tk.Entry(ListRoot, font=('Roman',12))
  NumMods = tk.Entry(ListRoot, font=('Roman',12))

# Entry locations
  VarNum.grid(row=0, column=2)
  NumMods.grid(row=1, column=2)

# Enter button
  btnStart=tk.Button(ListRoot,text='Enter',height=1, width=5,relief=RAISED, justify=CENTER,command=ListInput2, font=('Roman',12)).grid(row=2, column=2)

  ListRoot.mainloop()


# Grid Search input file update
def GridInput():

  global InDir, VarCount, data
  data = []
  VarCount = 0

# Input information about each variable
  def GridInput2():
    global InDir, VarTarget, VarCount, VarCountOld, Names_Used
    Names_Used = []
    VarTarget = VarNum.get()
    # VarCount = 0

# Update saved variable count
    def Count():
      global VarCount,VarCountOld
      VarCountOld = VarCount
      VarCount += 1

# Save button work
    def Save():

# Clear the Entry Widget Content
      def clear_text():
        VarName.delete(0, END)
        VarMin.delete(0, END)
        VarMax.delete(0, END)
        VarMods.delete(0, END)


      global InDir, Names_Used
# Check previously entered variables are correct
      fail = False

      Names = ['fFric','cFric','Biot','Byerly','aCreep_C','aCreep_M','bCreep_C','bCreep_M', # Approved variable list
               'cCreep_C','cCreep_M','dCreep_C','dCreep_M','eCreep','tAdiab','gradie','zBAsth','pltRef','iConve',
               'trHMax','tauMax','tauMax_S','tauMax_L','rhoH2O','rhoBar_C','rhoBar_M','rhoAst','gMean','oneKm',
               'radius','alphaT_C','alphaT_M','conduc_C','conduc_M','radio_C','radio_M','tSurf','temLim_C','temLim_M']

      if len(VarName.get()) != 0: # Check Variable name is allowed
        if VarName.get() not in Names:
          string = ('Variable name is not a Shells variable, choose from (case sensitive):\n'
                    'fFric, cFric, Biot, Byerly, aCreep_C, aCreep_M, bCreep_C, bCreep_M, cCreep_C, cCreep_M,'
                    'dCreep_C, dCreep_M, eCreep, tAdiab, gradie, zBAsth, pltRef, iConve, trHMax, tauMax,'
                    'tauMax_S, tauMax_L, rhoH2O, rhoBar_C, rhoBar_M, rhoAst, gMean, oneKm, radius, alphaT_C,'
                    'alphaT_M, conduc_C, conduc_M, radio_C, radio_M, tSurf, temLim_C, temLim_M')
          messagebox.showerror('Error',string)
          fail = True
# Check Variable name not used
        if VarName.get() in Names_Used:
          messagebox.showerror('Error','Variable has already been used')
          fail = True
      else:
        messagebox.showerror('Error','Nothing entered for variable name')
        fail = True
# Min - float
      if len(VarMin.get()) != 0 :
        try :
          float(VarMin.get())
        except :
          messagebox.showerror('Error','Variable minimum must be numeric')
          fail = True
      else:
        messagebox.showerror('Error','Nothing entered for variable minimum')
        fail = True
# Max - float
      if len(VarMax.get()) != 0 :
        try :
          float(VarMax.get())
        except :
          messagebox.showerror('Error','Variable maximum must be numeric')
          fail = True
      else:
        messagebox.showerror('Error','Nothing entered for variable maximum')
        fail = True
# NumMods - integer
      if not VarMods.get().isnumeric():
        messagebox.showerror('Error','Number of models must be an integer')
        fail = True

# Reset variable counter & go again
      if fail:
        global VarCount,VarCountOld
        VarCount = VarCountOld
      else:
# Updatee Names lists
        Names_Used.append(VarName.get())
        Names.remove(VarName.get())
# Save variable info to data
        global data
        if VarCount == 1: # only on first call
          data.append(VarNum.get() +','+ VarCells.get() +','+ VarLvls.get() +'\n')
          data.append('---------- \n')
          data.append(VarName.get()+' \n')
          data.append(VarMin.get()+' \n')
          data.append(VarMax.get()+' \n')
          data.append(VarMods.get()+' \n')
          data.append('---------- \n')
# Close old frame to avoid confusion
          GridRoot.quit()
          GridRoot.destroy()
# Clear input
          clear_text()
          Grid.title('Grid Search variable ' + str(VarCount+1))
        else:
          data.append(VarName.get()+' \n')
          data.append(VarMin.get()+' \n')
          data.append(VarMax.get()+' \n')
          data.append(VarMods.get()+' \n')
          data.append('---------- \n')
# Save to file & close
          if VarCount == int(VarTarget):
            with open(InDir+'/'+'GridInput.in', 'w') as file:
              file.writelines(data)
            Grid.quit()
            Grid.destroy()
          else:
# Clear input
            clear_text()
            Grid.title('Grid Search variable ' + str(VarCount+1))


# Check previously entered variables are correct
    fail = False
    if not VarNum.get().isnumeric():
      messagebox.showerror('Error','Number of variables must be an integer')
      fail = True
    if not VarCells.get().isnumeric():
      messagebox.showerror('Error','Number of cells kept must be an integer')
      fail = True
    if not VarLvls.get().isnumeric():
      messagebox.showerror('Error','Number of levels searched must be an integer')
      fail = True

# Ask for variable information
    if not fail:
# Initial frame
      Grid = tk.Toplevel()
      Grid.title('Grid Search variable ' + str(VarCount+1))

# Labels
      tk.Label(Grid, text='Variable Name', font=('Roman',12)).grid(row=0,sticky='w')
      tk.Label(Grid, text='Variable Minimum Value', font=('Roman',12)).grid(row=1,sticky='w')
      tk.Label(Grid, text='Variable Maximum Value', font=('Roman',12)).grid(row=2,sticky='w')
      tk.Label(Grid, text='Variable Number of Models', font=('Roman',12)).grid(row=3,sticky='w')

# Entries
      VarName = tk.Entry(Grid)
      VarMin = tk.Entry(Grid)
      VarMax = tk.Entry(Grid)
      VarMods = tk.Entry(Grid)

# Entry locations
      VarName.grid(row=0, column=2)
      VarMin.grid(row=1, column=2)
      VarMax.grid(row=2, column=2)
      VarMods.grid(row=3, column=2)

# Save button
      btnStart=tk.Button(Grid,text='Save',height=1, width=4,relief=RAISED, justify=CENTER,command=lambda: [Count(), Save()], font=('Roman',12)).grid(row=4, column=2)

      Grid.mainloop()


# Make Frame
  GridRoot = tk.Toplevel()
  GridRoot.title('Grid Search input file')

# Labels
  tk.Label(GridRoot, text='Number of variables', font=('Roman',12)).grid(row=0,sticky='w')
  tk.Label(GridRoot, text='Number of best cells to keep', font=('Roman',12)).grid(row=1,sticky='w')
  tk.Label(GridRoot, text='Number of levels to search', font=('Roman',12)).grid(row=2,sticky='w')

# Entries
  VarNum = tk.Entry(GridRoot)
  VarCells = tk.Entry(GridRoot)
  VarLvls = tk.Entry(GridRoot)

# Entry locations
  VarNum.grid(row=0, column=2)
  VarCells.grid(row=1, column=2)
  VarLvls.grid(row=2, column=2)

# Enter button
  btnStart=tk.Button(GridRoot,text='Enter',height=1, width=5,relief=RAISED, justify=CENTER,command=GridInput2, font=('Roman',12)).grid(row=3, column=2)

  GridRoot.mainloop()


# List Input button
def StartList():


  def GMInfo():
    string = ('Geometric Mean scoring option requires additional input:\n'
              '1) the number of misfit scores you wish to use (2-5)\n'
              '2) a list of the misfit scores selecting from GV, SSR, SD, FSR, SA (case sensitive)\n'
              'For example: 2,GV,SD \n'
              'Geometric mean is defined as:\n'
              '(X1 * X2 *...* Xn)^(1/n)')

    messagebox.showinfo('Geometric Mean',string)


  def StartList2():

    global MisRank, GMcla


    def NotReady():
      string = ('Not all required inputs completed:\n'
                'You must input at least the number of MPI workers and Iterations.')
      messagebox.showerror('Error',string)


    def Start():
      global GMcla

      MisRankOpt = ['GV','SSR','SD','FSR','SC','SA','GM']

# Check required inputs are complete
      if len(NP.get()) != 0  and  len(ITER.get())!= 0:
        fail = False

# Check integer & float values are correctly input
# MPI - iter
        if not NP.get().isnumeric():
          messagebox.showerror('Error','Number of MPI workers should be an integer')
          fail = True
# Iter
        if not ITER.get().isnumeric():
          messagebox.showerror('Error','Number of Iterations should be an integer')
          fail = True
# MC
        if len(MC1.get()) != 0  or  len(MC2.get()) != 0  or  len(MC3.get()) != 0:
# MC type - against list
          if MC1.get() not in MisRankOpt:
            string = ('Misfit Convergence misfit is not in OrbScore options, choose from (case sensitive):\n'
                      'GV, SSR, SD, FSR, SC, SA, GM')
            messagebox.showerror('Error',string)
            fail = True
# MC start iter - integer
          if not MC2.get().isnumeric():
            messagebox.showerror('Error','Misfit Convergence start iteration should be integer')
            fail = True
# MC conv rate - float
          try:
            float(MC3.get())
            if not ((float(MC3.get()) > 0) and (float(MC3.get()) < 1)) :
              messagebox.showerror('Error','Misfit Convergence range should be decimal from 0 to 1')
              fail = True
          except:
            messagebox.showerror('Error','Misfit Convergence range should be decimal from 0 to 1')
            fail = True
# Decide whether or not ShellSet should launch
        if not fail:
          import os
# Generate intitialisation line
          if var.get() ==1:
            Version = 'ShellSet.exe'
          if var.get() ==2:
            Version = 'ShellSet_OPT.exe'
          if var.get() ==3:
            Version = 'ShellSet_DB.exe'

          cmd = 'mpirun -np ' + NP.get() + ' ./' + Version + ' -InOpt List -Iter ' + ITER.get()
          if len(DIR.get()) != 0:
            cmd += ' -Dir ' + DIR.get()
          if len(MC1.get()) != 0  and  len(MC2.get()) != 0  and  len(MC3.get()) != 0:
            cmd += ' -MC ' + MC1.get().upper() + ',' + MC2.get() + ',' + MC3.get()
          if VLog.get():
            cmd += ' -V '
          if AEFLog.get():
            cmd += ' -AEF '
# Add any -GM definition
          if MisRank=='GM':
            cmd += ' -GM ' + GMcla.upper()

          print('Initialising ShellSet with following command line:')
          print(cmd)

          root.quit()
          root.destroy()
          os.system(cmd)

      else:
        NotReady()


# Make Frame
    List = tk.Toplevel()
    List.title('Listed models user options')

    VLog   = tk.BooleanVar(List)
    AEFLog = tk.BooleanVar(List)

# Labels
    tk.Label(List, text='Required List Search input:', font=('Roman',15)).grid(row=0,column=1)
    tk.Label(List, text='MPI workers', font=('Roman',12)).grid(row=1,sticky='w')
    tk.Label(List, text='Main loop max Iterations', font=('Roman',12)).grid(row=2,sticky='w')

    tk.Label(List, text='Suggested (Optional) input:', font=('Roman',15)).grid(row=3,column=1)
    tk.Label(List, text='Misfit Convergence', font=('Roman',12)).grid(row=4,sticky='w')
    tk.Label(List, text='Working Directory', font=('Roman',12)).grid(row=5,sticky='w')
    tk.Label(List, text='Verbose', font=('Roman',12)).grid(row=6,sticky='w')
    tk.Label(List, text='All Errors Fatal', font=('Roman',12)).grid(row=7,sticky='w')

# Entries
    NP   = tk.Entry(List, font=('Roman',12))
    ITER = tk.Entry(List, font=('Roman',12))
    MC1  = tk.Entry(List, font=('Roman',12))
    MC2  = tk.Entry(List, font=('Roman',12))
    MC3  = tk.Entry(List, font=('Roman',12))
    DIR  = tk.Entry(List, font=('Roman',12))

    if ScoreVar.get() ==1:
      MisRank = ''
    if ScoreVar.get() ==2:
      MisRank = 'GV'
    if ScoreVar.get() ==3:
      MisRank = 'SSR'
    if ScoreVar.get() ==4:
      MisRank = 'SD'
    if ScoreVar.get() ==5:
      MisRank = 'FSR'
    if ScoreVar.get() ==6:
      MisRank = 'SA'
    if ScoreVar.get() ==7:
      MisRank = 'SC'
    if ScoreVar.get() ==8:
      MisRank = 'GM'

    MC1.insert(0, MisRank)
    GMcla = GM.get()

# Check boxes
    V    = tk.Checkbutton(List, text='On/Off',variable=VLog, onvalue=True, offvalue=False, font=('Roman',12))
    AEF  = tk.Checkbutton(List, text='On/Off',variable=AEFLog, onvalue=True, offvalue=False, font=('Roman',12))

# Entry locations
    NP.grid(row=1, column=1)
    ITER.grid(row=2, column=1)
    MC1.grid(row=4, column=1)
    MC2.grid(row=4, column=2)
    MC3.grid(row=4, column=3)
    DIR.grid(row=5, column=1)
    V.grid(row=6, column=1)
    AEF.grid(row=7, column=1)

# Start button
    btnStart=tk.Button(List,text='Start',command=Start, font=('Roman',12)).grid(row=6, column=3)

# Close OrbScore option frame
    ScoreOpt.quit()
    ScoreOpt.destroy()

    List.mainloop()


# Make Frame
  ScoreOpt = tk.Toplevel()
  ScoreOpt.title('Misfit Convergence (optional)')

  GVLog  = tk.BooleanVar(ScoreOpt)
  SSRLog = tk.BooleanVar(ScoreOpt)
  SDLog  = tk.BooleanVar(ScoreOpt)
  FSRLog = tk.BooleanVar(ScoreOpt)
  SCLog  = tk.BooleanVar(ScoreOpt)
  SALog  = tk.BooleanVar(ScoreOpt)

# Labels
  tk.Label(ScoreOpt, text='Misfit Ranker for Misfit Convergence', font=('Roman',15)).grid(row=0,rowspan=3)

  ScoreVar = tk.IntVar()
# Radio buttons
  R1 = tk.Radiobutton(ScoreOpt, text='None', variable=ScoreVar, value=1, font=('Roman',12))
  R2 = tk.Radiobutton(ScoreOpt, text='GV  = Geodetic Velocity (misfit)', variable=ScoreVar, value=2, font=('Roman',12))
  R3 = tk.Radiobutton(ScoreOpt, text='SSR = Seafloor Spreading Rate (misfit)', variable=ScoreVar, value=3, font=('Roman',12))
  R4 = tk.Radiobutton(ScoreOpt, text='SD  = Stress Direction (misfit)', variable=ScoreVar, value=4, font=('Roman',12))
  R5 = tk.Radiobutton(ScoreOpt, text='FSR = Fault Slip Rate (misfit)', variable=ScoreVar, value=5, font=('Roman',12))
  R6 = tk.Radiobutton(ScoreOpt, text='SA  = Seismic Anisotropy (misfit)', variable=ScoreVar, value=6, font=('Roman',12))
  R7 = tk.Radiobutton(ScoreOpt, text='SC  = Seismic Correlation (score)', variable=ScoreVar, value=7, font=('Roman',12))
  R8 = tk.Radiobutton(ScoreOpt, text='GM  = Geometric Mean (misfit)', variable=ScoreVar, value=8, font=('Roman',12))

  ScoreVar.set(1)

# GM entry & info button
  GM = tk.Entry(ScoreOpt, font=('Roman',12))
  btnGM=tk.Button(ScoreOpt,text='?',command=GMInfo, font=('Roman', 8)).grid(row=10, column=2)

# Entry locations
  R1.grid(row=3, column=0,sticky='w')
  R2.grid(row=4, column=0,sticky='w')
  R3.grid(row=5, column=0,sticky='w')
  R4.grid(row=6, column=0,sticky='w')
  R5.grid(row=7, column=0,sticky='w')
  R6.grid(row=8, column=0,sticky='w')
  R7.grid(row=9, column=0,sticky='w')
  R8.grid(row=10, column=0,sticky='w')
  GM.grid(row=10, column=1,sticky='w')

# Enter button
  btnStart=tk.Button(ScoreOpt,text='Enter',height=1, width=5,relief=RAISED, justify=CENTER,command=StartList2, font=('Roman',12)).grid(row=11, column=1)

  ScoreOpt.mainloop()


# Grid Search button
def StartGrid():

  def GMInfo():
    string = ('Geometric Mean scoring option requires additional input:\n'
              '1) the number of misfit scores you wish to use (2-5)\n'
              '2) a list of the misfit scores selecting from GV, SSR, SD, FSR, SA (case sensitive)\n'
              'For example: 2,GV,SD \n'
              'Geometric mean is defined as:\n'
              '(X1 * X2 *...* Xn)^(1/n)')

    messagebox.showinfo('Geometric Mean',string)

  def StartGrid2():

    global MisRank, GMcla

    def NotReady():
      string = ('Not all required inputs completed:\n'
                'You must input at least the number of MPI workers, Iterations and one of Misfit Ranker or '
                'Misfit Convergence.')
      messagebox.showerror('Error',string)

    def MRInfo():
      string = ('Select a Misfit Ranker option from:\n'
                'GV, SSR, SD, FSR, SA, SC, GM \n'
                'in the 1st box. \n'
                'The GM option requires additional input in the 2nd box:\n'
                '1) the number of misfit scores you wish to use (2-5)\n'
                '2) a list of the misfit scores selecting from GV, SSR, SD, FSR, SA (case sensitive)\n'
                'For example: 2,GV,SD \n'
                'Geometric mean is defined as:\n'
                '(X1 * X2 *...* Xn)^(1/n)')

      messagebox.showinfo('Misfit Ranker',string)

    def Start():

      global GMcla

      MisRankOpt = ['GV','SSR','SD','FSR','SC','SA','GM']

# Check required inputs are complete
      if(len(NP.get())!=0 and len(ITER.get())!=0 and
      (len(MR.get())!=0 or (len(MC1.get())!=0 and len(MC2.get())!=0 and len(MC3.get())!=0)) ):
        fail = False

# Check integer & float values are correctly input
# MPI - iter
        if not NP.get().isnumeric():
          messagebox.showerror('Error','Number of MPI workers should be an integer')
          fail = True
# Iter
        if not ITER.get().isnumeric():
          messagebox.showerror('Error','Number of Iterations should be an integer')
          fail = True
# MC
        if len(MC1.get()) != 0  and  len(MC2.get()) != 0  and  len(MC3.get()) != 0:
# MC type - against list
          if MC1.get() not in MisRankOpt:
            string = ('Misfit Convergence misfit is not in OrbScore options, choose from (case sensitive):\n'
                      'GV, SSR, SD, FSR, SC, SA, GM')
            messagebox.showerror('Error',string)
            fail = True
# MC start iter - integer
          if not MC2.get().isnumeric():
            messagebox.showerror('Error','Misfit Convergence start iteration should be integer')
            fail = True
# MC conv rate - float
          try:
            float(MC3.get())
            if not ((float(MC3.get()) > 0) and (float(MC3.get()) < 1)) :
              messagebox.showerror('Error','Misfit Convergence range should be decimal from 0 to 1')
              fail = True
          except:
            messagebox.showerror('Error','Misfit Convergence range should be decimal from 0 to 1')
            fail = True
# MR - against list
        else:
          if MR.get() not in MisRankOpt:
            string = ('Misfit Ranker is not in OrbScore options, choose from (case sensitive):\n'
                      'GV, SSR, SD, FSR, SC, SA, GM')
            messagebox.showerror('Error',string)
            fail = True
# ML - float
        if len(ML.get()) != 0 :
          try :
            float(ML.get())
          except :
            messagebox.showerror('Error','Misfit Limit should be numeric')
            fail = True
# KL - float
        if len(KL.get()) != 0 :
          try :
            float(KL.get())
          except :
            messagebox.showerror('Error','Kill Limit should be numeric')
            fail = True

# Decide whether or not ShellSet should launch
        if not fail:
          import os
# Generate intitialisation line
          if var.get() ==1:
            Version = 'ShellSet.exe'
          if var.get() ==2:
            Version = 'ShellSet_OPT.exe'
          if var.get() ==3:
            Version = 'ShellSet_DB.exe'

          cmd = 'mpirun -np ' + NP.get() + ' ./' + Version + ' -InOpt Grid -Iter ' + ITER.get()
          if len(DIR.get()) != 0:
            cmd += ' -Dir ' + DIR.get()
          if len(MC1.get()) != 0  and  len(MC2.get()) != 0  and  len(MC3.get()) != 0:
            cmd += ' -MC ' + MC1.get().upper() + ',' + MC2.get() + ',' + MC3.get()
          if len(MR.get()) != 0  and  MisRank != 'GM'  and   (len(MC1.get()) == 0  or  len(MC2.get()) == 0  or  len(MC3.get()) == 0):
            if MR.get().upper() == 'GM':
              cmd += ' -GM ' + GM.get()
            else:
              cmd += ' -MR ' + MR.get().upper()
          if len(ML.get()) != 0:
            cmd += ' -ML ' + MLStr
          if len(KL.get()) != 0:
            cmd += ' -KL ' + KLStr
          if VLog.get():
            cmd += ' -V '
          if AEFLog.get():
            cmd += ' -AEF '
# Add any -GM definition
          if MisRank=='GM':
            cmd += ' -GM ' + GMcla.upper()

          print('Initialising ShellSet with following command line:')
          print(cmd)

          root.quit()
          root.destroy()
          os.system(cmd)

      else:
        NotReady()


# Make Frame
    Grid = tk.Toplevel()
    Grid.title('Grid Search user options')

    VLog    = tk.BooleanVar(Grid)
    AEFLog  = tk.BooleanVar(Grid)

# Labels
    tk.Label(Grid, text='Required Grid Search input:', font=('Roman',15)).grid(row=0,column=1)
    tk.Label(Grid, text='MPI workers', font=('Roman',12)).grid(row=1,sticky='w')
    tk.Label(Grid, text='Main loop max Iterations', font=('Roman',12)).grid(row=2,sticky='w')

    tk.Label(Grid, text='Must use 1 of the following 2 options:', font=('Roman',15)).grid(row=3,column=1)
    tk.Label(Grid, text='Misfit Ranker (select best models)', font=('Roman',12)).grid(row=4,sticky='w')
    tk.Label(Grid, text='Misfit Convergence', font=('Roman',12)).grid(row=5,sticky='w')

    tk.Label(Grid, text='Suggested (Optional) input:', font=('Roman',15)).grid(row=6,column=1)
    tk.Label(Grid, text='Working Directory', font=('Roman',12)).grid(row=7,sticky='w')
    tk.Label(Grid, text='Misfit Limit', font=('Roman',12)).grid(row=8,sticky='w')
    tk.Label(Grid, text='Kill Limit', font=('Roman',12)).grid(row=9,sticky='w')
    tk.Label(Grid, text='Verbose', font=('Roman',12)).grid(row=10,sticky='w')
    tk.Label(Grid, text='All Errors Fatal', font=('Roman',12)).grid(row=11,sticky='w')

# Entries
    NP   = tk.Entry(Grid, font=('Roman',12))
    ITER = tk.Entry(Grid, font=('Roman',12))
    MR   = tk.Entry(Grid, font=('Roman',12))
    GM   = tk.Entry(Grid, font=('Roman',12)) # GM entry
    btnGM = tk.Button(Grid,text='?',command=MRInfo, font=('Roman', 8)).grid(row=4, column=3) # info button
    MC1  = tk.Entry(Grid, font=('Roman',12))
    MC2  = tk.Entry(Grid, font=('Roman',12))
    MC3  = tk.Entry(Grid, font=('Roman',12))
    DIR  = tk.Entry(Grid, font=('Roman',12))
    ML   = tk.Entry(Grid, font=('Roman',12))
    KL   = tk.Entry(Grid, font=('Roman',12))

    if ScoreVar.get() ==1:
      MisRank = ''
    if ScoreVar.get() ==2:
      MisRank = 'GV'
    if ScoreVar.get() ==3:
      MisRank = 'SSR'
    if ScoreVar.get() ==4:
      MisRank = 'SD'
    if ScoreVar.get() ==5:
      MisRank = 'FSR'
    if ScoreVar.get() ==6:
      MisRank = 'SA'
    if ScoreVar.get() ==7:
      MisRank = 'SC'
    if ScoreVar.get() ==8:
      MisRank = 'GM'

    MR.insert(0, MisRank)
    MC1.insert(0, MisRank)
    GMcla = GM_MC.get()

# Check boxes
    V    = tk.Checkbutton(Grid, text='On/Off',variable=VLog, onvalue=True, offvalue=False, font=('Roman',12))
    AEF  = tk.Checkbutton(Grid, text='On/Off',variable=AEFLog, onvalue=True, offvalue=False, font=('Roman',12))

# Entry locations
    NP.grid(row=1, column=1)
    ITER.grid(row=2, column=1)
    MR.grid(row=4, column=1)
    GM.grid(row=4, column=2)
    MC1.grid(row=5, column=1)
    MC2.grid(row=5, column=2)
    MC3.grid(row=5, column=3)
    DIR.grid(row=7, column=1)
    ML.grid(row=8, column=1)
    KL.grid(row=9, column=1)
    V.grid(row=10, column=1)
    AEF.grid(row=11, column=1)

# Start button
    btnStart=tk.Button(Grid,text='Start',command=Start, font=('Roman',12),).grid(row=9, column=3)

# Close OrbScore option frame
    ScoreOpt.quit()
    ScoreOpt.destroy()

    Grid.mainloop()


# Make Frame
  ScoreOpt = tk.Toplevel()
  ScoreOpt.title('Misfit Convergence (optional)')

  GVLog  = tk.BooleanVar(ScoreOpt)
  SSRLog = tk.BooleanVar(ScoreOpt)
  SDLog  = tk.BooleanVar(ScoreOpt)
  FSRLog = tk.BooleanVar(ScoreOpt)
  SCLog  = tk.BooleanVar(ScoreOpt)
  SALog  = tk.BooleanVar(ScoreOpt)

# Labels
  tk.Label(ScoreOpt, text='Misfit Ranker for Misfit Convergence', font=('Roman',15)).grid(row=0,rowspan=3)

  ScoreVar = tk.IntVar()
# Radio buttons
  R1 = tk.Radiobutton(ScoreOpt, text='None', variable=ScoreVar, value=1, font=('Roman',12))
  R2 = tk.Radiobutton(ScoreOpt, text='GV  = Geodetic Velocity (misfit)', variable=ScoreVar, value=2, font=('Roman',12))
  R3 = tk.Radiobutton(ScoreOpt, text='SSR = Seafloor Spreading Rate (misfit)', variable=ScoreVar, value=3, font=('Roman',12))
  R4 = tk.Radiobutton(ScoreOpt, text='SD  = Stress Direction (misfit)', variable=ScoreVar, value=4, font=('Roman',12))
  R5 = tk.Radiobutton(ScoreOpt, text='FSR = Fault Slip Rate (misfit)', variable=ScoreVar, value=5, font=('Roman',12))
  R6 = tk.Radiobutton(ScoreOpt, text='SA  = Seismic Anisotropy (misfit)', variable=ScoreVar, value=6, font=('Roman',12))
  R7 = tk.Radiobutton(ScoreOpt, text='SC  = Seismic Correlation (score)', variable=ScoreVar, value=7, font=('Roman',12))
  R8 = tk.Radiobutton(ScoreOpt, text='GM  = Geometric Mean (misfit)', variable=ScoreVar, value=8, font=('Roman',12))

  ScoreVar.set(1)

# GM entry & info button
  GM_MC = tk.Entry(ScoreOpt, font=('Roman',12))
  btnGM=tk.Button(ScoreOpt,text='?',command=GMInfo, font=('Roman', 8)).grid(row=10, column=2)

# Entry locations
  R1.grid(row=3, column=0,sticky='w')
  R2.grid(row=4, column=0,sticky='w')
  R3.grid(row=5, column=0,sticky='w')
  R4.grid(row=6, column=0,sticky='w')
  R5.grid(row=7, column=0,sticky='w')
  R6.grid(row=8, column=0,sticky='w')
  R7.grid(row=9, column=0,sticky='w')
  R8.grid(row=10, column=0,sticky='w')
  GM_MC.grid(row=10, column=1,sticky='w')

# Enter button
  btnStart=tk.Button(ScoreOpt,text='Enter',height=1, width=5,relief=RAISED, justify=CENTER,command=StartGrid2, font=('Roman',12)).grid(row=11, column=1)

  ScoreOpt.mainloop()


# Close root window
def main_closing():

  def on_closing():
    close.quit()
    close.destroy()
    root.quit()
    root.destroy()


  if messagebox.askokcancel('Quit', 'Exit ShellSet GUI?'):
    close = Tk()
    close.geometry('650x120')
    close.title('ShellSet Copyright (C) 2023')
    close.protocol('WM_DELETE_WINDOW', on_closing)

    w = Label(close, text='Dr Jon Bryan May (INGV), Prof. Peter Bird (UCLA), Dr Michele Carafa (INGV).', font=('Roman',12))
    w.pack()
    w = Label(close, text='This program comes with ABSOLUTELY NO WARRANTY.', font=('Roman',12))
    w.pack(side=TOP, anchor=NW)
    w = Label(close, text='This is free software, you are welcome to redistribute it under conditions set out', font=('Roman',12))
    w.pack(side=TOP, anchor=NW)
    w = Label(close, text='in GNU General Public License version 3, or any later version.', font=('Roman',12))
    w.pack(side=TOP, anchor=NW)
    w = Label(close,text='See the included license file, or https://www.gnu.org/licenses/, for details.', font=('Roman',12))
    w.pack(side=TOP, anchor=NW)

    close.mainloop()


# Root Window
root=tk.Tk()
root.title('ShellSet')
root.geometry('500x270')
root.protocol('WM_DELETE_WINDOW', main_closing)

label = tk.Label(root, text='Select ShellSet version:', font=('Roman',15)).place(x=100,y=10)
var = tk.IntVar()
if Default == 'Standard':
  var.set(1)
elif Default == 'Optimal':
  var.set(2)
elif Default == 'Debug':
  var.set(3)

R1 = tk.Radiobutton(root, text='Standard', variable=var, value=1, font=('Roman',12)).place(x=100,y=35)
R2 = tk.Radiobutton(root, text='Optimal', variable=var, value=2, font=('Roman',12)).place(x=210,y=35)
R3 = tk.Radiobutton(root, text='Debug', variable=var, value=3, font=('Roman',12)).place(x=310,y=35)

label = tk.Label(root, text='Update/Create input files?:', font=('Roman',15)).place(x=100,y=70)
btnIn=tk.Button(root, text='Parameter input', relief=RAISED, justify=CENTER, command=ParamInput, font=('Roman',12), height=1, width=15).place(x=100,y=100)
btnIn=tk.Button(root, text='Input files list', relief=RAISED, justify=CENTER, command=InputFiles, font=('Roman',12), height=1, width=15).place(x=300,y=100)
btnGrIn=tk.Button(root, text='List input', relief=RAISED, justify=CENTER, command=ListInput, font=('Roman',12), height=1, width=10).place(x=100,y=140)
btnGrIn=tk.Button(root, text='Grid input', relief=RAISED, justify=CENTER, command=GridInput, font=('Roman',12), height=1, width=10).place(x=250,y=140)

label = tk.Label(root, text='Select desired model input type:', font=('Roman',15)).place(x=100,y=180)
btnLi=tk.Button(root, text='List', relief=RAISED, justify=CENTER, command=StartList, font=('Roman',12), height=1, width=4).place(x=100,y=210)
btnGr=tk.Button(root, text='Grid', relief=RAISED, justify=CENTER, command=StartGrid, font=('Roman',12), height=1, width=4).place(x=190,y=210)

if Default != 'Standard' and Default != 'Optimal' and Default != 'Debug':
  string = ('Default ShellSet version set in ShellSet_GUI.py not recognised. \n '
            'You will need to select one for this simulation. \n'
            'You are advised to update the option in ShellSet_GUI.py.')
  messagebox.showerror('Error',string)

root.attributes('-topmost',True)
root.mainloop()