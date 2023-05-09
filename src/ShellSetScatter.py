import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import pyplot as plt
import re
import os

# *************************************************************************
# ********************************* Modify ********************************

TestDir   = "EXAMPLES/ListEx1" # Name of directory where test is stored
ModInFile = "ListInput.in"     # Name of Grid or List Input file - for variable names
ResFile   = "Models.txt"       # Name of Results file - for location and results data
OutPath   = "Plots"            # Name of output directory - will be created in TestDir if it does not exist
OutType   = "png"              # Filetype that generated graphs will be saved as (eps, pdf, pgf, png, ps, raw, rgba, svg, svgz)
ErrLvl    = "error"            # matplotlib error level, options: "notset", "debug", "info", "warning", "error", "critical"

# *************************************************************************
# *************************************************************************
# ********************* You should not need to modify *********************

matplotlib.set_loglevel(ErrLvl)

ModInFile = TestDir+"/"+ModInFile
ResFile   = TestDir+"/"+ResFile
OutPath   = TestDir+"/"+OutPath+"/"

# Read model input type from ResFile
with open(ResFile, "r") as f:
  content = f.readline()
  Scores = f.readline()

content = content.split(" ")
for i in content:
  if i == "-InOpt":
    Input = content[content.index(i)+1]

# Find options for score plotting
if Input == "List":
  Scores = Scores.strip("\n").split(" ")[4:]
elif Input == "Grid":
  Scores = Scores.strip("\n").split(" ")[6:]

for i in range(0, len(Scores)-1):
  Scores[i] = Scores[i].strip(",")

print("Select one of "+" ".join(Scores)+" to plot:")
ScrPlt = input()

if ScrPlt not in Scores:
  exit("Wrong score option entered")

ScrPltNum = Scores.index(ScrPlt)+1

# Create Plots directory if it does not exits
if not os.path.exists(OutPath):
  os.mkdir(OutPath)



# *************************************************************************
# ******************************* List Input ******************************
if Input == "List":

# Read information from ModInFile
  with open(ModInFile, "r") as f:
    content = f.readlines()

  NumVar = int(re.split(" |,",content[0])[0])  # Number of Variables
  NumMod = int(re.split(" |,",content[0])[1])  # Number of Models

# ****************************** 1D Test ******************************
  if NumVar == 1:

    x1Label = content[1].split(" ")[0]          # X1 variable name
    data = np.genfromtxt(ResFile,dtype="str",skip_header=2,usecols=(2,2+ScrPltNum)) # Read as string to handle "Failed"

# Set "Failed" models to an unobtainable -1 misfit
    data[:,1][data[:,1] == "Failed"]     = "-1"
    data[:,1][data[:,1] == "to"]         = "-1"
    data[:,1][data[:,1] == "converge"]   = "-1"
    data[:,1][data[:,1] == "in"]         = "-1"
    data[:,1][data[:,1] == "required"]   = "-1"
    data[:,1][data[:,1] == "MKL"]        = "-1"
    data[:,1][data[:,1] == "iterations"] = "-1"
    data[:,1][data[:,1] == "input"]      = "-1"
    data[:,1][data[:,1] == "variable"]   = "-1"
    data[:,1][data[:,1] == "check"]      = "-1"
    data[:,1][data[:,1] == "-"]          = "-1"
    data[:,1][data[:,1] == "values"]     = "-1"

# Convert to float
    X1in    = data[:,0].astype(float)
    Mis     = data[:,1].astype(float)
    Mis[Mis == -1] = max(Mis) + max(Mis)/2       # Set "Failed" to be 50% higher than other maximum for plot

# ********** Plot **********
    x1  = X1in
    mis = Mis

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cmhot = plt.cm.get_cmap("gnuplot")
    l = ax.scatter(x1, mis, s=20, marker="o", c=mis, cmap=cmhot)
    fig.colorbar(l,label=ScrPlt+" Score")
    plt.title("All models scored using "+ScrPlt)
    ax.set_xlabel(x1Label)
    ax.set_ylabel("Misfit")
    # ax.set_xlim([min, max])
    # ax.set_ylim([min, max])
    fig.savefig(OutPath+ScrPlt+"_AllModels."+OutType)
    plt.show()


# ****************************** 2D Test ******************************
  elif NumVar == 2:

    x1Label = content[1].split(" ")[0]          # X1 variable name
    x2Label = content[2].split(" ")[0]          # X2 variable name
    data = np.genfromtxt(ResFile,dtype="str",skip_header=2,usecols=(2,3,3+ScrPltNum)) # Read as string to handle "Failed"

# Set "Failed" models to an unobtainable -1 misfit
    data[:,2][data[:,2] == "Failed"]     = "-1"
    data[:,2][data[:,2] == "to"]         = "-1"
    data[:,2][data[:,2] == "converge"]   = "-1"
    data[:,2][data[:,2] == "in"]         = "-1"
    data[:,2][data[:,2] == "required"]   = "-1"
    data[:,2][data[:,2] == "MKL"]        = "-1"
    data[:,2][data[:,2] == "iterations"] = "-1"
    data[:,2][data[:,2] == "input"]      = "-1"
    data[:,2][data[:,2] == "variable"]   = "-1"
    data[:,2][data[:,2] == "check"]      = "-1"
    data[:,2][data[:,2] == "-"]          = "-1"
    data[:,2][data[:,2] == "values"]     = "-1"

# Convert to float
    X1in    = data[:,0].astype(float)
    X2in    = data[:,1].astype(float)
    Mis     = data[:,2].astype(float)
    Mis[Mis == -1] = max(Mis) + max(Mis)/2       # Set "Failed" to be 50% higher than other maximum for plot

# ********** Plot **********
    x1  = X1in
    x2  = X2in
    mis = Mis

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cmhot = plt.cm.get_cmap("gnuplot")
    s = 10
    l = ax.scatter(x1, x2, s=20, marker="o", c=mis, cmap=cmhot)
    fig.colorbar(l,label=ScrPlt+" Score")
    plt.title("All models scored using "+ScrPlt)
    ax.set_xlabel(x1Label)
    ax.set_ylabel(x2Label)
    # ax.set_xlim([min, max])
    # ax.set_ylim([min, max])
    fig.savefig(OutPath+ScrPlt+"_AllModels."+OutType)
    plt.show()


# ****************************** 3D Test ******************************
  elif NumVar == 3:

    x1Label = content[1].split(" ")[0]          # X1 variable name
    x2Label = content[2].split(" ")[0]          # X2 variable name
    x3Label = content[3].split(" ")[0]          # X3 variable name
    data = np.genfromtxt(ResFile,dtype="str",skip_header=2,usecols=(2,3,4,4+ScrPltNum)) # Read as string to handle "Failed"

# Set "Failed" models to an unobtainable -1 misfit
    data[:,3][data[:,3] == "Failed"]     = "-1"
    data[:,3][data[:,3] == "to"]         = "-1"
    data[:,3][data[:,3] == "converge"]   = "-1"
    data[:,3][data[:,3] == "in"]         = "-1"
    data[:,3][data[:,3] == "required"]   = "-1"
    data[:,3][data[:,3] == "MKL"]        = "-1"
    data[:,3][data[:,3] == "iterations"] = "-1"
    data[:,3][data[:,3] == "input"]      = "-1"
    data[:,3][data[:,3] == "variable"]   = "-1"
    data[:,3][data[:,3] == "check"]      = "-1"
    data[:,3][data[:,3] == "-"]          = "-1"
    data[:,3][data[:,3] == "values"]     = "-1"

# Convert to float
    X1in    = data[:,0].astype(float)
    X2in    = data[:,1].astype(float)
    X3in    = data[:,2].astype(float)
    Mis     = data[:,3].astype(float)
    Mis[Mis == -1] = max(Mis) + max(Mis)/2       # Set "Failed" to be 50% higher than other maximum for plot

# ********** Plot **********
    x1  = X1in
    x2  = X2in
    x3  = X3in
    mis = Mis
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    cmhot = plt.cm.get_cmap("gnuplot")
    l = ax.scatter(x1, x2, x3, s=20, c=mis, marker="o", cmap=cmhot)
    fig.colorbar(l,label=ScrPlt+" Score")
    plt.title("All models scored using "+ScrPlt)
    ax.set_xlabel(x1Label)
    ax.set_ylabel(x2Label)
    ax.set_zlabel(x3Label)
    # ax.set_xlim([min, max])
    # ax.set_ylim([min, max])
    # ax.set_zlim([min, max])
    fig.savefig(OutPath+ScrPlt+"_AllModels."+OutType)
    plt.show()



# *************************************************************************
# ****************************** Grid Search ******************************
if Input == "Grid":

# Read information from ModInFile
  with open(ModInFile, "r") as f:
    content = f.readlines()

  NumVar = int(re.split(" |,",content[0])[0])  # Number of Variables
  KCell  = int(re.split(" |,",content[0])[1])  # Number of Kept Cells
  NumLvl = int(re.split(" |,",content[0])[2])  # Number of Levels

# ****************************** 1D Search ******************************
  if NumVar == 1:
    print("Plotting 1D grid results")
    print("X=Variable Y=Misfit")

    AdaptiveLayers = [KCell]*(NumLvl-1)
    AdaptiveLayers = [1] + AdaptiveLayers

    x1Label = content[2].split(" ")[0]          # X1 variable name
    x1Mods  = int(content[5].split(" ")[0])     # Number of models in X1

# ********** Repeated Models? **********

    if x1Mods % 2 != 0 and ScrPltNum != 1:
      print("Selected plot option was not the misfit ranker - copying required data")

      with open(ResFile, "r") as f:
        Rslts = f.readlines()[2:]

      for i in range(1,((NumLvl-1) * KCell)+1): # Update every repeated model
        for line in Rslts:
          if line.split()[1] == "0":
            index1 = [x for x in range(len(Rslts)) if line.split()[1] in Rslts[x].split()[1]] # find index of repeated model
            for line2 in Rslts[:index1[0]]:
              if line2.split()[4] == line.split()[4]:
                index = [x for x in range(len(Rslts)) if line2.split()[4] in Rslts[x].split()[4]] # find index of original model
                Rslts[index1[0]] = Rslts[index[0]] # update repeat with original results
      data = np.genfromtxt(Rslts,dtype="str",usecols=(4,4+ScrPltNum))

    else:
      data = np.genfromtxt(ResFile,dtype="str",skip_header=2,usecols=(4,4+ScrPltNum))

# Set "Failed" models to an unobtainable -1 misfit
    data[:,1][data[:,1] == "Failed"]     = "-1"
    data[:,1][data[:,1] == "to"]         = "-1"
    data[:,1][data[:,1] == "converge"]   = "-1"
    data[:,1][data[:,1] == "in"]         = "-1"
    data[:,1][data[:,1] == "required"]   = "-1"
    data[:,1][data[:,1] == "MKL"]        = "-1"
    data[:,1][data[:,1] == "iterations"] = "-1"
    data[:,1][data[:,1] == "input"]      = "-1"
    data[:,1][data[:,1] == "variable"]   = "-1"
    data[:,1][data[:,1] == "check"]      = "-1"
    data[:,1][data[:,1] == "-"]          = "-1"
    data[:,1][data[:,1] == "values"]     = "-1"

    NumMods = x1Mods

# Convert data to float
    X1in    = data[:,0].astype(float)
    Mis     = data[:,1].astype(float)
    Mis[Mis == -1] = max(Mis) + max(Mis)/2       # Set "Failed" to be 50% higher than other maximum for plot

# ********** Plot Levels **********
    for iPlot in range(1,NumLvl+1):
      x1  = X1in[sum(AdaptiveLayers[0:iPlot-1])*NumMods : sum(AdaptiveLayers[0:iPlot])*NumMods]
      mis = Mis[sum(AdaptiveLayers[0:iPlot-1])*NumMods : sum(AdaptiveLayers[0:iPlot])*NumMods]

      fig = plt.figure()
      ax = fig.add_subplot(111)
      cmhot = plt.cm.get_cmap("gnuplot")
      l = ax.scatter(x1, mis, s=20, marker="o", c=mis, cmap=cmhot)
      fig.colorbar(l,label=ScrPlt+" Score")
      plt.title("Level "+str(iPlot)+" - scored using "+ScrPlt)
      ax.set_xlabel(x1Label)
      ax.set_ylabel("Misfit")
      # ax.set_xlim([min, max])
      # ax.set_ylim([min, max])
      fig.savefig(OutPath+ScrPlt+"_Level" + str(iPlot) + "_1D."+OutType)

# ********** All models **********
    x1  = X1in
    mis = Mis

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cmhot = plt.cm.get_cmap("gnuplot")
    l = ax.scatter(x1, mis, s=20, marker="o", c=mis, cmap=cmhot)
    fig.colorbar(l,label=ScrPlt+" Score")
    plt.title("All models - scored using "+ScrPlt)
    ax.set_xlabel(x1Label)
    ax.set_ylabel("Misfit")
    # ax.set_xlim([min, max])
    # ax.set_ylim([min, max])
    fig.savefig(OutPath+ScrPlt+"_AllModels_1D."+OutType)
    plt.show()


# ****************************** 2D Search ******************************
  elif NumVar == 2:
    print("Plotting 2D grid results")
    print("X=Variable1 Y=Variable2 colour=Misfit")

    AdaptiveLayers = [KCell]*(NumLvl-1)
    AdaptiveLayers = [1] + AdaptiveLayers

    x1Label = content[2].split(" ")[0]          # X1 variable name
    x2Label = content[7].split(" ")[0]          # X2 variable name
    x1Mods  = int(content[5].split(" ")[0])     # Number of models in X1
    x2Mods  = int(content[10].split(" ")[0])    # Number of models in X2

# ********** Repeated Models? **********
    if x1Mods % 2 != 0 and x2Mods % 2 != 0 and ScrPltNum != 1:
      print("Selected plot option was not the misfit ranker - copying required data\n")

      with open(ResFile, "r") as f:
        Rslts = f.readlines()[2:]

      for i in range(1,((NumLvl-1) * KCell)+1): # Update every repeated model
        for line in Rslts:
          if line.split()[1] == "0":
            index1 = [x for x in range(len(Rslts)) if line.split()[1] in Rslts[x].split()[1]] # find index of repeated model
            for line2 in Rslts[:index1[0]]:
              if(line2.split()[4] == line.split()[4] and line2.split()[5] == line.split()[5]):
                index = [x for x in range(len(Rslts)) if (line2.split()[4] in Rslts[x].split()[4] and line2.split()[5] in Rslts[x].split()[5]) ] # find index of original model
                Rslts[index1[0]] = Rslts[index[0]] # update repeat with original results

      data = np.genfromtxt(Rslts,dtype="str",usecols=(4,5,5+ScrPltNum))

    else:
      data = np.genfromtxt(ResFile,dtype="str",skip_header=2,usecols=(4,5,5+ScrPltNum))

# Set "Failed" models to an unobtainable -1 misfit
    data[:,2][data[:,2] == "Failed"]     = "-1"
    data[:,2][data[:,2] == "to"]         = "-1"
    data[:,2][data[:,2] == "converge"]   = "-1"
    data[:,2][data[:,2] == "in"]         = "-1"
    data[:,2][data[:,2] == "required"]   = "-1"
    data[:,2][data[:,2] == "MKL"]        = "-1"
    data[:,2][data[:,2] == "iterations"] = "-1"
    data[:,2][data[:,2] == "input"]      = "-1"
    data[:,2][data[:,2] == "variable"]   = "-1"
    data[:,2][data[:,2] == "check"]      = "-1"
    data[:,2][data[:,2] == "-"]          = "-1"
    data[:,2][data[:,2] == "values"]     = "-1"

    NumMods = x1Mods*x2Mods

# Convert to float
    X1in    = data[:,0].astype(float)
    X2in    = data[:,1].astype(float)
    Mis     = data[:,2].astype(float)
    Mis[Mis == -1] = max(Mis) + max(Mis)/2       # Set "Failed" to be 50% higher than other maximum for plot

# ********** Plot Levels **********
    for iPlot in range(1,NumLvl+1):
      x1  = X1in[sum(AdaptiveLayers[0:iPlot-1])*NumMods : sum(AdaptiveLayers[0:iPlot])*NumMods]
      x2  = X2in[sum(AdaptiveLayers[0:iPlot-1])*NumMods : sum(AdaptiveLayers[0:iPlot])*NumMods]
      mis =  Mis[sum(AdaptiveLayers[0:iPlot-1])*NumMods : sum(AdaptiveLayers[0:iPlot])*NumMods]

      fig = plt.figure()
      ax = fig.add_subplot(111)
      cmhot = plt.cm.get_cmap("gnuplot")
      l = ax.scatter(x1, x2, s=20, c=mis, cmap=cmhot, marker="o")
      fig.colorbar(l,label=ScrPlt+" Score")
      plt.title("Level "+str(iPlot)+" - scored using "+ScrPlt)
      ax.set_xlabel(x1Label)
      ax.set_ylabel(x2Label)
      # ax.set_xlim([min, max])
      # ax.set_ylim([min, max])
      fig.savefig(OutPath+ScrPlt+"_Level" + str(iPlot) + "_2D."+OutType)

# ********** All models **********
    x1  = X1in
    x2  = X2in
    mis = Mis

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cmhot = plt.cm.get_cmap("gnuplot")
    l = ax.scatter(x1, x2, s=20, c=mis, cmap=cmhot, marker="o")
    fig.colorbar(l,label=ScrPlt+" Score")
    plt.title("All models - scored using "+ScrPlt)
    ax.set_xlabel(x1Label)
    ax.set_ylabel(x2Label)
    # ax.set_xlim([min, max])
    # ax.set_ylim([min, max])
    fig.savefig(OutPath+ScrPlt+"_AllModels_2D."+OutType)
    plt.show()


# ****************************** 3D Search ******************************
  elif NumVar == 3:
    print("Plotting 3D grid results")
    print("X=Variable1 Y=Variable2 Z=Variable3 colour=Misfit")

    AdaptiveLayers = [KCell]*(NumLvl-1)
    AdaptiveLayers = [1] + AdaptiveLayers

    x1Label = content[2].split(" ")[0]          # X1 variable name
    x2Label = content[7].split(" ")[0]          # X2 variable name
    x3Label = content[12].split(" ")[0]         # X3 variable name
    x1Mods  = int(content[5].split(" ")[0])     # Number of models in X1
    x2Mods  = int(content[10].split(" ")[0])    # Number of models in X2
    x3Mods  = int(content[15].split(" ")[0])    # Number of models in X2

# ********** Repeated Models? **********
    if x1Mods % 2 != 0 and x2Mods % 2 != 0 and x3Mods % 2 != 0 and ScrPltNum != 1:
      print("Selected plot option was not the misfit ranker - copying required data\n")

      with open(ResFile, "r") as f:
        Rslts = f.readlines()[2:]

      for i in range(1,((NumLvl-1) * KCell)+1): # Update every repeated model
        for line in Rslts:
          if line.split()[1] == "0":
            index1 = [x for x in range(len(Rslts)) if line.split()[1] in Rslts[x].split()[1]] # find index of repeated model
            for line2 in Rslts[:index1[0]]:
              if(line2.split()[4] == line.split()[4] and line2.split()[5] == line.split()[5] and line2.split()[6] == line.split()[6]):
                index = [x for x in range(len(Rslts)) if (line2.split()[4] in Rslts[x].split()[4] and line2.split()[5] in Rslts[x].split()[5] and line2.split()[6] in Rslts[x].split()[6]) ] # find index of original model
                Rslts[index1[0]] = Rslts[index[0]] # update repeat with original results

      data = np.genfromtxt(Rslts,dtype="str",usecols=(4,5,6,6+ScrPltNum))

    else:
      data = np.genfromtxt(ResFile,dtype="str",skip_header=2,usecols=(4,5,6,6+ScrPltNum))

# Set "Failed" models to an unobtainable -1 misfit
    data[:,3][data[:,3] == "Failed"]     = "-1"
    data[:,3][data[:,3] == "to"]         = "-1"
    data[:,3][data[:,3] == "converge"]   = "-1"
    data[:,3][data[:,3] == "in"]         = "-1"
    data[:,3][data[:,3] == "required"]   = "-1"
    data[:,3][data[:,3] == "MKL"]        = "-1"
    data[:,3][data[:,3] == "iterations"] = "-1"
    data[:,3][data[:,3] == "input"]      = "-1"
    data[:,3][data[:,3] == "variable"]   = "-1"
    data[:,3][data[:,3] == "check"]      = "-1"
    data[:,3][data[:,3] == "-"]          = "-1"
    data[:,3][data[:,3] == "values"]     = "-1"

    NumMods = x1Mods*x2Mods*x3Mods

# Convert to float
    X1in    = data[:,0].astype(float)
    X2in    = data[:,1].astype(float)
    X3in    = data[:,2].astype(float)
    Mis     = data[:,3].astype(float)
    Mis[Mis == -1] = max(Mis) + max(Mis)/2       # Set "Failed" to be 50% higher than other maximum for plot

# ********** Plot Levels **********
    for iPlot in range(1,NumLvl+1):
      x1  = X1in[sum(AdaptiveLayers[0:iPlot-1])*NumMods : sum(AdaptiveLayers[0:iPlot])*NumMods]
      x2  = X2in[sum(AdaptiveLayers[0:iPlot-1])*NumMods : sum(AdaptiveLayers[0:iPlot])*NumMods]
      x3  = X3in[sum(AdaptiveLayers[0:iPlot-1])*NumMods : sum(AdaptiveLayers[0:iPlot])*NumMods]
      mis =  Mis[sum(AdaptiveLayers[0:iPlot-1])*NumMods : sum(AdaptiveLayers[0:iPlot])*NumMods]

      fig = plt.figure()
      ax = fig.add_subplot(111, projection="3d")
      cmhot = plt.cm.get_cmap("gnuplot")
      l = ax.scatter(x1, x2, x3, s=20, c=mis, cmap=cmhot, marker="o")
      fig.colorbar(l,label=ScrPlt+" Score")
      plt.title("Level "+str(iPlot)+" - scored using "+ScrPlt)
      ax.set_xlabel(x1Label)
      ax.set_ylabel(x2Label)
      ax.set_zlabel(x3Label)
      # ax.set_xlim([min, max])
      # ax.set_ylim([min, max])
      # ax.set_zlim([min, max])
      fig.savefig(OutPath+ScrPlt+"_Level" + str(iPlot) + "_3D."+OutType)

# ********** All models **********
    x1  = X1in
    x2  = X2in
    x3  = X3in
    mis = Mis

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    cmhot = plt.cm.get_cmap("gnuplot")
    l = ax.scatter(x1, x2, x3, s=20, c=mis, cmap=cmhot, marker="o")
    fig.colorbar(l,label=ScrPlt+" Score")
    plt.title("All models - scored using "+ScrPlt)
    ax.set_xlabel(x1Label)
    ax.set_ylabel(x2Label)
    ax.set_zlabel(x3Label)
    # ax.set_xlim([min, max])
    # ax.set_ylim([min, max])
    # ax.set_zlim([min, max])
    fig.savefig(OutPath+ScrPlt+"_AllModels_3D."+OutType)
    plt.show()
