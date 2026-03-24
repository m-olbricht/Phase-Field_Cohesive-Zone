#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#
# create user subroutine object files and shared libraries
#
# Stephan Roth, TU Bergakademie Freiberg, 25.04.2024
#
#####################################################################################################
#
# prefix
prefix="CZUELlib_PFF_CZ_240425"
#
#####################################################################################################
#
# all subroutines
#
# array of all files, consider file order!
declare -a fileArray
# general and common used modules
fileArray[1]="ABQinterface.f90"
fileArray[2]="UEL_lib_PF_240228.f90"
fileArray[3]="BasicModulesUELlib_PF_240125.f90"
fileArray[4]="TensorModule_240405.f90"
fileArray[5]="DegradationModule_PFF_240425.f90"
fileArray[6]="InterfaceEnergyModule.f90"
# phase-field
fileArray[7]="PhaseField_Parameters_220316.f90"
fileArray[8]="FreeEnergyModule_PFF_240103.f90"
fileArray[9]="PhaseFieldC_modules_PFF_240123.f90"
fileArray[10]="PF_UELlibmodules_PFF_240103.f90"
fileArray[11]="PFUEL_lib_PFF_240228.f90"
# cohesive zone
fileArray[12]="CZ_Parameters_240103.f90"
fileArray[13]="CohesiveEnergyModule_PFF_240126.f90"
fileArray[14]="ContactEnergyModule_PFF_240126.f90"
fileArray[15]="PenaltyEnergyModule_PFF_240109.f90"
fileArray[16]="UMAT_TSL_UELlib_240126.f90"
fileArray[17]="CZUELlibmodules_240126.f90"
fileArray[18]="CZUEL_lib_240126.f90"
# main file
fileArray[19]="USER_MAIN_240103.f90"
#
# number of all files
numFiles=${#fileArray[*]} 
#
# combine all files in userfile-B.f
rm userfile-B.f
for (( filenum = 1; filenum <= $numFiles; filenum++ ))
  do
  cat ${fileArray[$filenum]} >> userfile-B.f
done
#
compilerOptions="-c -fPIC -qopenmp -extend-source -nostandard-realloc-lhs -WB -I%I"
ifort $compilerOptions userfile-B.f
mv userfile-B.o "$prefix".o
#
rm userfile-B.f *.mod
#
#####################################################################################################
