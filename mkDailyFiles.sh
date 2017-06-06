#!/bin/bash

#===========================================================================================
#Program Name: mkDailyFiles.sh

#Description: Wrapper to use with getWindHeightLevels.py and process all input files.
#Language: bash

#Programmer Oscar Jurado (ojurado@ciencias.unam.mx)
#Creation date: June-2017

#------------------Use------------------------------------------------------------------
#chmod +x mkDailyFiles.sh
#./mkDailyFiles.sh

#------------------Requisites------------------------------------------------------------
#WRF output data in ../../../salidas/
#pluma data in ../../../pluma
#getWindHeightLevels.py in same folder with requisites fulfilled.
#-----------------Version---------------------------------------------------------------
#v1.0 June/17 Program is created

#----------------Known issues-----------------------------------------------------------


#----------------Dependencies and Libraries----------------------------------------------


#---------------Import-------------------------------------------------------------------


#-----------------Local directories----------------------------------------------------- 
WRF_dir='../../../salidas/'
FILES='/home/esli/Documentos/WRF/experimentos/Toluca/pluma'

for FILE in `ls $FILES`
do

    awk -F'\t' '{if (NR != 1) {print $2}}' $FILE > hours.txt 
    awk -F'\t' '{if (NR == 3) {print $2}}' $FILE > date.txt 
done











