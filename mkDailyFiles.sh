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
FILES='../../../pluma'

if [ ! -d "latlonpairs/" ]
then
	mkdir latlonpairs
else
        rm -rf latlonpairs
        mkdir latlonpairs
fi

for FILE in `ls $FILES`
do
    #We want the first and last hour record to get starting and ending hour.
    #First record (no header)
    awk -F'\t' '{if (NR == 2) {print $2}}' $FILES/$FILE > firstTime.txt 
    #Last record
    awk -F'\t' END'{print $2}' $FILES/$FILE > lastTime.txt
    #Now manipulate the strings to get only the hours
    FIRST=`head firstTime.txt`
    LAST=`head lastTime.txt`
    HOUR_INIT=${FIRST:0:2}
    HOUR_FIN=${LAST:0:2}
    #Add one hour, since we need the whole thing
    HOUR_FIN=$(($HOUR_FIN+1))
    
    #Since day should not change, get any of the lines that's not a header. Say 3.
    awk -F'\t' '{if (NR == 3) {print $1}}' $FILES/$FILE > date.txt 
    #Get the day and month from the date string
    DATE=`head date.txt`
    DAY=${DATE:8}
    MONTH=${DATE:5:2}
    
    #Now comes the interesting part, using sed
    sed 's:'MONTH':'${MONTH}':g' getWindHeightLevels.py.template > mkDailyFile.py
    sed 's:'DAY':'${DAY}':g' mkDailyFile.py > mkDailyFile2.py
    sed 's:'PLUMTEXT':'${FILE}':g' mkDailyFile2.py > mkDailyFile.py
    sed 's:'HOUR_INIT':'${HOUR_INIT}':g' mkDailyFile.py > mkDailyFile2.py
    sed 's:'HOUR_FIN':'${HOUR_FIN}':g' mkDailyFile2.py > mkDailyFile.py
    
    mkdir latlonpairs/${FILE}
    python mkDailyFile.py
    cat latlonpairs/???.csv > latlonpairs/${FILE}/transecto_d${DAY}_m${MONTH}.txt
    
done











