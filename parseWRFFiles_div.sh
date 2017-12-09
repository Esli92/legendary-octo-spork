#mkStationReadsSeasons.sh
#Programa para autogenerar scripts que generan diagramas de Taylor.
#Programador Oscar Jurado
#Fecha de creacion: Nov/2016

#------------------Requisitos------------------------------------------------------------
#Directorio de estaciones a usar
#Archivos con datos estacionales en pares obs/pron para cada estacion

#-----------------Versiones---------------------------------------------------------------
#v1.0 Se crea el programa. 


#----------------Problemas Conocidos-----------------------------------------------------


#----------------Directorios Locales, cambiar si es necesario----------------------------
DIRECTORIO_SALIDAS=../dataFiles


rm -rf netcdf_div
mkdir -p netcdf_div


# figlet TAYLOR DIAGRAM
# echo "GENERATING R SCRIPTS FOR THE TAYLOR DIAGRAMS, PLEASE WAIT!"
rm ../figures/*.png

for FILE in `ls $DIRECTORIO_SALIDAS`
do
        sed 's/'FILE'/'${FILE}'/g' getDivergenceB.py.template > getDivergenceB.py
        python getDivergenceB.py
        mkdir -p ../figures/${FILE}
        mv ../figures/*.png ../figures/${FILE}
done


figlet DONE

