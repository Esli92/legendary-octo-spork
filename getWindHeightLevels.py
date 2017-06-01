#!/bin/bash

#split_day_in_hours.sh

#This program takes daily files from ADCIRC output and then splits the file into several files, one for each timestep. 
#This program is meant to be used to preprocess data to use with oilspill.m 
#Programmer Oscar Jurado (ojurado@ciencias.unam.mx)
#Creation date: 23-May-2017

#------------------Requisites------------------------------------------------------------
#NCO installed to use ncks function. 

#-----------------Version---------------------------------------------------------------
#v1.0 23/May/17 Program is created

#----------------Known issues-----------------------------------------------------------
#Requires that user knows number of values per timestep and manualy changes for file. 

#-----------------Local directories----------------------------------------------------- 
SCRIPT_DIR=`pwd`
DATA_DIR=./caso_01

#----------------Dependencies used-----------------------------------------------------
#ncks

#-----------------BEGIN PROGRAM--------------------------------------------------------
