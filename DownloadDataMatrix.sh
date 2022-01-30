#!/bin/bash

seriesName=$1
fileName="${seriesName}_series_matrix.txt.gz"
if [ -z "$1" ]
then

	echo "Introduce series name"
	exit 1

fi

if [ -z $2 ]
then

	outPutDir="$(pwd)"
else

	outPutDir=$2

fi


#check if file exists
if [ -f "${outPutDir}/${fileName}" ]
then
	echo "${outPutDir}/${fileName} already exists"
else

	#creating the url
	#after series, the folder name is the GSE up to the last three characters
	#then nnn. 
	#${var:start:stop} is for string slicing in bash
	#${#var} gets you character length

	url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${seriesName:0:$((${#seriesName}-3))}nnn/${seriesName}/matrix/${fileName}"
	echo "Looking for $url"
	curl "$url" -o "$outPutDir/$fileName"

fi
