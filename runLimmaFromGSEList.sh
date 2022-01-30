#!/bin/bash

gseList=$1
if [ -z "$1" ]
then

        echo "Introduce gse list"
        exit 1

fi

python3 create_DM.py $gseList
while read line
do
       echo $line	
	./DownloadDataMatrix.sh "GSE"$line "microarray_CEL/"
	Rscript runLimma.R "./DM/""$line""_DM.tsv" "microarray_CEL/GSE""$line""_series_matrix.txt.gz" "./DM/""$line""_node_ids.tsv" 

done < $gseList
