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

	Rscript runLimma.R "./DM/""$line""_DM.tsv" "./HTS_gene_counts/GSE""$line""_counts.tsv" "./DM/""$line""_node_ids.tsv" 

done < $gseList
