#!/bin/bash

download_by_id () {

	local id=$( basename $1 | cut -f 1 -d "_" )
	#TODO: check if there is job id
	local job_id=$(grep -o "\"job_id\":\"[^\"]*\"" $1 | cut -f 4 -d '"')
	curl -s -X GET --header "Accept: application/json" --header "user_key:e896c326dff0f857e0ecc8a91323433b" "https://api.clue.io/api/jobs/findByJobId/${job_id}" > .tmpCMAP_responses/${id}_response.tmp
	#TODO: check if there is download url
	download_url=$(grep -o "\"download_url\":\"[^\"]*\"" .tmpCMAP_responses/${id}_response.tmp | cut -f 4 -d '"')
	curl -X GET "http:"${download_url} > ${outPutDir}"/"${id}"_CMAP_data.tar.gz"
}

mkdir .tmpCMAP_responses

if [ -p /dev/stdin ]; then

	echo "Data was piped to this script!"
	
	mkdir "$(pwd)"/CMAP_data/ 
	outPutDir="$(pwd)"/CMAP_data/
        
	while IFS= read line; do
		download_by_id $line
        done
else

	echo "No input was found on stdin, skipping!"
	response_dir=$1
	
	if [ -z "$1" ]
	then

		echo "Introduce CMAP response directory"
		rm -r .tmpCMAP_responses
		exit 1

	fi

	if [ -z $2 ]
	then

		mkdir "$(pwd)"/CMAP_data/ 
		outPutDir="$(pwd)"/CMAP_data/
	else

		outPutDir=$2

	fi

	for i in $( ls $response_dir/* )
	do	
		download_by_id $i
	done


fi

rm -r .tmpCMAP_responses

