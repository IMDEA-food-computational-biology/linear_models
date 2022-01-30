#!/bin/bash
while IFS= read line; do
	job_id=$(tail -n 1 ${line} | grep -o "job_id\":\"[^\"]*" | cut -d '"' -f 3)
done
