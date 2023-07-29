#!/bin/bash

jobId=($(jq -r '.jobId' job_name.json))
url="https://rest.uniprot.org/idmapping/results/$jobId"
echo $url
curl -s $url | json_pp -json_opt pretty,canonical > ids_mapped.json
