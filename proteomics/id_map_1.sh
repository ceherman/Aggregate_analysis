#!/bin/bash

readarray -t ids_arr < ./generated_tables/master_accn.txt
printf -v ids_list '%s,' "${ids_arr[@]}"

curl --request POST 'https://rest.uniprot.org/idmapping/run' --form 'ids='"$ids_list"'' --form 'from="RefSeq_Protein"' --form 'to="UniProtKB"' | json_pp -json_opt pretty,canonical > job_name.json

jobId=($(jq -r '.jobId' job_name.json))
url="https://rest.uniprot.org/idmapping/status/$jobId"
curl -i $url
