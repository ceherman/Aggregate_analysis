#!/bin/bash

# Returns the curl request that you submitted - just a check

jobId=($(jq -r '.jobId' job_name.json))
url="https://rest.uniprot.org/idmapping/details/$jobId"
curl -i $url
