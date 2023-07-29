#!/bin/bash

# Returns the status

jobId=($(jq -r '.jobId' job_name.json))
url="https://rest.uniprot.org/idmapping/status/$jobId"
curl -i $url
