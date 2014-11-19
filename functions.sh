#!/bin/bash

set -o nounset
set -o errexit

function get_random_id {
    cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1   
}

function list_files {
    local DELIMITER=$1
    shift
    local FILES=$@
        
    LIST=$(ls -1 $FILES | tr '\n' "$DELIMITER")
    echo ${LIST%$DELIMITER}    
}
