#!/bin/bash

set -o nounset
set -o errexit

# Returns a random 32-character ID
function get_random_id {
    cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1   
}

# Returns files matching the specified pattern(s) separated by a particular
# delimiter character.
function list_files {
    local DELIMITER=$1
    shift
    local FILES=$@
        
    LIST=$(ls -1 $FILES | tr '\n' "$DELIMITER")
    echo ${LIST%$DELIMITER}    
}
