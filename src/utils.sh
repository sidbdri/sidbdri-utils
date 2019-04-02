#!/usr/bin/env bash

# Returns files matching the specified pattern(s) separated by a particular
# delimiter character.
function list_files {
    local DELIMITER=$1
    shift
    local FILES=$@

    LIST=$(ls -1 $FILES | tr '\n' "$DELIMITER")
    echo ${LIST%$DELIMITER}
}
