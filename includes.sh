#!/usr/bin/env bash

export BASE_DIR
BASE_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )

set -a
source "${BASE_DIR}/src/ensembl_utils.sh"
source "${BASE_DIR}/src/ncbi_utils.sh"
source "${BASE_DIR}/src/utils.sh"
set +a