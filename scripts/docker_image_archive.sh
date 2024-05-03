#!/bin/bash
#for image in `docker images --format "{{.Repository}}:{{.Tag}}"`; do
#    echo "bash /home/xinhe/Projects/sidbdri-utils/scripts/docker_image_archive.sh -s -i ${image}"
#done | xargs -t -n 1 -P 3 -I % bash -c "%"


function usage {
    echo""
    echo""
    echo "Usage: $0 (-s | -l) -f <file> [-i <input_image>]"
    echo "Options:"
    echo "  -s: Save the Docker image to the output file."
    echo "      -i <input_image>: Specify the Docker image to save or load."
    echo "  -l: Load the Docker image from the input file."
    echo "  -f <file>: Specify the file to save to or load from."

    exit 1
}

# Default values
input_image=""
file=""
action=""

# Read parameters
while getopts 'i:f:sl' opt; do
    case ${opt} in
        i) input_image="${OPTARG}" ;;
        f) file="${OPTARG}" ;;
        s) action="save" ;;
        l) action="load" ;;
        *) usage ;;
    esac
done

# Check if action is provided and not both -s and -l
if [[ -z "${action}" || (${action} != "save" && ${action} != "load") ]]; then
    echo "Error: Action must be provided and must be either 'save' or 'load'."
    usage
fi

# Check parameters based on action
if [[ ${action} == "save" ]]; then
    if [[ -z "${input_image}" ]]; then
        echo "Error: Input image must be provided."
        usage
    fi
    if [[ -z "${file}" ]]; then
        echo "File not provided, using input image name as file name."
        # replace '/' with '_' in the image name
        file=$(echo "${input_image}" | tr '/' '_' ).tar.gz
        echo "Output file: ${file}"
    fi

elif [[ ${action} == "load" ]]; then
    if [[ -z "${file}" ]]; then
        echo "Error: Input file must be provided."
        usage
    fi
fi

# Perform the specified action
case "${action}" in
    "save")
        docker save "${input_image}" | gzip > "${file}"
        echo "Docker image '${input_image}' saved to '${file}'."
        ;;
    "load")
        docker image load -i "${file}"
        echo "Docker image loaded from '${file}'."
        ;;
esac


function list_all_images(){
    docker images --format "{{.Repository}}:{{.Tag}}"
}