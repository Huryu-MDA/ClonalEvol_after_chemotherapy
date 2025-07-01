#!/bin/bash

FilePathList=$1 #e.g. VCFList
Tag=$2 # e.g. Tag

FilePath=$(cat ${FilePathList} | grep -e "${Tag}")
File_count=$(cat ${FilePathList} | grep -e "${Tag}" | wc -l)

# echo ${File_count}

if [[ ${File_count} -eq 1 ]]
then
  exit 0
else
  echo ${Tag}
  echo "FileNum_picked_up_by_tag_is_${File_count}." > /dev/stderr
  echo ${FilePath} | xargs -d " " -i echo {} > /dev/stderr
  exit 1
fi


