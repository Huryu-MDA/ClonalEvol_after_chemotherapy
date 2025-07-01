#!/bin/bash

FilePathList=$1 #e.g. VCFList
Tag=$2 # e.g. Tag

FilePath=$(cat ${FilePathList} | grep ${Tag})
File_count=$(cat ${FilePathList} | grep ${Tag} | wc -l)

# echo ${File_count}

if [[ ${File_count} -eq 1 ]]
then
  echo ${FilePath}
else
  echo "FileNum_picked_up_by_tag_is_${File_count}."
  echo ${FilePath} | xargs -d " " -i echo {}
  exit 1
fi


