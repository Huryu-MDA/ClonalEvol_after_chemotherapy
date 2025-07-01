#!/bin/bash

JobIDList=$1

cat ${JobIDList} | grep "^Job" |cut -f 2 -d "<" | cut -f 1 -d ">" | xargs -i printf "done({})&&" | rev | cut -f 3- -d "&" | rev
