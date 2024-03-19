#!/bin/bash

prefix=$1

tar czvf ${prefix}.tar.gz ${prefix}/
rm -r ${prefix}/
mkdir -p targz_storage && mv ${prefix}.tar.gz targz_storage/.
