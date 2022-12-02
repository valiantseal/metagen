#!/usr/bin/env bash

sed 's/ 1:N:0/\/3 1:N:0/g' merged_prep_temp.fastq  > merged_prep_temp-tag.fastq