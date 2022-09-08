#!/bin/bash

#Just create the data wher to attach the info.
mkdir $1 && cd $_
touch SRR_Acc_List.txt
touch miniMetadata.txt
touch table_$1.csv
cp /data/users/mbotos/PhDdataSets/Datasets/getSRRdata.sh .
