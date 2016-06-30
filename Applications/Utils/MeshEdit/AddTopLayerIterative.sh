#!/bin/bash

set -x
src_file=$1
dest=/tmp/01.vtu
mv $src_file $dest

OGS_PATH=bin

for i in `seq 1 5`
do
    $OGS_PATH/AddTopLayer -i /tmp/0${i}.vtu -o /tmp/0${i}_A.vtu -t 10.712595
    $OGS_PATH/reviseMesh -i /tmp/0${i}_A.vtu -o /tmp/0${i}_AC.vtu -s -e 0.001
    $OGS_PATH/editMaterialID -i /tmp/0${i}_AC.vtu -o /tmp/0$(( ${i}+1 )).vtu -r -m 717 -n 716
    $OGS_PATH/CheckTopLayerConnectivity -i /tmp/0$(( ${i}+1 )).vtu
done
