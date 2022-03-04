#!/bin/sh

# $1 should be folder with simulated parasite datasets
#  
#
#  Created by Laura Mulvey on 15.02.21.
#  
mkdir pales
files=$(ls $1)
counter=1
for i in $files
do
echo $i
docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve $PWD/$1/$i
mv $1/*.tre.ale pales
let counter=counter+1
done

mkdir probs

files=$(ls pales)
counter=1
for i in $files
do
docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated $PWD/data/h_mcc.tre $PWD/pales/$i
mv $i.uml_rec probs
rm $i.uTs
rm $i.ucons_tree
#mv $counter_ales.ul..rec probs
let counter=counter+1
done


files=$(ls probs)
for i in $files
do
sed -n '114p' probs/$i >> probs.csv
done

