#!/bin/sh

#  aledated.sh
#  
#
#  Created by Laura Mulvey on 12.07.21.


#sample=$1

counter=1
#for i in $sample
#do
#docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve $PWD/output/$i
docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve $PWD/data/parasites.trees
#let counter=counter+1
#done

#  
#docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEml_undated $PWD/host.trees $PWD/parasites.trees.ale sample=1





#docker run -v $PWD:$PWD -w $PWD boussau/alesuite ALEobserve $PWD/$1/$i