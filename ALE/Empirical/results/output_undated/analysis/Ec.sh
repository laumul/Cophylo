files=$(ls DTL)

echo 'Duplication Transfer Losses Speciations' > output.csv

for i in $files
do

cd DTL/$i

FOO=$(cat parasites_READY.trees.ale.uml_rec |  sed -n '1014p')

 printf '%s\n' "${FOO//Total/}" >> ../../output.csv

cd ../../

done
