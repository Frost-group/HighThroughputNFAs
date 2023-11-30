# Run Espaloma on a list of SDFs
# Messing around with filenames here rather than within the Python

for i
do 
 obabel $i -O ${i%.*}.pdb
 cp -a $i input.sdf
 cp -a ${i%.*}.pdb input.pdb
 python3 run_espaloma.py
 mv output.pdb ${i%.*}_espaloma_md.pdb
done

