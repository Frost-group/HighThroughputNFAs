IFS=$'\n'

for l in ` cat ../master.tsv | grep -v "^#" `
do
    name=` echo $l | cut -f1 -d$'\t' `
    echo $name
    smile=` echo $l | cut -f2 -d$'\t' `
    
    echo "$smile" > $name.smi

    ../bin/smi2sdf -p ../bin/mmxconst.prm $name.smi -o $name.sdf 
done

