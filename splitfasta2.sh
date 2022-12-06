while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        file=$(echo ${line#>} | cut -d' ' -f 1)
        outfile=${file}.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < input.fasta
