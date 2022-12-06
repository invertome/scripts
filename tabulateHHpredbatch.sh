for file in ./*.hhr ; do
        python2 ~/bin/tabulateHHpred.py -i $file -o ${file}.tsv
done

cat *uniref30.hhr.tsv > uniref30_results.tsv 

cat *pdb70.hhr.tsv > pdb70_results.tsv 

cat *scop95.hhr.tsv > scop95_results.tsv

cat *pfam.hhr.tsv > pfam_results.tsv 

cat *_results.tsv > all_results.tsv
