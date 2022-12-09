#UNIREF30

for file in ./*.fa ; do
       hhblits -B 1 -b 1 -p 80 -Z 1 -E 1E-10 -hide_cons -cpu 10 -maxmem 32 -i $file -d /work/y0rgan/databases/hhsuite/UniRef30_2020_06 -o ${file}.uniref30.hhr
done


#scop95

for file in ./*.fa ; do
       hhblits -B 1 -b 1 -p 80 -Z 1 -E 1E-10 -hide_cons -cpu 10 -maxmem 32 -i $file -d /work/y0rgan/databases/hhsuite/scop95 -o ${file}.scop95.hhr
done



#pdb70

for file in ./*.fa ; do
       hhblits -B 1 -b 1 -p 80 -Z 1 -E 1E-10 -hide_cons -cpu 10 -maxmem 32 -i $file -d /work/y0rgan/databases/hhsuite/pdb70 -o ${file}.pdb70.hhr
done


#pfam

for file in ./*.fa ; do
       hhblits -B 1 -b 1 -p 80 -Z 1 -E 1E-10 -hide_cons -cpu 10 -maxmem 32 -i $file -d /work/y0rgan/databases/hhsuite/pfam -o ${file}.pfam.hhr
done
