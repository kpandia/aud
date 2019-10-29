awk '{if(NF==2){print prevname" "NR-1" "NR-prevno-1;prevno=NR;prevname=$1}}' ../featDir/mfccDir/feats.txt >temp.txt
#remove ]  manually
less temp.txt |parallel -k --colsep ' ' "head -{2} ../featDir/mfccDir/feats.txt |tail -{3} > featDir/{1}.mfc"
ls featDir/*.mfc |parallel -k "awk 'NR==FNR{sum++;next;}{if(FNR==1){print \"39 \"sum};print \$0}' {1} {1} > {1.}.mfcc"
less r.sh |parallel -k
less feat.lst |parallel -k "less outDir/{1/.}.out" >mat.txt

seq 1 13748 |parallel -k "awk -v n={1} 'NR==n{print \$0;exit}' mat.txt |tr ' ' '\n' |awk '{print \$1\" \"NR}' |sort -nk1 |awk '{print \$2}' |head -5|tr '\n' ' '" > tempDel1.txt
seq 1 13748 |parallel -k "python -c \"print('{1} ' * 5)\"" >tempDel2.txt

awk '{print NR" "$1}' feat.lst > map.txt
less map.txt |parallel -k -j1 --colsep ' ' "sed -i 's/ {1} / {2} /' cluster_3_4_mapped.txt"
