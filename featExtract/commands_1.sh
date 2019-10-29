cut -d' ' -f1 0_01.txt |parallel -k "grep -w {1} 0_01.txt |awk '{print \$2\" \"\$3; prev=\$3; for(i=4;i<=NF;i++){print prev\" \"\$i; prev=\$i}}' > bndDir/{1}.bnd"
ls bndDir/*.bnd |parallel -k "awk -v file={1} '{printf \"%s %3.2f %3.2f\n\",file,\$1,\$2}' {1}" |awk '{printf "file_%05d %s \n",NR,$0}' > segments_all
awk '{if($4-$3 > 0.3){print $0}}' segments_all > segments
	## handling the segments at the END of each file 
awk -v p2="faem0_si1392" '{if($2!=p2){print p1" "p2" "p3" "p4-0.01}else{print p1" "p2" "p3" "p4}; p1=$1; p2=$2; p3=$3; p4=$4}' segments > segments1
extract-feature-segments --frame-length=25 --frame-shift=10 --min-segment-length=0.075 scp:mfccDir/train.scp segments1 ark,scp:mfccDir/feats.ark,mfccDir/feats.scp

apply-cmvn --utt2spk=ark:utt2spk scp:../../featExtract/mfccDir/cmvn_train.scp scp:mfccDir/feats.scp ark:- | add-deltas  ark:- ark,t:mfccDir/feats.txt
