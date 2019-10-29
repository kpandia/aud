ls /SpeakerID-RIC/karthik/j1/timit/bndDir/train/*.bnd |parallel -k "awk -v file={1} '{printf \"%s %3.2f %3.2f\n\",file,(\$1-6)/100,(\$2-6)/100}' {1}" |awk '{printf "file_%05d %s \n",NR,$0}' > segments_all
awk '{if($4-$3 > 0.3){print $0}}' segments_all > segments

compute-mfcc-feats --verbose=2 --config=conf/mfcc.conf scp,p:wav.scp ark:- | copy-feats --compress=true ark:- ark,scp:mfccDir/train.ark,mfccDir/train.scp
compute-cmvn-stats --spk2utt=ark:mfccDir/spk2utt scp:mfccDir/train.scp ark,scp:mfccDir/cmvn_train.ark,mfccDir/cmvn_train.scp
extract-feature-segments --frame-length=25 --frame-shift=10 --min-segment-length=0.05 scp:mfccDir/train.scp segments_kw ark,scp:mfccDir/kw.ark,mfccDir/kw.scp

