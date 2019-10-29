nphone=35
cp ../$nphone/path.sh .
cp ../$nphone/decode.sh .
cp ../$nphone/iter* .
cp ../$nphone/grammar.txt .
cp ../$nphone/prepare_dict.sh .
ln -s /common/software/kaldi_new/egs/timit/s5/local local
ln -s /common/software/kaldi_new/egs/timit/s5/utils utils
ln -s /common/software/kaldi_new/egs/timit/s5/steps steps
cp ../i_201/cluster_$nphone .

less cluster_$nphone |sed 's/^/ /' |sed 's/$/ /' > cluster_35_new
less ../i_201/data/lang/phones.txt |parallel -k -j1 --colsep ' ' "sed -i 's/ {2} / {1} /' cluster_35_new"
sed -i 's/^ //' cluster_35_new
sed -i 's/ $//' cluster_35_new
awk '{for(i=1;i<=NF;i++){printf "%s C%02d",$i,NR;printf "\n"}}' cluster_35_new > map.txt
cp ../i_201/data/Mono/temp_15.txt text_old
cp text_old text_new 
less map.txt |parallel -k -j1 --colsep ' ' "sed -i 's/ {1} / {2} /g' text_new"
sed -i 's/ $//' text_new


mkdir -p data/train
cp text_new data/train/text
bash prepare_dict.sh
bash utils/prepare_lang.sh --sil-prob 0.0 --position-dependent-phones false --num-sil-states 3 data/dict "sil" data/lang_tmp data/lang
fstcompile --isymbols=data/lang/words.txt --osymbols=data/lang/words.txt --keep_isymbols=false --keep_osymbols=false grammar.txt data/lang/G.fst
mkdir data/mono
cp ../i_201/data/Mono/feats.scp data/mono/feats.scp
cp ../i_201/data/Mono/cmvn.scp data/mono/cmvn.scp
cp ../i_201/data/Mono/utt2spk data/mono/utt2spk
cp ../i_201/data/Mono/spk2utt data/mono/spk2utt
cp text_new data/mono/text
cp text_new data/mono/temp_00.txt
iter=0
. ./path.sh
bash steps/train_mono.sh --nj 20 data/mono data/lang exp/mono_$iter
utils/mkgraph.sh data/lang exp/mono_$iter exp/mono_$iter/graph
cp exp/mono_$iter/graph/HCLG.fst exp/mono_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight exp/mono_$iter/graph/HCLG_wi.fst exp/mono_$iter/graph/HCLG.fst
cp exp/mono_$iter/graph/HCLG.fst exp/mono_$iter/graph/HCLG_wo.fst
decode.sh --nj 20 exp/mono_$iter/graph data/mono exp/mono_$iter/decode
seq 1 20 |parallel -k "grep '^file' exp/mono_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt
seq 1 15 |awk '{print $1-1" "$1}' | parallel -k -j1 --colsep ' ' "bash iter1.sh {1} {2}"
cp temp_15.txt data/mono/

mkdir -p Data/mono
cp ../../../featExtract/mfccDir/train.scp Data/mono/feats.scp
cp ../../../featExtract/mfccDir/cmvn_train.scp Data/mono/cmvn.scp
cp ../../../featExtract/spk2utt_train Data/mono/spk2utt
cp ../../../featExtract/utt2spk_train Data/mono/utt2spk

cp -r data/lang Data/
rm Data/lang/G.fst
gunzip -c data/lm/lm_phone_bg.arpa.gz | egrep -v '<s> <s>|</s> <s>|</s> </s>' | arpa2fst - | fstprint |utils/eps2disambig.pl | utils/s2eps.pl | fstcompile --isymbols=data/lang/words.txt --osymbols=data/lang/words.txt  --keep_isymbols=false --keep_osymbols=false | fstrmepsilon | fstarcsort --sort_type=ilabel > Data/lang/G.fst

ln -s /p2/sre/karthik/jasa_rebuttal/main_expt/gd/1state/$nphone/exp/mono_15/40.mdl Exp/mono_0/final.mdl
ln -s /p2/sre/karthik/jasa_rebuttal/main_expt/gd/1state/$nphone/exp/mono_15/tree Exp/mono_0/tree
utils/mkgraph.sh Data/lang Exp/mono_0 Exp/mono_0/graph

rm temp_*
iter=0
cp Exp/mono_$iter/graph/HCLG.fst Exp/mono_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight Exp/mono_$iter/graph/HCLG_wi.fst Exp/mono_$iter/graph/HCLG.fst
cp Exp/mono_$iter/graph/HCLG.fst Exp/mono_$iter/graph/HCLG_wo.fst
decode.sh --nj 20 Exp/mono_$iter/graph Data/mono Exp/mono_$iter/decode
seq 1 20 |parallel -k "grep '^f\|^m' Exp/mono_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt
seq 1 15 |awk '{print $1-1" "$1}' | parallel -k -j1 --colsep ' ' "bash iter2_mono.sh {1} {2}"

cp temp_15.txt Data/mono/text
steps/align_si.sh --boost-silence 1.25 --nj 20 Data/mono Data/lang Exp/mono_15 Exp/mono_15_ali

cp -r Data/mono/ Data/tri1
rm -rf Data/tri1/split20/
rm -rf Data/tri1/temp_*  
rm -rf Data/tri1/text  
cp temp_15.txt Data/tri1/temp_00.txt
cp temp_15.txt Data/tri1/text
steps/train_deltas.sh 2500 15000 Data/tri1 Data/lang Exp/mono_15_ali Exp/tri1_$iter
utils/mkgraph.sh Data/lang Exp/tri1_$iter Exp/tri1_$iter/graph
cp Exp/tri1_$iter/graph/HCLG.fst Exp/tri1_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight Exp/tri1_$iter/graph/HCLG_wi.fst Exp/tri1_$iter/graph/HCLG.fst
cp Exp/tri1_$iter/graph/HCLG.fst Exp/tri1_$iter/graph/HCLG_wo.fst
decode.sh --nj 20 Exp/tri1_$iter/graph Data/tri1 Exp/tri1_$iter/decode
rm temp_*.txt
seq 1 20 |parallel -k "grep '^f\|^m' Exp/tri1_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt
seq 1 15 |awk '{print $1-1" "$1}' | parallel -k -j1 --colsep ' ' "bash iter2_tri1.sh {1} {2}"
steps/align_si.sh --nj 20 Data/tri1 Data/lang Exp/tri1_15 Exp/tri1_15_ali

cp -r Data/tri1 Data/tri2
rm -rf Data/tri2/split20/
rm -rf Data/tri2/temp_*  
rm -rf Data/tri2/text  

cp temp_15.txt Data/tri2/temp_00.txt
cp temp_15.txt Data/tri2/text
steps/train_lda_mllt.sh --splice-opts "--left-context=3 --right-context=3" 2500 15000 Data/tri2 Data/lang Exp/tri1_15_ali Exp/tri2_$iter
utils/mkgraph.sh Data/lang Exp/tri2_$iter Exp/tri2_$iter/graph
cp Exp/tri2_$iter/graph/HCLG.fst Exp/tri2_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight Exp/tri2_$iter/graph/HCLG_wi.fst Exp/tri2_$iter/graph/HCLG.fst
cp Exp/tri2_$iter/graph/HCLG.fst Exp/tri2_$iter/graph/HCLG_wo.fst
decode.sh --nj 20 Exp/tri2_$iter/graph Data/tri2 Exp/tri2_$iter/decode
rm temp_*.txt
seq 1 20 |parallel -k "grep '^f\|^m' Exp/tri2_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt
seq 1 15 |awk '{print $1-1" "$1}' | parallel -k -j1 --colsep ' ' "bash iter2_tri2.sh {1} {2}"
steps/align_si.sh --nj 20 --use-graphs true Data/tri2 Data/lang Exp/tri2_15 Exp/tri2_15_ali

cp -r Data/tri2 Data/tri3
rm -rf Data/tri3/split20/
rm -rf Data/tri3/temp_*
rm -rf Data/tri3/text

cp temp_15.txt Data/tri3/temp_00.txt
cp temp_15.txt Data/tri3/text
steps/train_sat.sh 2500 15000 Data/tri3 Data/lang Exp/tri2_15 Exp/tri3_$iter
utils/mkgraph.sh Data/lang Exp/tri3_$iter Exp/tri3_$iter/graph
cp Exp/tri3_$iter/graph/HCLG.fst Exp/tri3_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight Exp/tri3_$iter/graph/HCLG_wi.fst Exp/tri3_$iter/graph/HCLG.fst
cp Exp/tri3_$iter/graph/HCLG.fst Exp/tri3_$iter/graph/HCLG_wo.fst
decode_fmllr.sh --nj 20 Exp/tri3_$iter/graph Data/tri3 Exp/tri3_$iter/decode
rm temp_*.txt
seq 1 20 |parallel -k "grep '^f\|^m' Exp/tri3_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt
seq 1 15 |awk '{print $1-1" "$1}' | parallel -k -j1 --colsep ' ' "bash iter2_tri3.sh {1} {2}"


