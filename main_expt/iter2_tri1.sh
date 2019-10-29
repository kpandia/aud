iterp=$1
iter=$2
rm -rf Data/tri1/split20/
cp temp_$iterp.txt Data/tri1/text 
cp temp_$iterp.txt Data/tri1/
steps/train_deltas.sh 2500 15000 Data/tri1 Data/lang Exp/tri1_$iterp Exp/tri1_$iter
utils/mkgraph.sh Data/lang Exp/tri1_$iter Exp/tri1_$iter/graph
cp Exp/tri1_$iter/graph/HCLG.fst Exp/tri1_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight Exp/tri1_$iter/graph/HCLG_wi.fst Exp/tri1_$iter/graph/HCLG.fst
cp Exp/tri1_$iter/graph/HCLG.fst Exp/tri1_$iter/graph/HCLG_wo.fst
decode.sh --nj 20 Exp/tri1_$iter/graph Data/tri1 Exp/tri1_$iter/decode
seq 1 20 |parallel -k "grep '^f\|^m' Exp/tri1_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt

#cp temp_15.txt Data/mono/text
#steps/align_si.sh --boost-silence 1.25 --nj 20 Data/mono Data/lang Exp/mono_15 Exp/mono_15_ali

