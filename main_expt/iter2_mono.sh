iterp=$1
iter=$2
rm -rf Data/mono/split20/
cp temp_$iterp.txt Data/mono/text
cp temp_$iterp.txt Data/mono/
bash steps/train_mono.sh --nj 20 Data/mono Data/lang Exp/mono_$iter
utils/mkgraph.sh Data/lang Exp/mono_$iter Exp/mono_$iter/graph
cp Exp/mono_$iter/graph/HCLG.fst Exp/mono_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight Exp/mono_$iter/graph/HCLG_wi.fst Exp/mono_$iter/graph/HCLG.fst
cp Exp/mono_$iter/graph/HCLG.fst Exp/mono_$iter/graph/HCLG_wo.fst
decode.sh --nj 20 Exp/mono_$iter/graph Data/mono Exp/mono_$iter/decode
seq 1 20 |parallel -k "grep '^f\|^m' Exp/mono_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt
