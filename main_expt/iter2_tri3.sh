iterp=$1
iter=$2
rm -rf Data/tri3/split20/
cp temp_$iterp.txt Data/tri3/text 
cp temp_$iterp.txt Data/tri3/
steps/train_sat.sh 2500 15000 Data/tri3 Data/lang Exp/tri3_$iterp Exp/tri3_$iter
utils/mkgraph.sh Data/lang Exp/tri3_$iter Exp/tri3_$iter/graph
cp Exp/tri3_$iter/graph/HCLG.fst Exp/tri3_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight Exp/tri3_$iter/graph/HCLG_wi.fst Exp/tri3_$iter/graph/HCLG.fst
cp Exp/tri3_$iter/graph/HCLG.fst Exp/tri3_$iter/graph/HCLG_wo.fst
decode_fmllr.sh --nj 20 Exp/tri3_$iter/graph Data/tri3 Exp/tri3_$iter/decode
seq 1 20 |parallel -k "grep '^f\|^m' Exp/tri3_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt
#steps/align_si.sh --nj 20 --use-graphs true Data/tri2 Data/lang Exp/tri2_15 Exp/tri2_15_ali
