iterp=$1
iter=$2

rm -rf Data/tri2/split20/
cp temp_$iterp.txt Data/tri2/text 
cp temp_$iterp.txt Data/tri2/
steps/train_lda_mllt.sh --splice-opts "--left-context=3 --right-context=3" 2500 15000 Data/tri2 Data/lang Exp/tri2_$iterp Exp/tri2_$iter
utils/mkgraph.sh Data/lang Exp/tri2_$iter Exp/tri2_$iter/graph
cp Exp/tri2_$iter/graph/HCLG.fst Exp/tri2_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight Exp/tri2_$iter/graph/HCLG_wi.fst Exp/tri2_$iter/graph/HCLG.fst
cp Exp/tri2_$iter/graph/HCLG.fst Exp/tri2_$iter/graph/HCLG_wo.fst
decode.sh --nj 20 Exp/tri2_$iter/graph Data/tri2 Exp/tri2_$iter/decode
seq 1 20 |parallel -k "grep '^f\|^m' Exp/tri2_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt

#steps/align_si.sh --nj 20 Data/tri1 Data/lang Exp/tri1_15 Exp/tri1_15_ali
