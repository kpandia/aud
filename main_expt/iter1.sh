iterp=$1
iter=$2
rm -rf data/mono/split20/
cp temp_$iterp.txt data/mono/text
cp temp_$iterp.txt data/mono/
bash steps/train_mono.sh --nj 20 data/mono data/lang exp/mono_$iter
utils/mkgraph.sh data/lang exp/mono_$iter exp/mono_$iter/graph
cp exp/mono_$iter/graph/HCLG.fst exp/mono_$iter/graph/HCLG_wi.fst
fstmap --map_type=rmweight exp/mono_$iter/graph/HCLG_wi.fst exp/mono_$iter/graph/HCLG.fst
cp exp/mono_$iter/graph/HCLG.fst exp/mono_$iter/graph/HCLG_wo.fst
decode.sh --nj 20 exp/mono_$iter/graph data/mono exp/mono_$iter/decode
seq 1 20 |parallel -k "grep '^file' exp/mono_$iter/decode/log/decode.{1}.log |grep -v WARN |grep -v Log" > temp_$iter.txt
