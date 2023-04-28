wgsim -e 0 -N 8192 -r 0 -R 0 -X 0 -S 1 -d 0 -A 0 viral_20230330/index/viral.fna viral_20230330/shit1.fq viral_20230330/shit2.fq
./fqt transpose -i viral_20230330/shit1.fq -o viral_20230330/shit1 -l 70


make && ./main -l 70 -i viral_20230330/shit1 -x viral_20230330/index/viral.fna -t ./ -o viral_20230330/output.txt -n viral_20230330/taxonomy/nodes.dmp