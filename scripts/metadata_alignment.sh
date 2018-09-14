
for file in *final.out; do
sample=${file%%.*}
cat $file | grep -vE '%$' | grep -vE ':$' > $sample.txt
done
