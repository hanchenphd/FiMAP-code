echo "Input gzipped file $1..."
echo "Output file prefix $2..."
echo "Chunk number: $3..."

nlines=`zcat $1 | wc -l`
nperbatch=$(( $nlines/$3 ))
for (( i=1; i<$3; i++ ))
do
    zcat $1 | awk -v begin=$(( ($i-1)*nperbatch+1 )) -v end=$(( $i*nperbatch )) 'NR>=begin && NR<=end {print}' > $2.part$i
done
zcat $1 | awk -v begin=$(( ($3-1)*nperbatch+1 )) -v end=$nlines 'NR>=begin && NR<=end {print}' > $2.part$3

