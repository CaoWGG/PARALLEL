SIZE=1150
for NAME in rand10k #rgb rand10k rand100k pattern snowsingle biglittle
do
./render -b 0:1 -r cuda -f cuda_$NAME $NAME -s $SIZE
./render -b 0:1 -f ref_$NAME $NAME -s $SIZE
done


