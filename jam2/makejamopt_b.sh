seed=$(( $1 + 0 ))

echo "Main:numberOfEvents = 2000  ! number of events to generate" > jamopt.inp
echo "Random:seed $seed" >> jamopt.inp

cat input/besex_7p7.inp >> jamopt.inp 
