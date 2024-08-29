seed=$(( $1 + 0 ))

echo "Main:numberOfEvents = 2000  ! number of events to generate" > jamopt.inp
echo "Random:seed $seed" >> jamopt.inp

cat input/besex.inp >> jamopt.inp 
#cat input/4p5_Cascade_rmSpectator.inp >> jamopt.inp 
