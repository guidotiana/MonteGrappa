#!/bin/zsh

rm temp_plot.gp
touch temp_plot.gp

for i in $(seq 0 99); do



name="moves/$i.dat"
echo "set terminal png" > temp_plot.gp
echo "set output 'plotmoves/"$i".png'" >> temp_plot.gp
echo "set xrange[0:170]" >> temp_plot.gp

echo "set title 'mpivot"$i"'" >> temp_plot.gp
echo "p '"$name"' u 1:2 w p" >>temp_plot.gp
gnuplot temp_plot.gp

	
done
