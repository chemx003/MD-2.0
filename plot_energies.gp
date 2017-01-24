#set terminal png
#set output "energy.png" 

set xlabel "Timestep"
set ylabel "Energy"

set yrange [-10:10]

plot "energy.dat" using 1:2 title 'V' with lines, \
	 "energy.dat" using 1:3 title 'K' with lines, \
	 "energy.dat" using 1:4 title 'E' with lines

pause -1
