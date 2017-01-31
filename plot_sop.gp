#set terminal png
#set output "energy.png" 

set xlabel "Timestep"
set ylabel "Scalar order parameter"
set yrange [-1:1]

plot "sop.dat" using 1:2 title 'x' with lines, \
	 "sop.dat" using 1:3 title 'y' with lines, \
	 "sop.dat" using 1:4 title 'z' with lines

pause -1
