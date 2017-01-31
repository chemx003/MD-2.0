set xlabel "timestep"
set ylabel "temp"

plot "temp.dat" u 1:2 with lines

pause -1
