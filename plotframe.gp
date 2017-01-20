set terminal qt size 900, 900
set view 85, 0

set xrange [-2:10]
set yrange [-2:10]
set zrange [-2:10]

splot "vector.dat" using 1:2:3:4:5:6 with vectors nohead

pause -1
