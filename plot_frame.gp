set terminal qt size 900, 900
set view 85, 0
set palette model RGB defined (0 'green', 1 'red')
set cbrange [0:1]

set xrange [-2:10]
set yrange [-2:10]
set zrange [-2:10]

splot "vector.dat" u 1:2:3:4:5:6 with vector nohead 

pause -1
