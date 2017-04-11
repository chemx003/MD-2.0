set xrange [-1:50]
set yrange [-1:20]

set size ratio -1

plot "director.dat" using ($1-$3/2):($2-$4/2):3:4 with vector nohead

pause -1 
