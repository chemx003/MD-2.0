set xrange [0:50]
set yrange [0:20]
set zrange [6:11]

set size ratio -1

splot "director.dat" using ($1-$4/2):($2-$5/2):($3-$6/2) \
	:4:5:6 with vector nohead

pause -1 
