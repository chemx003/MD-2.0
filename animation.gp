stats "vector.dat" nooutput

set xrange [-2: 50]
set yrange [-2: 20]
set zrange [-2: 20]

set size ratio -1

do for [i = 1: int(STATS_blocks) - 1]{
	splot "vector.dat" index (i-1) using 1:2:3:4:5:6 with vector nohead 
	print i
	#pause 0.0001
}

pause -1
