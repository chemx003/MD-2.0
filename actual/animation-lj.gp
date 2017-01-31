stats "vector.dat" nooutput

set xrange [-2: 10]
set yrange [-2: 10]
set zrange [-2: 10]

do for [i = 1: int(STATS_blocks) - 1]{
	splot "vector.dat" index (i-1) using 1:2:3
	pause 0.0001
}

pause -1
