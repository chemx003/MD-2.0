stats "vector.dat" nooutput

set xrange [-2: 10]
set yrange [-2: 10]

do for [i = 1: int(STATS_blocks) - 1]{
	splot "vector.dat" index (i-1) using 1:2:3:4:5:6 with vectors nohead
	pause 0.0001
}
