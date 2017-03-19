#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "o(r)"
set xrange [0:11]
set yrange [0:1]
#set logscale y

#m=1; a=1;
#fit(x) = a*exp(-b*(x-m))
#fit [0:9] fit(x) "ocf.dat" using 1:2 via m, a, b

plot "ocf.dat" with lines#, fit(x)

pause -1
