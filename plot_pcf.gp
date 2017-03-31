#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "g(r)"
set yrange [0:2.5]

#m=1; a=1;
#f(x) = a*exp(-b*(x-m))
#fit [5:18] f(x) "pcf.dat" using 1:2 via m, a, b

plot "pcf.dat" with lines #, f(x)

pause -1
