#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "g(r)"
set yrange [0:1.5]

m=1; a=1;
gauss(x) = a*exp(-b*(x-m))
fit [1:5] gauss(x) "pcf.dat" using 1:2 via m, a, b

plot "pcf.dat" with lines, gauss(x)

pause -1
