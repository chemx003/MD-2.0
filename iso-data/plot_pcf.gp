#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "g(r)"
set logscale y

m=1; a=1;
fit(x) = a*exp(-(x-m))
fit [1:5] fit(x) "pcf.dat" using 1:2 via m, a

plot "pcf.dat" with lines, fit(x)

pause -1
