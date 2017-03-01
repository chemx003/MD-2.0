#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "g(r)"
#set logscale y

m=1; a=1;
fit(x) = a*exp(-b*(x-m))
fit [5:10] fit(x) "pcf.dat" using 1:2 via m, a, b

plot "pcf.dat" with lines, fit(x)

pause -1
