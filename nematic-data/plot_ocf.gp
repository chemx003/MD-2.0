#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "o(r)"
set logscale y

m=1; a=1;
fit(x) = a*exp(-b*(x-m))
fit [1:18] fit(x) "ocf.dat" using 1:2 via m, a, b

plot "ocf.dat" with lines, fit(x)

pause -1
