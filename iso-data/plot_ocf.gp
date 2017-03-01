#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "o(r)"
set logscale y

m=1; a=1; b=1;
f(x) = a*exp(-b*(x-m))
fit [0:5] f(x) "ocf.dat" using 1:2 via m, a,  b

plot "ocf.dat" with lines, f(x)

pause -1
