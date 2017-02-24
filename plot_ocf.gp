#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "o(r)"

m=1; a=1;
gauss(x) = a*exp(-b*(x-m))
fit [1:5] gauss(x) "ocf.dat" using 1:2 via m, a, b 

plot "ocf.dat" with lines, gauss(x)

pause -1
