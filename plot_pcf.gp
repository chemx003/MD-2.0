#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "g(r)"

plot "pcf.dat" with lines

pause -1
