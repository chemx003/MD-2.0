#set terminal png
#set output "pcf.png"

set xlabel "r"
set ylabel "o(r)"

plot "ocf.dat" with lines

pause -1
