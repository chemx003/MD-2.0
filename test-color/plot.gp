set palette model RGB defined (0 'red', 1 'green')
splot "plot.dat" using 1:2:3:($4 == 0 ? 0:1) with points palette
pause -1
