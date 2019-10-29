set linetype 1 lw 2 lc rgb "blue" pointtype 6
set yrange [0:15]
set xrange [0:15]
plot "co.dat" using 1:2 with points
pause 0.1
reread 