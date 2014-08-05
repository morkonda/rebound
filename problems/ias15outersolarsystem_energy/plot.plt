#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 3in,2in
set xlabel "timestep [days]"
set ylabel "relative energy after 100 orbits"
set logscale xy
set autoscale fix
set yrange [1e-16:0.1]
set key right bottom

set st d lp

plot \
"energy_wh.txt" t "WH (REBOUND)",  \
"energy_ias15.txt" t "IAS15", \
1e-3*(x/1000.)**2 t "dt^{2}", \
1e-12*(x/1000.)**15 t "dt^{15}", \
