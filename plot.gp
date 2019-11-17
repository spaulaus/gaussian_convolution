reset
load '~/.gnuplot'

file='temp.dat'

s1=1.0
s2=2.0

amp=1.0
mu=0.0
sigma=sqrt(s1*s1+s2*s2)

plot 'temp.dat' u 1:2 ps 2  lw 3 t 'g1', '' u 1:3 ps 2 lw 3 t 'g2', '' u 1:4 ps 2 lw 3 t 'g1(*)g2', gaussian(x) lw 3 t 'theory'

set terminal postscript eps enhanced color solid "Helvetica, 20"
set output 'conv.eps'
replot