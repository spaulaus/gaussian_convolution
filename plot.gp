reset
load '~/.gnuplot'

file=temp.dat

s1=1.0
s2=2.0

amp=1.0
mu=0.0
sigma=s1*s1+s2*s2


plot 'temp.dat' u 1:2, '' u 1:3, '' u 1:4, gaussian(x)