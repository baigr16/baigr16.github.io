set encoding iso_8859_1
#set terminal  postscript enhanced color font "TimesNewRoman, 11" size 5, 4
set terminal  pngcairo truecolor enhanced lw 5.0 font "TimesNewRoman, 44" size 1920, 1680
set palette rgbformulae 22, 13, -31
# # set palette rgbformulae 7,5,15
set output 'band.png'
set border
# unset colorbox
#set title "C\\_pz" offset 0, -0.8 font "TimesNewRoman, 54"
set style data linespoints
unset ztics
unset key
# # set key outside top vertical center
# # set pointsize 0.3
set view 0,0
set xtics font "TimesNewRoman, 44"
set xtics offset 0, 0.3
set ytics font "TimesNewRoman, 44"
set ytics -10, 5, 10
set ylabel font "TimesNewRoman, 48"
set ylabel offset 1.0, 0
#set xrange [0:5.4704]
set ylabel "Energy (eV)"
set yrange [-8:8]
#set xtics ("G" 0.00000, "K" 1.695, "M" 3.938, "G" 5.407)
#plot 'BAND.dat' u 1:2 w lines lw 1.5 lc 'blue'
plot 'band_fermi.dat' u 1:2 w lines lw 1.5 lc 'blue'
