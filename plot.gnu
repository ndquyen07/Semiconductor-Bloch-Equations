# Gnuplot script for Semiconductor Bloch Equations results
# Run with: gnuplot plot.gnu

set terminal pngcairo size 1600,1200 font 'Arial,14'
set output 'gnuplot_results.png'

set multiplot layout 2,2 title "Semiconductor Bloch Equations - Results"

# Plot 1: Population dynamics
set xlabel "Time"
set ylabel "Population"
set grid
plot 'data/population.dat' using 1:2 with lines lw 2 lc rgb "blue" title 'Total Population'

# Plot 2: Polarization (Real and Imaginary)
set xlabel "Time"
set ylabel "Polarization"
set grid
plot 'data/polarization.dat' using 1:2 with lines lw 2 lc rgb "red" title 'Re[P]', \
     'data/polarization.dat' using 1:3 with lines lw 2 lc rgb "blue" title 'Im[P]'

# Plot 3: Absorption spectrum comparison
set xlabel "Energy"
set ylabel "Absorption"
set grid
set key top right
plot 'data/absorption_comparison_with.dat' using 1:2 with lines lw 2 lc rgb "blue" title 'With Coulomb', \
     'data/absorption_comparison_without.dat' using 1:2 with lines lw 2 lc rgb "red" title 'Without Coulomb'

# Plot 4: Occupation f_n for selected levels (n=1, 2, 3, 4, 5)
set xlabel "Time"
set ylabel "Occupation f_n"
set grid
set key top left
plot 'data/fn.dat' using 1:2 with lines lw 2 title 'n=1', \
     'data/fn.dat' using 1:3 with lines lw 2 title 'n=2', \
     'data/fn.dat' using 1:4 with lines lw 2 title 'n=3', \
     'data/fn.dat' using 1:5 with lines lw 2 title 'n=4', \
     'data/fn.dat' using 1:6 with lines lw 2 title 'n=5'

unset multiplot

print "Plot saved to: gnuplot_results.png"
