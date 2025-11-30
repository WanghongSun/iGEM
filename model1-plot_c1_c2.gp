set title "Population Dynamics of C1 and C2"
set xlabel "Time"
set ylabel "Concentration"
set xrange [0:1000]
set yrange [0:*]
set grid
set style line 1 lc rgb "blue" lw 2
set style line 2 lc rgb "red" lw 2
set terminal pngcairo size 800,600
set output "c1_c2_dynamics.png"
plot "qs_model_results.txt" using 1:2 with lines linestyle 2 title "C1", \
     "qs_model_results.txt" using 1:3 with lines linestyle 1 title "C2"

