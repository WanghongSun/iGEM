set terminal pngcairo enhanced font 'SimHei,12' size 600,600
set output 'N_vs_time.png'
set title "大肠杆菌数量随时间变化"
set xlabel "时间（小时）"
set ylabel "大肠杆菌数量（CFU）"
set xrange [-5:105]
set yrange [0:110000]
set grid
plot "model2_igem_whu_results.txt" using 1:2 with lines lw 2 lc rgb 'purple' title 'N(t)'