set terminal pngcairo enhanced font 'SimHei,12' size 600,600
set output 'X_vs_time.png'
set title "二肽浓度随时间变化"
set xlabel "时间（小时）"
set ylabel "二肽浓度（μg/mL）"
set xrange [-5:105]
set yrange [-80:2150]  
set grid
plot "model2_igem_whu_results.txt" using 1:3 with lines lw 2 lc rgb 'blue' title 'X(t)'