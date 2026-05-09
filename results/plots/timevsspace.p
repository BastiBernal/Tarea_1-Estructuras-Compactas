set datafile separator ","
set encoding utf8

set terminal pdfcairo enhanced size 14cm,8cm
set output "draws/text_size_benchmark.pdf"

unset title

set xlabel "Tamaño del texto (MB)"
set ylabel "Tiempo promedio (μs)"

set grid ytics
set key outside bottom center horizontal

set xtics (50,100,200)

# estilos claros
set style line 1 pt 7 ps 1.5 lw 2 lc rgb "#1f77b4"
set style line 2 pt 5 ps 1.5 lw 2 lc rgb "#ff7f0e"
set style line 3 pt 9 ps 1.5 lw 2 lc rgb "#2ca02c"

plot \
"< tail -n +2 ../busqueda_texto_var_benchmark.csv.csv | grep 'wt,'" \
using (strstrt(strcol(2), "50MB") ? 50 : strstrt(strcol(2), "100MB") ? 100 : 200):4 \
with points pt 7 ps 1.5 title "Pointerless", \
\
"< tail -n +2 ../busqueda_texto_var_benchmark.csv.csv | grep 'wt_blc,'" \
using (strstrt(strcol(2), "50MB") ? 50 : strstrt(strcol(2), "100MB") ? 100 : 200):4 \
with points pt 5 ps 1.5 title "Balanced", \
\
"< tail -n +2 ../busqueda_texto_var_benchmark.csv.csv | grep 'sdsl,'" \
using (strstrt(strcol(2), "50MB") ? 50 : strstrt(strcol(2), "100MB") ? 100 : 200):4 \
with points pt 9 ps 1.5 title "Huffman"