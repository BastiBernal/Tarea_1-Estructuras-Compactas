set datafile separator ","
set encoding utf8

set terminal pdfcairo enhanced size 12cm,8cm
set output "draws/text_var1.pdf"

set title "Espacio vs Tiempo de respuesta"
set xtics rotate by -45
set xlabel "Espacio de la estructura (MB)"
set ylabel "Tiempo promedio (μs)"

#set logscale y
#set yrange [7e6:9e7]
set grid
set key outside


# Lineas mas visibles
set style line 1 lw 2 pt 7 ps 1.2
set style line 2 lw 2 pt 5 ps 1.2
set style line 3 lw 2 pt 9 ps 1.2

set format x "%.0s%c"
set format y "%.0f"
set xtics offset 0,-1

plot \
    "< tail -n +2 ../busqueda_texto_var_benchmark.csv.csv | grep 'wt,' | sort -t, -k3,3n" \
    using 8:4 \
    with points pt 7 ps 1 title "Pointerless", \
\
    "< tail -n +2 ../busqueda_texto_var_benchmark.csv.csv | grep 'wt_blc,' | sort -t, -k3,3n" \
    using 8:4 \
    with points pt 5 ps 1 title "Balanced", \
\
    "< tail -n +2 ../busqueda_texto_var_benchmark.csv.csv | grep 'sdsl,' | sort -t, -k3,3n" \
    using 8:4 \
    with points pt 9 ps 1 title "Huffman", \
\
    "< tail -n +2 ../busqueda_texto_var_benchmark.csv.csv | grep 'fuerza_' | sort -t, -k3,3n" \
    using 8:4 \
    with points pt 9 ps 1 title "Fuerza bruta"