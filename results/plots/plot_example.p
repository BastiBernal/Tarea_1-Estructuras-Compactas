set datafile separator ","
set encoding utf8

set terminal pdfcairo enhanced size 12cm,8cm
set output "draws/creacion.pdf"

set title "Tiempo de creacion vs tamaño"
set xtics rotate by -45
set xlabel "Largo del archivo (Bytes)"
set ylabel "Tiempo promedio (μs)"

set grid
set key outside

set format x "%.0s%cB"
set format y "%.0f"
set xtics offset 0,-1

# Lineas mas visibles
set style line 1 lw 2 pt 7 ps 1.2
set style line 2 lw 2 pt 5 ps 1.2
set style line 3 lw 2 pt 9 ps 1.2

plot \
    "< tail -n +2 ../creacion_benchmark.csv.csv | grep 'wt,' | sort -t, -k3,3n" using 3:4:7 with yerrorlines pt 7 ps 0.4 title "Pointerless", \
    "< tail -n +2 ../creacion_benchmark.csv.csv | grep 'wt_blc,' | sort -t, -k3,3n" using 3:4:7 with yerrorlines pt 7 ps 0.4 title "Balanced", \
    "< tail -n +2 ../creacion_benchmark.csv.csv | grep 'fuerza_bruta,' | sort -t, -k3,3n" using 3:4:7 with yerrorlines pt 7 ps 0.4 title "Fuerza Bruta",\
    "< tail -n +2 ../creacion_benchmark.csv.csv | grep 'sdsl,' | sort -t, -k3,3n" using 3:4:7 with yerrorlines pt 7 ps 0.4 title "Huffman"