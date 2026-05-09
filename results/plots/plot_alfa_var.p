set datafile separator ","
set encoding utf8

set terminal pdfcairo enhanced size 14cm,8cm

set xlabel "Largo del patron"
set ylabel "Tiempo promedio (μs)"

set grid ytics
set key outside bottom center horizontal

set xtics rotate by 45 right
#set format x "%.0e"

# ---------------- DNA ----------------

set output "draws/dna_benchmark.pdf"

set title "Benchmark sobre textos DNA"

plot \
    "< tail -n +2 ../busqueda_patron_var_benchmark.csv.csv | grep 'dna' | grep 'wt,'" \
    using 3:4 with linespoints lw 2 pt 7 ps 1 title "Pointerless", \
\
    "< tail -n +2 ../busqueda_patron_var_benchmark.csv.csv | grep 'dna' | grep 'wt_blc,'" \
    using 3:4 with linespoints lw 2 pt 5 ps 1 title "Balanced", \
\
    "< tail -n +2 ../busqueda_patron_var_benchmark.csv.csv | grep 'dna' | grep 'sdsl,'" \
    using 3:4 with linespoints lw 2 pt 9 ps 1 title "Huffman"


# ---------------- XML ----------------

set output "draws/xml_benchmark.pdf"

set title "Benchmark sobre textos XML"

plot \
    "< tail -n +2 ../busqueda_patron_var_benchmark.csv.csv | grep 'xml' | grep 'wt,'" \
    using 3:4 with linespoints lw 2 pt 7 ps 1 title "Pointerless", \
\
    "< tail -n +2 ../busqueda_patron_var_benchmark.csv.csv | grep 'xml' | grep 'wt_blc,'" \
    using 3:4 with linespoints lw 2 pt 5 ps 1 title "Balanced", \
\
    "< tail -n +2 ../busqueda_patron_var_benchmark.csv.csv | grep 'xml' | grep 'sdsl,'" \
    using 3:4 with linespoints lw 2 pt 9 ps 1 title "Huffman"

# ---------------- SOURCES ----------------

set output "draws/sources_benchmark.pdf"

set title "Benchmark sobre textos Sources"

plot \
    "< tail -n +2 ../busqueda_patron_var_benchmark.csv.csv | grep 'sources' | grep 'wt,'" \
    using 3:4 with linespoints lw 2 pt 7 ps 1 title "Pointerless", \
\
    "< tail -n +2 ../busqueda_patron_var_benchmark.csv.csv | grep 'sources' | grep 'wt_blc,'" \
    using 3:4 with linespoints lw 2 pt 5 ps 1 title "Balanced", \
\
    "< tail -n +2 ../busqueda_patron_var_benchmark.csv.csv | grep 'sources' | grep 'sdsl,'" \
    using 3:4 with linespoints lw 2 pt 9 ps 1 title "Huffman"