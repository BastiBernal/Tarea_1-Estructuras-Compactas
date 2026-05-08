#include <cassert>
#include <variant>
#include <iostream>
#include <fstream>
#include <string>

#include "bench-lib/benchmark.hpp"
#include "fm-index/FM-index.hpp"
#include "utils-p/utils.hpp"
#include "fm-index/anti_wavelet.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "wavelet/wavelet_tree.hpp"

#include "fib-lib/fib_tabulated.hpp"
#include "fib-lib/fib_memoized.hpp"
#include "fib-lib/fib_recursive.hpp" 
using namespace std;

const string PATH = "benchmarks/textos/";

int main() {
  AbstractFM* fm_index = nullptr;
  vector<int> sizes= {32,64,128,256,512};
  vector<string> archivos = {"dna.50MB", "dblp.xml.50MB", "sources.50MB"};
  vector<string> estructuras = {"fm-index-sdsl","wt","fuerza_bruta","wt_blc"};
  int prev = NULL;
  for(size_t i = 0; i < archivos.size(); ++i) {
    for (size_t k = 0; k < estructuras.size(); ++k){
      if (estructuras[k] == "fm-index-sdsl"){
            fm_index = new FMIndexSDSL<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127> >, (1<<30), (1<<30)>>();
          } else if (estructuras[k] == "wt"){
            fm_index = new FMIndex<WaveletTreeBinary>;
          } else if (estructuras[k]== "fuerza_bruta"){
            fm_index = new FMIndex<AntiWavelet>;
          } else if (estructuras[k] == "wt_blc"){
            fm_index = new FMWaveletSDSL<sdsl::wt_blcd<>>;
          }
          fm_index->construct(PATH + archivos[i]);
      for (size_t j = 0; j < sizes.size(); j++){
          //prev = NULL;
          string pattern = build_pattern(PATH + archivos[i], sizes[j]);
          int val = fm_index->count(pattern); // Hacer un acceso para asegurar que la estructura se ha construido completamente y no hay costos de construcción ocultos en la medición de tiempo de búsqueda
          assert(val > 0); // Asegurarse de que el patrón se encuentra en el texto
          //if (prev != NULL) assert(val == prev);
          //prev = val; // Asegurarse de que todas las estructuras devuelven el mismo resultado
          {
            BenchLib::Benchmark bench;
            bench.add("Busqueda en" + estructuras[k], [&fm_index, &pattern]() {
              return fm_index->count(pattern);
            }).set_input_size(pattern.size()).set_label(archivos[i]).set_size_in_megabytes(fm_index->size_in_bytes() / (1024.0 * 1024.0));
            bench.run(30,10);

            if (j == 0 && i == 0 && k == 0) bench.write_csv("busqueda_patron_var_benchmark2.csv");
            else bench.append_csv("busqueda_patron_var_benchmark2.csv");
          }
        }
        delete fm_index;
        fm_index = nullptr;

    }
    
  }
  return 0;
}