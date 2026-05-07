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
  vector<string> archivos = {"dna.50MB", "xml.50MB", "sources.50MB"};
  vector<string> estructuras = {"fb","fm","wt","wthf"};

  for(size_t i = 0; i < archivos.size(); ++i) {
    for (size_t j = 0; j < sizes.size(); j++){
       string pattern = build_pattern(PATH + archivos[i], sizes[j]);
       for (size_t k = 0; k < estructuras.size(); ++k){
          if (estructuras[j] == "fm-index-sdsl"){
            cout << "Construyendo FM-index de SDSL" << endl;
            fm_index = new FMIndexSDSL<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127> >, (1<<30), (1<<30)>>();
          } else if (estructuras[j] == "wt"){
            fm_index = new FMIndex<WaveletTreeBinary>;
          } else if (estructuras[j]== "fuerza_bruta"){
            fm_index = new FMIndex<AntiWavelet>;
          } else if (estructuras[j] == "wt_blc"){
            fm_index = new FMWaveletSDSL<sdsl::wt_blcd<>>;
          }
          fm_index->construct(PATH + archivos[i]);
          
          {
            BenchLib::Benchmark bench;
            bench.add("Busqueda en" + estructuras[k], [&fm_index, &pattern]() {
              return fm_index->count(pattern);
            }).set_input_size(pattern.size()).set_label(archivos[i]);
            bench.run(50,20);

            if (j == 0 && i == 0) bench.write_csv("busqueda_patron_var_benchmark.csv");
            else bench.append_csv("busqueda_patron_var_benchmark.csv");
          }

          delete fm_index;
          fm_index = nullptr;
        }

    }
    
  }
  return 0;
}