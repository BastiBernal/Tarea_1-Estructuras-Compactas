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
#include "wavelet/wavelet_binary.hpp"
#include "sdsl/suffix_arrays.hpp"

#include "fib-lib/fib_tabulated.hpp"
#include "fib-lib/fib_memoized.hpp"
#include "fib-lib/fib_recursive.hpp" 
using namespace std;

const string PATH = "benchmarks/textos/";

int main(){
  AbstractFM* fm_index = nullptr;

  vector<string> archivos = {"dna.50MB","dna.100MB","dna.200MB","dblp.xml.50MB",
    "dblp.xml.100MB","dblp.xml.200MB", "sources.50MB","sources.100MB","sources.200MB"};

  vector<string> estructuras = {"fm-index-sdsl","wt_blc","fuerza_bruta","wt"};

  int prev = NULL;
  for(size_t i = 0; i < archivos.size(); ++i) {
    prev = NULL;
    string pattern = build_pattern(PATH + archivos[i],8); // Sacar un patrón pequeño del inicio del texto
    for (size_t j = 0; j < estructuras.size(); ++j){
      std:: cout << "Benchmarking " << estructuras[j] << " con el archivo " << archivos[i] << std::endl;
      // Se define el tipo de solución 
      if (estructuras[j] == "fm-index-sdsl"){
        fm_index = new FMIndexSDSL<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127> >, (1<<30), (1<<30)>>();
      } else if (estructuras[j] == "wt"){
        fm_index = new FMIndex<WaveletTreePointerless>;
      } else if (estructuras[j]== "fuerza_bruta"){
        fm_index = new FMIndex<AntiWavelet>;
      } else if (estructuras[j] == "wt_blc"){
        fm_index = new FMWaveletSDSL<sdsl::wt_blcd<>>;
      }

      {
        BenchLib::Benchmark bench;
        bench.add("Construcción de " + estructuras[j], [&fm_index, &archivos, i]() {
        fm_index->construct(PATH + archivos[i]);
        }).set_input_size(file_size(PATH + archivos[i])).set_label(archivos[i]);
   
        bench.run(20,5);
      
        if (j == 0 && i == 0) bench.write_csv("creacion_benchmark.csv");
        else bench.append_csv("creacion_benchmark.csv");
      } // Benchmark de construcción
  
      int val = fm_index->count(pattern); // Hacer un acceso para asegurar que la estructura se ha construido completamente y no hay costos de construcción ocultos en la medición de tiempo de búsqueda
      cout << "Valor de count para el patrón: " << val << endl;
      assert(val > 0); // Asegurarse de que el patrón se encuentra en el texto
      if (prev != NULL) assert(val == prev);
      prev = val; // Asegurarse de que todas las estructuras devuelven el mismo resultado

      {
        BenchLib::Benchmark bench2;
        bench2.add("Busqueda en " + estructuras[j], [&fm_index, &pattern]() {
          return fm_index->count(pattern);
        }).set_input_size(pattern.size()).set_label(archivos[i]).set_size_in_megabytes(fm_index->size_in_bytes() / (1024.0 * 1024.0));
        bench2.run(30,10);

        if (j == 0 && i == 0) bench2.write_csv("busqueda_texto_var_benchmark.csv");
        else bench2.append_csv("busqueda_texto_var_benchmark.csv");
      }
      //bench2 : calcular tiempo de busqueda del patron pequeño

      delete fm_index;
      fm_index = nullptr;
    }
  }  
  return 0;
}