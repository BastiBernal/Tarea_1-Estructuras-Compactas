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

#include "fib-lib/fib_tabulated.hpp"
#include "fib-lib/fib_memoized.hpp"
#include "fib-lib/fib_recursive.hpp" 
using namespace std;

const string PATH = "benchmarks/textos/";

int main(){
  AbstractFM* fm_index = nullptr;

  vector<string> archivos = {"dna.50MB"};//,"dna.100MB","dna.200MB","dblp.xml.50MB",
    //"dblp.xml.100MB","dblp.xml.200MB", "sources.50MB","sources.100MB","sources.200MB"};

  vector<string> estructuras = {/*fm-index-sdsl",/*"wt",*/"wt_blc"/*,"fuerza_bruta"*/};

  for(size_t i = 0; i < archivos.size(); ++i) {
    string pattern = build_pattern(PATH + archivos[i], 16); // Sacar un patrón pequeño del inicio del texto
    for (size_t j = 0; j < estructuras.size(); ++j){
      std:: cout << "Benchmarking " << estructuras[j] << " con el archivo " << archivos[i] << std::endl;
      if (estructuras[j] == "fm-index-sdsl"){
        cout << "Construyendo FM-index de SDSL" << endl;
        fm_index = new FMIndexSDSL<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<127> >, (1<<30), (1<<30)>>();
      } else if (estructuras[j] == "wt"){
        continue;
        //fm_index = new FMIndex<>;
      } else if (estructuras[j]== "fuerza_bruta"){
        fm_index = new FMIndex<AntiWavelet>;
      } else if (estructuras[j] == "wt_blc"){
        fm_index = new FMWaveletSDSL<sdsl::wt_blcd<>>;
      }
      cout << "Construyendo la estructura..." << endl;
      {
        BenchLib::Benchmark bench;
        bench.add("Construcción de " + estructuras[j], [&fm_index, &archivos, i]() {
        fm_index->construct(PATH + archivos[i]);
        }).set_input_size(file_size(PATH + archivos[i])).set_label(archivos[i]);
        cout << "Ejecutando benchmark..." << endl;
        bench.run(1,1);
        cout << "Escribiendo resultados..." << endl;
        if (j == 0 && i == 0) bench.write_csv("creacion_benchmark.csv");
        else bench.append_csv("creacion_benchmark.csv");
      }
      
      {
        BenchLib::Benchmark bench2;
        bench2.add("Busqueda en" + estructuras[j], [&fm_index, &pattern]() {
          return fm_index->count(pattern);
        }).set_input_size(pattern.size()).set_label(archivos[i]).set_size_in_megabytes(fm_index->size_in_bytes() / (1024.0 * 1024.0));
        bench2.run(50,20);

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