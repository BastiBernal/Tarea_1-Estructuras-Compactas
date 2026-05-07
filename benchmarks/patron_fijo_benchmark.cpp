#include <cassert>
#include <variant>
#include <iostream>
#include <fstream>
#include <string>

#include "bench-lib/benchmark.hpp"
#include "fm-index/FM-index.hpp"
#include "utils-p/utils.hpp"

#include "fib-lib/fib_tabulated.hpp"
#include "fib-lib/fib_memoized.hpp"
#include "fib-lib/fib_recursive.hpp" 
using namespace std;

const string PATH = "textos/";

int main(){
  AbstractFM* fm_index = nullptr;

  vector<string> archivos = {"dna.50MB","dna.100MB","dna.200MB","dblp.xml.50MB",
    "dblp.xml.100MB","dblp.xml.200MB", "sources.50MB","sources.100MB","sources.200MB"};

  vector<string> estructuras = {"fuerza_bruta","fm-index-sdsl","wt","wthf"};

  for(size_t i = 0; i < archivos.size(); ++i) {
    string pattern = build_pattern(PATH + archivos[i], 16); // Sacar un patrón pequeño del inicio del texto
    for (size_t j = 0; j < estructuras.size(); ++j){
      if (estructuras[j] == "fm-index-sdsl"){
        fm_index = new FMIndexSDSL<sdsl::csa_wt<wt_huff<rrr_vector<127> >, 512, 1024>>();
      } else if (estructuras[j] == "wt"){
        fm_index = new FMWaveletSDSL<>;
      } else if (estructuras[j]== "fm"){
        fm_index = new FMIndex<>;
      } else if (estructuras[j] == ""){

      } else if (estructuras[j]== ""){

      }
      {
        BenchLib::Benchmark bench;
        bench.add("Construcción de " + estructuras[j], [&fm_index, &archivos, i]() {
        fm_index->construct(PATH + archivos[i]);
        }).set_input_size(archivos[i].size()).set_label(archivos[i]);

        bench.run(30,10);
        if (j == 0 && i == 0) bench.write_csv("creacion_benchmark.csv");
        else bench.append_csv("creacion_benchmark.csv");
      }
      
      {
        BenchLib::Benchmark bench2;
        bench2.add("Busqueda en" + estructuras[j], [&fm_index, &pattern]() {
          return fm_index->count(pattern);
        }).set_input_size(pattern.size()).set_label(archivos[i]);
        bench2.run(50,20);

        if (j == 0 && i == 0) bench2.write_csv("busqueda_texto_var_benchmark.csv");
        else bench2.append_csv("busqueda_texto_var_benchmark.csv");
      }
      //bench2 : calcular tiempo de busqueda del patron pequeño

      delete fm_index;
      fm_index = nullptr;
    }
  return 0;
}