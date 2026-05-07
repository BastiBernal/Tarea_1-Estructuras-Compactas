#include "utils-p/utils.hpp"
#include <iostream>
#include <fstream>

std::string build_pattern(std::string filename,int large){
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        std::cerr << "Error al abrir el archivo\n";
        return NULL;
    }

    std::string pattern(large, '\0');

    in.read(&pattern[0], large);
    pattern.resize(in.gcount()); // ajusta si el archivo tiene menos de N

    std::cout << pattern << std::endl;
    return pattern;
}