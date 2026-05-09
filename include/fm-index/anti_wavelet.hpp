#pragma once

#include "sdsl/suffix_arrays.hpp"

// Implementación de un FM-index que no utiliza un wavelet tree, sino que almacena el BWT directamente
// y hace un conteo de ocurrencias lineal para cada símbolo.
class AntiWavelet{
private:
    sdsl::int_vector<> text;
public:
    void build(sdsl::int_vector<> &bwt);
    int rank(int i, int c);
    int size_in_bytes();
};