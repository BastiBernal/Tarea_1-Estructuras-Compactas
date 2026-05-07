#pragma once

#include "sdsl/suffix_arrays.hpp"

class AntiWavelet{
private:
    sdsl::int_vector<> text;
public:
    void build(sdsl::int_vector<> &bwt);
    int rank(int i, int c);
    int size_bytes();
};