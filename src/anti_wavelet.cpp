#include "fm-index/anti_wavelet.hpp"
#include "sdsl/suffix_arrays.hpp"

void AntiWavelet::build(sdsl::int_vector<> &bwt){
    this->text = bwt;
}

int AntiWavelet::rank(int i, int c){
    int count = 0;
    for (int j = 0; j < i; j++){
        if (text[j] == c) count++;
    }
    return count;
}

int AntiWavelet::size_in_bytes(){
    return size_in_mega_bytes(this->text) * 1024 * 1024;
}