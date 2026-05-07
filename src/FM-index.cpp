#include "fm-index/FM-index.hpp"
#include "sdsl/suffix_arrays.hpp"
#include <algorithm>


//using namespace sdsl;
using namespace std;

sdsl::int_vector<> build_bwt(string filename){
    sdsl::int_vector<> seq;
    int32_t n;
    {
    sdsl::load_vector_from_file(seq, filename, 1);
    n = seq.size();
    seq.resize(n+1);
    n = seq.size();
    seq[n-1] = 0; // Representa el final de texto. Suele representarse por el
    // símbolo $
    }
    //cout << "Construyendo el Suffix array ..." << endl;
    sdsl::int_vector<> sa(1, 0, sdsl::bits::hi(n)+1);
    sa.resize(n);
    sdsl::algorithm::calculate_sa((const unsigned char*)seq.data(), n, sa);
    //cout << "Construyendo la BWT ..." << endl;
    sdsl::int_vector<> bwt(1, 0, 8);
    bwt.resize(n);
    int32_t to_add[2] = {(int32_t)-1,n-1};
    for (int32_t i=0; i < n; ++i)
        bwt[i] = seq[ sa[i]+to_add[sa[i]==0] ];
    //Constructor de las clases
    return bwt;
}

template <typename SDSL_structure>
int count_pattern(const SDSL_structure &wavelet, map<int, int> &C, string &pattern){
    int m = pattern.size();
    int i = m - 1;
    int c = pattern[m-1];

    int sp = C[c]+1;
    auto it = C.upper_bound(c);
    int ep;
    if (it == C.end()) ep = sp; //Si c era la ultima letra del alfabeto, entonces ep = sp
    else ep = it->second; //Si c no era la ultima letra del alfabeto, entonces ep = C[c+1]

    while ((sp<=ep) && (i>=1)){
        c = pattern[i-1];
        sp = C[c] + wavelet.rank(sp-1, c)+1;
        ep = C[c] + wavelet.rank(ep, c);
        i = i - 1;
    }
    if (sp > ep) return 0;
    return ep - sp + 1;
}

template <typename WaveletTree>
void FMIndex<WaveletTree>::construct(string filename){
    sdsl::int_vector<> bwt = build_bwt(filename);

    map<int, int> freq;
    for (int i = 0; i < bwt.size(); i++) freq[bwt[i]]++;

    int sum = 0;
    for (auto &[c, f] : freq) {
        this->C[c] = sum;
        sum += f;
    }
    this->wt.build(bwt);
}

template <typename SDSL_structure>
void FMIndexSDSL<SDSL_structure>::construct(string filename){
    sdsl::construct(this->fm, filename,1);
}

template <typename SDSL_wavelet>
void FMWaveletSDSL<SDSL_wavelet>::construct(string filename){
    sdsl::int_vector<> bwt = build_bwt(filename);

    map<int, int> freq;
    for (int i = 0; i < bwt.size(); i++) freq[bwt[i]]++;

    int sum = 0;
    for (auto &[c, f] : freq) {
        this->C[c] = sum;
        sum += f;
    }

    sdsl::construct_im(this->wavelet, bwt,1);
}

template <typename WaveletTree>
int FMIndex<WaveletTree>::count(string &pattern){
    return count_pattern(this->wt, this->C, pattern);
}

template <typename SDSL_structure>
int FMIndexSDSL<SDSL_structure>::count(string &pattern){
    return sdsl::count(this->fm, pattern.begin(), pattern.end());
}

template <typename SDSL_wavelet>
int FMWaveletSDSL<SDSL_wavelet>::count(string &pattern){
    return count_pattern(this->wavelet, this->wavelet.C, pattern);
}

template <typename WaveletTree>
int FMIndex<WaveletTree>::size_in_bytes(){
    return this->wt.size_in_bytes() + this->C.size() * (sizeof(int) * 2);
}

template <typename SDSL_structure>
int FMIndexSDSL<SDSL_structure>::size_in_bytes(){
    return sdsl::size_in_mega_bytes(this->fm) * 1024 * 1024;
}

template <typename SDSL_wavelet>
int FMWaveletSDSL<SDSL_wavelet>::size_in_bytes(){
    return sdsl::size_in_mega_bytes(this->wavelet) * 1024 * 1024 + this->C.size() * (sizeof(int) * 2);
}