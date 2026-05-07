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

int count_pattern(const sdsl::int_vector<> &bwt, const vector<size_t> &C, const string &pattern){
    int m = pattern.size();
    int c = pattern[m-1];
    int sp = C[c];
    int ep = C[c+1]-1;
    for (int i=m-2; i>=0 && sp <= ep; --i){
        c = pattern[i];
        sp = C[c] + sdsl::rank_support_v<1>(&bwt)(sp-1, c);
        ep = C[c] + sdsl::rank_support_v<1>(&bwt)(ep, c) - 1;
    }
    if (sp > ep) return 0;
    return ep - sp + 1;
}

void FMIndex::construct(string filename){
    sdsl::int_vector<> bwt = build_bwt(filename);

    map<int, int> freq;
    for (int i = 0; i < n; i++) freq[bwt[i]]++;

    int sum = 0;
    for (auto &[c, f] : freq) {
        C[c] = sum;
        sum += f;
    }
    //Constructor del wavelet
}

void FMIndexSDSL::construct(string filename){
    sdsl::construct(this->fm, filename,1);
}

void FMWaveletSDSL::construct(string filename){
    sdsl::int_vector<> bwt = build_bwt(filename);
    sdsl::construct_im(this->wavelet, bwt,1);
}