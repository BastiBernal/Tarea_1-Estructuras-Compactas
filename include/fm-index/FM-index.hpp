#pragma once

#include <string>
#include <map>
#include "utils-p/utils.hpp"
#include "sdsl/suffix_arrays.hpp"
#include <algorithm>

//Interfaz común para encapsular los distintos wavelet tree
class AbstractFM{
public:
    virtual void construct(std::string filename) = 0;
    virtual int count(std::string &pattern) = 0;
    virtual int size_in_bytes() = 0;
    virtual ~AbstractFM() = default ;
};

// FMIndex con wavelet tree propio
template <typename WaveletTree>
class FMIndex: public AbstractFM{

private:
    WaveletTree wt;
    std::map<uint8_t, int> C; 
    int size;

public:
    void construct(std::string filename) override;
    int count(std::string &pattern) override;
    int size_in_bytes() override;

};

// FMIndex que solo encapsula el FMIndex de SDSL
template <typename SDSL_structure>
class FMIndexSDSL: public AbstractFM{
    
private:
    SDSL_structure fm ;

public:
    void construct(std::string filename) override;
    int count(std::string &pattern) override;
    int size_in_bytes() override;
};

// FMIndex con wavelet tree de SDSL pero con la interfaz de nuestro FMIndex
template <typename SDSL_wavelet>
class FMWaveletSDSL: public AbstractFM{

private:
    SDSL_wavelet wavelet;
    std::map<uint8_t, int> C;
    int size ;

public:
    void construct(std::string filename) override;
    int count(std::string &pattern) override;
    int size_in_bytes() override;

};


/* Implementación de las clases , se hace dentro del .hpp para evitar problemas con el template al compilar*/

template <typename WaveletTree>
void FMIndex<WaveletTree>::construct(std::string filename){
    sdsl::int_vector<> bwt = build_bwt(filename);

    std::map<uint8_t, int> freq;
    for (size_t i = 0; i < bwt.size(); ++i) freq[(uint8_t)bwt[i]]++;

    int sum = 0;
    for (auto &[c, f] : freq) {
        this->C[c] = sum;
        sum += f;
    }
    this->size = bwt.size();
    this->wt.build(bwt);
}

template <typename SDSL_structure>
void FMIndexSDSL<SDSL_structure>::construct(std::string filename){
    sdsl::construct(this->fm, filename,1);
}

template <typename SDSL_wavelet>
void FMWaveletSDSL<SDSL_wavelet>::construct(std::string filename){
    sdsl::int_vector<> bwt = build_bwt(filename);

    std::map<uint8_t, int> freq;
    for (size_t i = 0; i < bwt.size(); ++i) freq[(uint8_t)bwt[i]]++;

    int sum = 0;
    for (auto &[c, f] : freq) {
        this->C[c] = sum;
        sum += f;
    }
    this->size = bwt.size();
    sdsl::construct_im(this->wavelet, bwt);
}

template <typename WaveletTree>
int FMIndex<WaveletTree>::count(std::string &pattern){
    return count_pattern<WaveletTree>(this->wt, this->C, pattern, this->size);
}

template <typename SDSL_structure>
int FMIndexSDSL<SDSL_structure>::count(std::string &pattern){
    return sdsl::count(this->fm, pattern.begin(), pattern.end());
}

template <typename SDSL_wavelet>
int FMWaveletSDSL<SDSL_wavelet>::count(std::string &pattern){
    return count_pattern_sdsl<SDSL_wavelet>(this->wavelet, this->C, pattern, this->size);
}

template <typename WaveletTree>
int FMIndex<WaveletTree>::size_in_bytes(){
    return this->wt.size_in_bytes() + this->C.size() * (sizeof(int) * 2) + sizeof(int);
}

template <typename SDSL_structure>
int FMIndexSDSL<SDSL_structure>::size_in_bytes(){
    return sdsl::size_in_mega_bytes(this->fm) * 1024 * 1024;
}

template <typename SDSL_wavelet>
int FMWaveletSDSL<SDSL_wavelet>::size_in_bytes(){
    return sdsl::size_in_mega_bytes(this->wavelet) * 1024 * 1024 + this->C.size() * (sizeof(int) * 2) + sizeof(int);
}