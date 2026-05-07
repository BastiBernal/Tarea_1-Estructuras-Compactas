#pragma once

#include <string>


class AbstractFM{
public:
    virtual void construct(std::string filename) = 0;
    virtual int count(std::string pattern) = 0;
    virtual ~AbstractFM() = default ;
};

template <typename WaveletTree>
class FMIndex: public AbstractFM{

private:
    WaveletTree wt;
    std::map<int,int> C; 

public:
    void construct(std::string filename) override;
    int count(std::string pattern) override;
    //void load_wt(string file, string text);

};

template <typename SDSL_structure>
class FMIndexSDSL: public AbstractFM{
    
private:
    SDSL_structure fm ;

public:
    void construct(std::string filename) override;
    int count(std::string pattern) override;
};

template <typename SDSL_wavelet>
class FMWaveletSDSL: public AbstractFM{

private:
    SDSL_wavelet wavelet;

public:
    void construct(std::string filename) override;
    int count(std::string pattern) override;

};