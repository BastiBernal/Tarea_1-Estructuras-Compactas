#pragma once 

#include "bitvector/bitvector.hpp"
#include "sdsl/suffix_arrays.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

std::string build_pattern(std::string filename,int large);


struct ResultadoEntropia {
	std::unordered_map<char, double> distribucion;
	double entropia = 0.0;
};

// Calcula la distribución de probabilidades de cada caracter en el texto
// y la entropía de Shannon (en bits).
//
// Si el texto está vacío, retorna distribucion vacía y entropia = 0.0.
ResultadoEntropia CalcularEntropia(const std::string& texto);

// ---------------- BitVectors: selector automático ----------------

enum class BitVectorTipo {
	Plain,
	SD,
	RRR,
};

struct BitVectorDecision {
	BitVectorTipo tipo = BitVectorTipo::Plain;
	double entropia = 0.0; // bits/símbolo (Shannon), sobre {0,1}
	double p1 = 0.0;       // proporción de 1s
	uint64_t n = 0;
	uint64_t ones = 0;
};

// Entropía de Shannon de un bitvector
// de {0,1}: H(p) = -p log2 p - (1-p) log2(1-p).
//
// Para p=0 o p=1 retorna 0.
double CalcularEntropiaBits(double p1);

// Decide qué estructura conviene según una heurística basada en entropía
// y densidad (proporción de 1s).
BitVectorDecision ElegirEstructuraBitVector(const std::vector<uint8_t>& bits);

// Construye y retorna un BitVector (Plain/SD/RRR) según la decisión anterior.
// Si `outDecision` no es nullptr, también devuelve la decisión tomada.
std::unique_ptr<BitVector> ConstruirBitVectorAuto(
	const std::vector<uint8_t>& bits,
	BitVectorDecision* outDecision = nullptr);

sdsl::int_vector<> build_bwt(std::string filename);

template <typename SDSL_structure>
int count_pattern(SDSL_structure &wavelet,std::map<uint8_t, int> &C,std::string &pattern,int size){
    if(pattern.empty()) return 0;

    int n = size;

    unsigned char c = pattern.back();

    auto it = C.find(c);

    if(it == C.end())
        return 0;

    int sp = it->second;

    auto next = C.upper_bound(c);

    int ep =
        (next == C.end())
        ? n
        : next->second;

    for(int i = pattern.size()-2; i >= 0; --i)
    {
        c = (unsigned char)pattern[i];

        auto it2 = C.find(c);

        if(it2 == C.end())
            return 0;
            
        std::cout << "sp: " << sp << ", ep: " << ep << ", c: " << c << std::endl;
        std::cout << "rank1: " << wavelet.rank(sp,c) << ", rank2: " << wavelet.rank(ep,c) << std::endl;    
        sp = it2->second + wavelet.rank(sp,c);
        ep = it2->second + wavelet.rank(ep,c);

        if(sp >= ep)
            return 0;
    }

    return ep - sp;

}

template <typename SDSL_structure>
int count_pattern_sdsl(SDSL_structure &wavelet,std::map<uint8_t, int> &C,std::string &pattern,int size){
    if(pattern.empty()) return 0;

    int n = size;

    unsigned char c = pattern.back();

    auto it = C.find(c);

    if(it == C.end())
        return 0;

    int sp = it->second;

    auto next = C.upper_bound(c);

    int ep =
        (next == C.end())
        ? n
        : next->second;

    for(int i = pattern.size()-2; i >= 0; --i)
    {
        c = (unsigned char)pattern[i];

        auto it2 = C.find(c);

        if(it2 == C.end())
            return 0;
        sp = it2->second + wavelet.rank(sp+1,c);
        ep = it2->second + wavelet.rank(ep+1,c);

        if(sp >= ep)
            return 0;
    }

    return ep - sp;

}
size_t file_size(const std::string& filename);
