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
BitVectorDecision ElegirEstructuraBitVector(const std::vector<bool>& bits);

// Construye y retorna un BitVector (Plain/SD/RRR) según la decisión anterior.
// Si `outDecision` no es nullptr, también devuelve la decisión tomada.
std::unique_ptr<BitVector> ConstruirBitVectorAuto(
	const std::vector<bool>& bits,
	BitVectorDecision* outDecision = nullptr);

sdsl::int_vector<> build_bwt(std::string filename);

template <typename SDSL_structure>
int count_pattern(SDSL_structure &wavelet, std::map<int, int> &C, std::string &pattern){
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
