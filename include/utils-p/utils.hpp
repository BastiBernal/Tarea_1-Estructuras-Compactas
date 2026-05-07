#pragma once 

#include "bitvector/bitvector.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

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

