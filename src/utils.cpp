#include "utils-p/utils.hpp"
#include "bitvector/plain_bitvector.hpp"
#include "bitvector/rrr_array.hpp"
#include "bitvector/sd_array.hpp"
#include "sdsl/suffix_arrays.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <unordered_map>
#include <iostream>
#include <fstream>

std::string build_pattern(std::string filename,int large){
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        std::cerr << "Error al abrir el archivo\n";
        return std::string();
    }

    std::string pattern(large, '\0');

    in.read(&pattern[0], large);
    pattern.resize(in.gcount()); // ajusta si el archivo tiene menos de N

    return pattern;
}

ResultadoEntropia CalcularEntropia(const std::string& texto) {
	ResultadoEntropia resultado;

	const std::size_t n = texto.size();
	if (n == 0) {
		return resultado;
	}

	// Pasada 1: Conteo de frecuencias
	std::unordered_map<char, std::size_t> frecuencias;
	frecuencias.reserve(texto.size());
	for (char c : texto) {
		++frecuencias[c];
	}

	// Pasada 2: Cálculo de probabilidades y Shannon
	resultado.distribucion.reserve(frecuencias.size());
	double entropia = 0.0;
	for (const auto& [caracter, conteo] : frecuencias) {
		const double p = static_cast<double>(conteo) / static_cast<double>(n);
		resultado.distribucion[caracter] = p;
		// p nunca es 0 aquí porque viene de un conteo > 0.
		entropia -= p * std::log2(p);
	}

	resultado.entropia = entropia;
	return resultado;
}

double CalcularEntropiaBits(double p1) {
	if (p1 <= 0.0 || p1 >= 1.0) return 0.0;
	const double p0 = 1.0 - p1;
	return -(p1 * std::log2(p1) + p0 * std::log2(p0));
}

BitVectorDecision ElegirEstructuraBitVector(const std::vector<bool>& bits) {
	BitVectorDecision d;
	d.n = static_cast<uint64_t>(bits.size());
	if (d.n == 0) {
		d.tipo = BitVectorTipo::Plain;
		return d;
	}

	uint64_t ones = 0;
	for (bool b : bits) if (b) ++ones;
	d.ones = ones;
	d.p1 = static_cast<double>(ones) / static_cast<double>(d.n);
	d.entropia = CalcularEntropiaBits(d.p1);

	// Heurística simple:
	// - SDArray suele convenir cuando hay muy pocos 1s (sparse) (p1 muy chico).
	// - PlainBitVector conviene cuando la entropía es alta (casi aleatorio, p~0.5).
	// - RRRArray es una opción robusta cuando la entropía es baja/moderada o la
	//   densidad es alta (p cercano a 1) y no queremos penalizar SD por ser denso.
	const double p = d.p1;
	const double H = d.entropia;

	if (p <= 0.07) {
		d.tipo = BitVectorTipo::SD;
		return d;
	}
	if (p >= 0.93) {
		d.tipo = BitVectorTipo::RRR;
		return d;
	}
	if (H >= 0.98) {
		d.tipo = BitVectorTipo::Plain;
		return d;
	}

	// Zona intermedia: si aún es relativamente sparse, SD; en caso contrario RRR.
	if (p <= 0.15) {
		d.tipo = BitVectorTipo::SD;
	} else {
		d.tipo = BitVectorTipo::RRR;
	}
	return d;
}

std::unique_ptr<BitVector> ConstruirBitVectorAuto(
	const std::vector<bool>& bits,
	BitVectorDecision* outDecision) {
	BitVectorDecision d = ElegirEstructuraBitVector(bits);
	if (outDecision) *outDecision = d;

	switch (d.tipo) {
		case BitVectorTipo::Plain: {
			auto bv = std::make_unique<PlainBitVector>();
			bv->build(bits);
			return bv;
		}
		case BitVectorTipo::SD: {
			auto bv = std::make_unique<SDArray>();
			bv->build(bits);
			return bv;
		}
		case BitVectorTipo::RRR: {
			auto bv = std::make_unique<RRRArray>();
			bv->build(bits);
			return bv;
		}
	}

	// Fallback (no debería ocurrir)
	auto bv = std::make_unique<PlainBitVector>();
	bv->build(bits);
	return bv;
}

sdsl::int_vector<> build_bwt(std::string filename){
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