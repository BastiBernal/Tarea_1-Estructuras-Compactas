
#pragma once

#include "bitvector/bitvector.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// Wavelet Tree para secuencias de bytes (std::string interpretado como unsigned char).
//
// Semántica (consistente con BitVector):
// - rank(c, i) = # de ocurrencias de c en el prefijo [0, i)
// - select(c, j) con j 1-based: posición (0-based) del j-ésimo c, o size() si no existe
class WaveletTree {
public:
	WaveletTree() = default;
	explicit WaveletTree(const std::string& text);

	void build(const std::string& text);

	uint64_t size() const { return n_; }
	// Aproximación del uso de memoria en RAM (en bytes), incluyendo bitvectors.
	// Nota: no incluye overhead del allocator.
	uint64_t size_bytes() const;
	bool empty() const { return n_ == 0; }

	char access(uint64_t i) const;
	uint64_t rank(char c, uint64_t i) const;
	uint64_t select(char c, uint64_t j) const;

private:
	struct Node {
		unsigned char lo = 0;
		unsigned char hi = 0;
		uint64_t n = 0;
		std::unique_ptr<BitVector> bv;
		std::unique_ptr<Node> left;
		std::unique_ptr<Node> right;
	};

	std::unique_ptr<Node> root_;
	uint64_t n_ = 0;

	static std::unique_ptr<Node> build_rec(
		const std::vector<unsigned char>& data,
		unsigned char lo,
		unsigned char hi);
};

// Implementación "pointerless/levelwise": en vez de guardar un BitVector por nodo,
// guarda 1 BitVector por nivel con la concatenación de los nodos de ese nivel.
//
// Convención de bits (como pediste):
// - 1 => hijo izquierdo
// - 0 => hijo derecho
//
// Además precomputa, por nodo, el rango [start,end) dentro del bitvector del nivel
// para poder ubicar inmediatamente dónde termina cada nodo.
class WaveletTreeBinary {
public:
	WaveletTreeBinary() = default;
	explicit WaveletTreeBinary(const std::string& text);

	void build(const std::string& text);

	uint64_t size() const { return n_; }
	// Aproximación del uso de memoria en RAM (en bytes), incluyendo bitvectors.
	// Nota: no incluye overhead del allocator.
	uint64_t size_bytes() const;
	bool empty() const { return n_ == 0; }

	char access(uint64_t i) const;
	uint64_t rank(char c, uint64_t i) const;
	uint64_t select(char c, uint64_t j) const;

private:
	struct NodeInfo {
		unsigned char lo = 0;
		unsigned char hi = 0;
		unsigned char mid = 0;
		uint64_t start = 0; // inicio (incl.) en el bitvector del nivel
		uint64_t end = 0;   // fin (excl.)  en el bitvector del nivel
		int32_t left = -1;  // índice de nodo hijo (nivel+1), o -1 si es hoja/ausente
		int32_t right = -1;
		uint64_t rank1_before = 0; // #1s en [0, start)
		uint64_t rank0_before = 0; // #0s en [0, start)
	};

	struct Level {
		std::unique_ptr<BitVector> bv;
		std::vector<NodeInfo> nodes;
	};

	std::vector<Level> levels_;
	uint64_t n_ = 0;
	unsigned char min_c_ = 0;
	unsigned char max_c_ = 0;

	static int32_t build_rec(
		std::vector<std::vector<bool>>& level_bits,
		std::vector<std::vector<NodeInfo>>& level_nodes,
		const std::vector<unsigned char>& data,
		unsigned char lo,
		unsigned char hi,
		uint32_t depth);
};

