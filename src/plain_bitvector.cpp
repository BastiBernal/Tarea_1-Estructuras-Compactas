#include "bitvector/plain_bitvector.hpp"

#include <algorithm>
#include <cstddef>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace {
	static inline uint32_t popcount64(uint64_t x) {
#if defined(_MSC_VER)
		return static_cast<uint32_t>(__popcnt64(x));
#else
		return static_cast<uint32_t>(__builtin_popcountll(x));
#endif
	}

	// Retorna la posición (0..63) del k-ésimo 1 (k es 1-based) en `word`.
	// Precondición: 1 <= k <= popcount(word).
	static inline uint32_t select1_in_word(uint64_t word, uint64_t k) {
		while (k > 1) {
			word &= (word - 1);
			--k;
		}
#if defined(_MSC_VER)
		unsigned long idx = 0;
		_BitScanForward64(&idx, word);
		return static_cast<uint32_t>(idx);
#else
		return static_cast<uint32_t>(__builtin_ctzll(word));
#endif
	}
}

uint64_t PlainBitVector::rank1(uint64_t i) const {
	if (i > nbits) i = nbits;
	const uint64_t word_index = i / 64;
	const uint32_t bit_offset = static_cast<uint32_t>(i % 64);
	const uint64_t word_count = data.size();
	if (word_index >= word_count) {
		return ones;
	}
	const uint64_t sb = word_index / 8;
	uint64_t acc = static_cast<uint64_t>(superblocks[static_cast<std::size_t>(sb)])
		+ static_cast<uint64_t>(blocks[static_cast<std::size_t>(word_index)]);
	if (bit_offset == 0) {
		return acc;
	}
	const uint64_t mask = (bit_offset == 64) ? ~uint64_t{0} : ((uint64_t{1} << bit_offset) - 1);
	acc += popcount64(data[static_cast<std::size_t>(word_index)] & mask);
	return acc;
}

uint64_t PlainBitVector::select1(uint64_t j) const {
	if (j == 0 || j > ones) return nbits;

	// Buscar superblock
	const uint64_t target = j;
	const auto it = std::upper_bound(superblocks.begin(), superblocks.end(), static_cast<uint32_t>(target - 1));
	uint64_t sb = (it == superblocks.begin()) ? 0 : static_cast<uint64_t>(std::distance(superblocks.begin(), it) - 1);
	uint64_t base_ones = superblocks[static_cast<std::size_t>(sb)];
	uint64_t word_start = sb * 8;
	uint64_t word_end = std::min<uint64_t>(word_start + 8, data.size());

	for (uint64_t wi = word_start; wi < word_end; ++wi) {
		const uint32_t wpop = popcount64(data[static_cast<std::size_t>(wi)]);
		if (base_ones + wpop >= target) {
			const uint64_t k = target - base_ones; // 1-based dentro de la palabra
			const uint32_t bitpos = select1_in_word(data[static_cast<std::size_t>(wi)], k);
			const uint64_t idx = wi * 64 + bitpos;
			return (idx < nbits) ? idx : nbits;
		}
		base_ones += wpop;
	}

	return nbits;
}

uint64_t PlainBitVector::select0(uint64_t j) const {
	const uint64_t zeros = nbits - ones;
	if (j == 0 || j > zeros) return nbits;

	uint64_t seen = 0;
	for (uint64_t wi = 0; wi < data.size(); ++wi) {
		uint64_t word = ~data[static_cast<std::size_t>(wi)];
		if (wi + 1 == data.size() && (nbits % 64) != 0) {
			const uint32_t used = static_cast<uint32_t>(nbits % 64);
			const uint64_t used_mask = (uint64_t{1} << used) - 1;
			word &= used_mask;
		}
		const uint32_t wpop = popcount64(word);
		if (seen + wpop >= j) {
			const uint64_t k = j - seen;
			const uint32_t bitpos = select1_in_word(word, k);
			const uint64_t idx = wi * 64 + bitpos;
			return (idx < nbits) ? idx : nbits;
		}
		seen += wpop;
	}
	return nbits;
}

void PlainBitVector::build(const std::vector<bool>& bits) {
	nbits = static_cast<uint64_t>(bits.size());
	const uint64_t word_count = (nbits + 63) / 64;
	data.assign(static_cast<std::size_t>(word_count), 0);

	for (uint64_t i = 0; i < nbits; ++i) {
		if (bits[static_cast<std::size_t>(i)]) {
			data[static_cast<std::size_t>(i / 64)] |= (uint64_t{1} << (i % 64));
		}
	}

	// Construir estructuras de rank (superblock=512 bits, block=64 bits)
	const uint64_t superblock_count = (word_count + 7) / 8;
	superblocks.assign(static_cast<std::size_t>(superblock_count + 1), 0);
	blocks.assign(static_cast<std::size_t>(word_count), 0);

	uint64_t acc = 0;
	for (uint64_t sb = 0; sb < superblock_count; ++sb) {
		superblocks[static_cast<std::size_t>(sb)] = static_cast<uint32_t>(acc);
		const uint64_t wstart = sb * 8;
		const uint64_t wend = std::min<uint64_t>(wstart + 8, word_count);
		for (uint64_t wi = wstart; wi < wend; ++wi) {
			blocks[static_cast<std::size_t>(wi)] = static_cast<uint16_t>(acc - superblocks[static_cast<std::size_t>(sb)]);
			acc += popcount64(data[static_cast<std::size_t>(wi)]);
		}
	}
	superblocks[static_cast<std::size_t>(superblock_count)] = static_cast<uint32_t>(acc);
	ones = acc;
}

bool PlainBitVector::access(uint64_t i) const {
	if (i >= nbits) return false;
	return (data[static_cast<std::size_t>(i / 64)] >> (i % 64)) & 1ULL;
}

uint64_t PlainBitVector::rank(bool bit, uint64_t i) const {
	if (i > nbits) i = nbits;
	const uint64_t r1 = rank1(i);
	return bit ? r1 : (i - r1);
}

uint64_t PlainBitVector::select(bool bit, uint64_t j) const {
	return bit ? select1(j) : select0(j);
}

uint64_t PlainBitVector::size() const {
	return nbits;
}

uint64_t PlainBitVector::bytes() const {
	uint64_t total = sizeof(PlainBitVector);
	total += static_cast<uint64_t>(data.capacity()) * sizeof(uint64_t);
	total += static_cast<uint64_t>(superblocks.capacity()) * sizeof(uint32_t);
	total += static_cast<uint64_t>(blocks.capacity()) * sizeof(uint16_t);
	return total;
}
