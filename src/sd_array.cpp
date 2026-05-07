#include "sd_array.hpp"

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
}

uint64_t SDArray::high_rank1(uint64_t i) const {
	if (i > high_nbits) i = high_nbits;
	const uint64_t word_index = i / 64;
	const uint32_t bit_offset = static_cast<uint32_t>(i % 64);
	if (word_index >= high_bits.size()) {
		return ones;
	}
	const uint64_t sb = word_index / 8;
	uint64_t acc = static_cast<uint64_t>(high_superblocks[static_cast<std::size_t>(sb)])
		+ static_cast<uint64_t>(high_blocks[static_cast<std::size_t>(word_index)]);
	if (bit_offset == 0) return acc;
	const uint64_t mask = (bit_offset == 64) ? ~uint64_t{0} : ((uint64_t{1} << bit_offset) - 1);
	acc += popcount64(high_bits[static_cast<std::size_t>(word_index)] & mask);
	return acc;
}

uint64_t SDArray::high_select1(uint64_t j) const {
	// j es 1-based. Retorna posición (0-based) del j-ésimo 1 en high_bits.
	if (j == 0 || j > ones) return high_nbits;
	uint64_t lo = 0;
	uint64_t hi = high_nbits;
	while (lo < hi) {
		const uint64_t mid = lo + (hi - lo) / 2;
		const uint64_t cnt = high_rank1(mid + 1);
		if (cnt >= j) hi = mid;
		else lo = mid + 1;
	}
	return lo;
}

uint64_t SDArray::get_low(uint64_t idx) const {
	if (low_bit_length == 0) return 0;
	const uint64_t bitpos = idx * static_cast<uint64_t>(low_bit_length);
	const uint64_t word_index = bitpos / 64;
	const uint32_t offset = static_cast<uint32_t>(bitpos % 64);
	const uint64_t mask = (low_bit_length == 64) ? ~uint64_t{0} : ((uint64_t{1} << low_bit_length) - 1);
	uint64_t value = 0;
	if (word_index < low_bits.size()) {
		value = low_bits[static_cast<std::size_t>(word_index)] >> offset;
		if (offset + low_bit_length > 64 && word_index + 1 < low_bits.size()) {
			value |= (low_bits[static_cast<std::size_t>(word_index + 1)] << (64 - offset));
		}
	}
	return value & mask;
}

uint64_t SDArray::value_at(uint64_t k) const {
	// k es 0-based en [0, ones)
	const uint64_t pos = high_select1(k + 1); // hi + k
	const uint64_t hi = pos - k;
	const uint64_t lo = get_low(k);
	return (hi << low_bit_length) | lo;
}

uint64_t SDArray::rank1(uint64_t i) const {
	if (i > nbits) i = nbits;
	// Número de 1s en [0, i) => contar posiciones < i.
	uint64_t lo = 0;
	uint64_t hi = ones;
	while (lo < hi) {
		const uint64_t mid = lo + (hi - lo) / 2;
		if (value_at(mid) < i) lo = mid + 1;
		else hi = mid;
	}
	return lo;
}

void SDArray::build(const std::vector<bool>& bits) {
	nbits = static_cast<uint64_t>(bits.size());
	std::vector<uint64_t> positions;
	positions.reserve(bits.size());
	for (uint64_t i = 0; i < nbits; ++i) {
		if (bits[static_cast<std::size_t>(i)]) positions.push_back(i);
	}
	ones = static_cast<uint64_t>(positions.size());
	low_bits.clear();
	high_bits.clear();
	high_superblocks.clear();
	high_blocks.clear();
	low_bit_length = 0;
	high_nbits = 0;

	if (ones == 0) {
		return;
	}

	// Elias-Fano: l = floor(log2(n / m))
	uint64_t ratio = nbits / ones;
	uint8_t l = 0;
	while ((l + 1) < 63 && (uint64_t{1} << (l + 1)) <= ratio) {
		++l;
	}
	low_bit_length = l;

	const uint64_t high_len = (nbits >> low_bit_length) + ones + 1;
	high_nbits = high_len;
	const uint64_t high_words = (high_len + 63) / 64;
	high_bits.assign(static_cast<std::size_t>(high_words), 0);

	// Empaquetar low bits
	const uint64_t low_total_bits = ones * static_cast<uint64_t>(low_bit_length);
	const uint64_t low_words = (low_total_bits + 63) / 64;
	low_bits.assign(static_cast<std::size_t>(low_words), 0);

	const uint64_t low_mask = (low_bit_length == 0) ? 0 : ((uint64_t{1} << low_bit_length) - 1);

	for (uint64_t idx = 0; idx < ones; ++idx) {
		const uint64_t v = positions[static_cast<std::size_t>(idx)];
		const uint64_t hi = v >> low_bit_length;
		const uint64_t lo = (low_bit_length == 0) ? 0 : (v & low_mask);

		// set bit (hi + idx) en high_bits
		const uint64_t hp = hi + idx;
		high_bits[static_cast<std::size_t>(hp / 64)] |= (uint64_t{1} << (hp % 64));

		// escribir low bits
		if (low_bit_length != 0) {
			const uint64_t bitpos = idx * static_cast<uint64_t>(low_bit_length);
			const uint64_t wi = bitpos / 64;
			const uint32_t off = static_cast<uint32_t>(bitpos % 64);
			low_bits[static_cast<std::size_t>(wi)] |= (lo << off);
			if (off + low_bit_length > 64) {
				low_bits[static_cast<std::size_t>(wi + 1)] |= (lo >> (64 - off));
			}
		}
	}

	// Rank structures sobre high_bits (superblock=512b)
	const uint64_t superblock_count = (high_words + 7) / 8;
	high_superblocks.assign(static_cast<std::size_t>(superblock_count + 1), 0);
	high_blocks.assign(static_cast<std::size_t>(high_words), 0);
	uint64_t acc = 0;
	for (uint64_t sb = 0; sb < superblock_count; ++sb) {
		high_superblocks[static_cast<std::size_t>(sb)] = static_cast<uint32_t>(acc);
		const uint64_t wstart = sb * 8;
		const uint64_t wend = std::min<uint64_t>(wstart + 8, high_words);
		for (uint64_t wi = wstart; wi < wend; ++wi) {
			high_blocks[static_cast<std::size_t>(wi)] = static_cast<uint16_t>(acc - high_superblocks[static_cast<std::size_t>(sb)]);
			acc += popcount64(high_bits[static_cast<std::size_t>(wi)]);
		}
	}
	high_superblocks[static_cast<std::size_t>(superblock_count)] = static_cast<uint32_t>(acc);
	// acc == ones
}

bool SDArray::access(uint64_t i) const {
	if (i >= nbits) return false;
	const uint64_t r = rank1(i + 1);
	if (r == 0) return false;
	return value_at(r - 1) == i;
}

uint64_t SDArray::rank(bool bit, uint64_t i) const {
	if (i > nbits) i = nbits;
	const uint64_t r1 = rank1(i);
	return bit ? r1 : (i - r1);
}

uint64_t SDArray::select(bool bit, uint64_t j) const {
	if (bit) {
		if (j == 0 || j > ones) return nbits;
		return value_at(j - 1);
	}

	// select0 vía búsqueda binaria sobre posiciones usando rank0
	const uint64_t zeros = nbits - ones;
	if (j == 0 || j > zeros) return nbits;
	uint64_t lo = 0;
	uint64_t hi = nbits;
	while (lo < hi) {
		const uint64_t mid = lo + (hi - lo) / 2;
		const uint64_t cnt0 = rank(false, mid + 1);
		if (cnt0 >= j) hi = mid;
		else lo = mid + 1;
	}
	return lo;
}

uint64_t SDArray::size() const {
	return nbits;
}

uint64_t SDArray::bytes() const {
	uint64_t total = sizeof(SDArray);
	total += static_cast<uint64_t>(high_bits.capacity()) * sizeof(uint64_t);
	total += static_cast<uint64_t>(low_bits.capacity()) * sizeof(uint64_t);
	total += static_cast<uint64_t>(high_superblocks.capacity()) * sizeof(uint32_t);
	total += static_cast<uint64_t>(high_blocks.capacity()) * sizeof(uint16_t);
	return total;
}
