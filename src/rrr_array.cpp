#include "bitvector/rrr_array.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

// ─── Utilidades internas ────────────────────────────────────────────────────

namespace {

static constexpr uint8_t kBlockSize = 15;

static inline uint32_t popcount16(uint16_t x) {
#if defined(_MSC_VER)
    return static_cast<uint32_t>(__popcnt(static_cast<unsigned int>(x)));
#else
    return static_cast<uint32_t>(__builtin_popcount(static_cast<unsigned int>(x)));
#endif
}

// Tabla de combinatorios C(n,k) para n,k <= 15 — inicializada una sola vez
static uint32_t comb_table[16][16];

static const uint32_t* init_comb_table() {
    for (uint32_t i = 0; i <= 15; ++i) {
        comb_table[i][0] = 1;
        comb_table[i][i] = 1;
        for (uint32_t j = 1; j < i; ++j)
            comb_table[i][j] = comb_table[i-1][j-1] + comb_table[i-1][j];
    }
    return &comb_table[0][0];
}

// Garantiza inicialización thread-safe
static const uint32_t* comb_table_init = init_comb_table();

static inline uint32_t comb(uint32_t n, uint32_t k) {
    if (k > n) return 0;
    return comb_table[n][k];
}

// Índice colex de una máscara con exactamente k bits prendidos
static inline uint16_t rank_colex(uint16_t mask, uint8_t k) {
    if (k == 0) return 0;
    uint16_t r     = 0;
    uint8_t  count = 0;
    for (uint8_t pos = 0; pos < kBlockSize; ++pos) {
        if ((mask >> pos) & 1U) {
            r = static_cast<uint16_t>(r + comb(pos, ++count));
            if (count == k) break;
        }
    }
    return r;
}

// Inversa: reconstruir máscara desde clase k y offset rank
static inline uint16_t unrank_colex(uint8_t k, uint16_t rank) {
    if (k == 0) return 0;  // única combinación válida es máscara vacía
    uint16_t mask      = 0;
    uint32_t r         = rank;
    uint8_t  remaining = k;
    uint8_t  pos       = kBlockSize;
    while (remaining > 0) {
        uint8_t p = static_cast<uint8_t>(pos - 1);
        for (;; --p) {
            const uint32_t c = comb(p, remaining);
            if (c <= r) {
                r    -= c;
                mask |= static_cast<uint16_t>(uint16_t{1} << p);
                pos   = p;
                --remaining;
                break;
            }
            if (p == 0) break; // no debería ocurrir con inputs válidos
        }
    }
    return mask;
}

} // namespace

// ─── Implementación RRRArray ────────────────────────────────────────────────

uint16_t RRRArray::block_mask(uint64_t block_index) const {
    const auto& b = blocks[static_cast<std::size_t>(block_index)];
    return unrank_colex(b.class_c, b.offset);
}

// rank1(i): cantidad de 1s en [0, i)
uint64_t RRRArray::rank1(uint64_t i) const {
    if (i == 0)      return 0;
    if (i >= nbits)  return ones;

    const uint64_t block_index = i / kBlockSize;
    const uint32_t offset      = static_cast<uint32_t>(i % kBlockSize);

    uint64_t acc = partial_sums[static_cast<std::size_t>(block_index)];
    if (offset == 0) return acc;

    const uint16_t mask   = block_mask(block_index);
    const uint16_t prefix = mask & static_cast<uint16_t>((uint16_t{1} << offset) - 1);
    return acc + popcount16(prefix);
}

// ─── Build templatizado ──────────────────────────────────────────────────────

template<typename BitSeq>
void RRRArray::build_impl(const BitSeq& bits) {
    nbits = static_cast<uint64_t>(bits.size());
    const uint64_t block_count = (nbits + kBlockSize - 1) / kBlockSize;

    blocks.assign(static_cast<std::size_t>(block_count), RRRBlock{0, 0});
    partial_sums.assign(static_cast<std::size_t>(block_count + 1), 0);

    uint32_t acc = 0;
    for (uint64_t bi = 0; bi < block_count; ++bi) {
        partial_sums[static_cast<std::size_t>(bi)] = acc;
        uint16_t mask = 0;
        for (uint8_t t = 0; t < kBlockSize; ++t) {
            const uint64_t idx = bi * kBlockSize + t;
            if (idx < nbits && bits[static_cast<std::size_t>(idx)])
                mask |= static_cast<uint16_t>(uint16_t{1} << t);
        }
        const uint8_t  c   = static_cast<uint8_t>(popcount16(mask));
        const uint16_t off = rank_colex(mask, c);
        blocks[static_cast<std::size_t>(bi)] = RRRBlock{c, off};
        acc += c;
    }
    partial_sums[static_cast<std::size_t>(block_count)] = acc;
    ones = acc;
}

//void RRRArray::build(const std::vector<uint8_t>&    bits) { build_impl(bits); }
void RRRArray::build(const std::vector<uint8_t>& bits) { build_impl(bits); }

// ─── Operaciones públicas ────────────────────────────────────────────────────

uint8_t RRRArray::access(uint64_t i) const {
    if (i >= nbits) return false;
    const uint16_t mask = block_mask(i / kBlockSize);
    return (mask >> (i % kBlockSize)) & 1U;
}

uint64_t RRRArray::rank(uint8_t bit, uint64_t i) const {
    if (i > nbits) i = nbits;
    const uint64_t r1 = rank1(i);
    return bit ? r1 : (i - r1);
}

uint64_t RRRArray::select(uint8_t bit, uint64_t j) const {
    if (j == 0) return nbits;

    if (bit) {
        if (j > ones) return nbits;
        const auto it = std::upper_bound(
            partial_sums.begin(), partial_sums.end(),
            static_cast<uint32_t>(j - 1));
        const uint64_t bi = static_cast<uint64_t>(
            std::distance(partial_sums.begin(), it) - 1);
        const uint64_t k = j - partial_sums[static_cast<std::size_t>(bi)];
        assert(k >= 1);

        uint16_t w  = block_mask(bi);
        uint64_t kk = k;
        while (kk > 1) { w &= static_cast<uint16_t>(w - 1); --kk; }

#if defined(_MSC_VER)
        unsigned long idx = 0;
        _BitScanForward(&idx, w);
        const uint32_t bitpos = static_cast<uint32_t>(idx);
#else
        const uint32_t bitpos = static_cast<uint32_t>(
            __builtin_ctz(static_cast<unsigned int>(w)));
#endif
        const uint64_t pos = bi * kBlockSize + bitpos;
        return (pos < nbits) ? pos : nbits;
    }

    // select0: búsqueda binaria sobre rank0
    const uint64_t zeros = nbits - ones;
    if (j > zeros) return nbits;
    uint64_t lo = 0, hi = nbits;
    while (lo < hi) {
        const uint64_t mid  = lo + (hi - lo) / 2;
        const uint64_t cnt0 = rank(false, mid + 1);
        if (cnt0 >= j) hi = mid;
        else           lo = mid + 1;
    }
    return lo;
}

uint64_t RRRArray::size() const { return nbits; }

uint64_t RRRArray::size_in_bytes() const {
    return sizeof(RRRArray)
         + blocks.size()        * sizeof(RRRBlock)   // 3 bytes c/u con packed
         + partial_sums.size()  * sizeof(uint32_t);
}