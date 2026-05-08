#include "bitvector/sd_array.hpp"
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
} // namespace

// ─── Rank sobre high_bits ────────────────────────────────────────────────────

uint64_t SDArray::high_rank1(uint64_t i) const {
    if (i == 0)          return 0;
    if (i >= high_nbits) return ones;

    const uint64_t wi  = i / 64;
    const uint32_t off = static_cast<uint32_t>(i % 64);

    const uint64_t sb  = wi / 8;
    uint64_t acc = static_cast<uint64_t>(high_superblocks[sb])
                 + static_cast<uint64_t>(high_blocks[wi]);
    if (off == 0) return acc;

    const uint64_t mask = (uint64_t{1} << off) - 1; // off nunca es 64
    return acc + popcount64(high_bits[wi] & mask);
}

// ─── Select sobre high_bits usando superblocks ───────────────────────────────

uint64_t SDArray::high_select1(uint64_t j) const {
    if (j == 0 || j > ones) return high_nbits;

    // 1. Localizar superbloque con búsqueda binaria
    uint64_t sb_lo = 0;
    uint64_t sb_hi = static_cast<uint64_t>(high_superblocks.size()) - 1;
    while (sb_lo + 1 < sb_hi) {
        const uint64_t mid = sb_lo + (sb_hi - sb_lo) / 2;
        if (high_superblocks[mid] < j) sb_lo = mid;
        else sb_hi = mid;
    }

    // 2. Escanear words del superbloque
    const uint64_t w_start = sb_lo * 8;
    const uint64_t w_end   = std::min<uint64_t>(w_start + 8,
                                 static_cast<uint64_t>(high_bits.size()));
    uint64_t acc = high_superblocks[sb_lo];
    uint64_t wi  = w_start;
    while (wi + 1 < w_end) {
        const uint64_t cnt = popcount64(high_bits[wi]);
        if (acc + cnt >= j) break;
        acc += cnt;
        ++wi;
    }

    // 3. Localizar bit dentro del word
    uint64_t word   = high_bits[wi];
    uint64_t remain = j - acc;
    while (remain > 1) { word &= word - 1; --remain; }

#if defined(_MSC_VER)
    unsigned long idx = 0;
    _BitScanForward64(&idx, word);
    return wi * 64 + static_cast<uint64_t>(idx);
#else
    return wi * 64 + static_cast<uint64_t>(__builtin_ctzll(word));
#endif
}

// ─── Low bits ────────────────────────────────────────────────────────────────

uint64_t SDArray::get_low(uint64_t idx) const {
    if (low_bit_length == 0) return 0;
    const uint64_t bitpos = idx * static_cast<uint64_t>(low_bit_length);
    const uint64_t wi     = bitpos / 64;
    const uint32_t off    = static_cast<uint32_t>(bitpos % 64);
    const uint64_t mask   = (uint64_t{1} << low_bit_length) - 1;
    uint64_t value = low_bits[wi] >> off;
    if (off + low_bit_length > 64 && wi + 1 < low_bits.size())
        value |= low_bits[wi + 1] << (64 - off);
    return value & mask;
}

uint64_t SDArray::value_at(uint64_t k) const {
    const uint64_t pos = high_select1(k + 1);
    return ((pos - k) << low_bit_length) | get_low(k);
}

// ─── Rank1 sobre el bitvector completo ──────────────────────────────────────

uint64_t SDArray::rank1(uint64_t i) const {
    if (i == 0)     return 0;
    if (i >= nbits) return ones;
    // Primer índice k tal que value_at(k) >= i → k es la cantidad de 1s en [0,i)
    uint64_t lo = 0, hi = ones;
    while (lo < hi) {
        const uint64_t mid = lo + (hi - lo) / 2;
        if (value_at(mid) < i) lo = mid + 1;
        else                   hi = mid;
    }
    return lo;
}

// ─── Build ───────────────────────────────────────────────────────────────────

template<typename BitSeq>
void SDArray::build_impl(const BitSeq& bits) {
    nbits = static_cast<uint64_t>(bits.size());
    std::vector<uint64_t> positions;
    positions.reserve(nbits);
    for (uint64_t i = 0; i < nbits; ++i)
        if (bits[static_cast<std::size_t>(i)]) positions.push_back(i);

    ones           = static_cast<uint64_t>(positions.size());
    low_bit_length = 0;
    high_nbits     = 0;
    low_bits.clear(); high_bits.clear();
    high_superblocks.clear(); high_blocks.clear();

    if (ones == 0) return;

    // l = floor(log2(n/m))
    uint8_t l = 0;
    uint64_t ratio = nbits / ones;
    while ((l + 1) < 63 && (uint64_t{1} << (l + 1)) <= ratio) ++l;
    low_bit_length = l;

    high_nbits = (nbits >> low_bit_length) + ones + 1;
    high_bits.assign((high_nbits + 63) / 64, 0);

    const uint64_t low_words = (ones * low_bit_length + 63) / 64;
    low_bits.assign(low_words == 0 ? 0 : low_words, 0);

    const uint64_t low_mask = (low_bit_length == 0) ? 0
                            : ((uint64_t{1} << low_bit_length) - 1);

    for (uint64_t idx = 0; idx < ones; ++idx) {
        const uint64_t v  = positions[idx];
        const uint64_t hi = v >> low_bit_length;
        const uint64_t lo = v & low_mask;          // low_mask==0 → lo==0

        const uint64_t hp = hi + idx;
        high_bits[hp / 64] |= uint64_t{1} << (hp % 64);

        if (low_bit_length != 0) {
            const uint64_t bitpos = idx * low_bit_length;
            const uint64_t wi     = bitpos / 64;
            const uint32_t off    = static_cast<uint32_t>(bitpos % 64);
            low_bits[wi] |= lo << off;
            if (off + low_bit_length > 64)
                low_bits[wi + 1] |= lo >> (64 - off);
        }
    }

    // Rank support: superblock=512 bits (8 words de 64b)
    const uint64_t hw    = static_cast<uint64_t>(high_bits.size());
    const uint64_t sb_n  = (hw + 7) / 8;
    high_superblocks.assign(sb_n + 1, 0);
    high_blocks.assign(hw, 0);

    uint64_t acc = 0;
    for (uint64_t sb = 0; sb < sb_n; ++sb) {
        high_superblocks[sb] = static_cast<uint32_t>(acc);
        const uint64_t wend  = std::min(sb * 8 + 8, hw);
        for (uint64_t wi = sb * 8; wi < wend; ++wi) {
            high_blocks[wi] = static_cast<uint16_t>(
                acc - high_superblocks[sb]);
            acc += popcount64(high_bits[wi]);
        }
    }
    high_superblocks[sb_n] = static_cast<uint32_t>(acc);
}

//void SDArray::build(const std::vector<uint8_t>&    bits) { build_impl(bits); }
void SDArray::build(const std::vector<uint8_t>& bits) { build_impl(bits); }

// ─── Operaciones públicas ────────────────────────────────────────────────────

uint8_t SDArray::access(uint64_t i) const {
    if (i >= nbits) return false;
    const uint64_t r = rank1(i + 1);
    return r > 0 && value_at(r - 1) == i;
}

uint64_t SDArray::rank(uint8_t bit, uint64_t i) const {
    if (i > nbits) i = nbits;
    const uint64_t r1 = rank1(i);
    return bit ? r1 : (i - r1);
}

uint64_t SDArray::select(uint8_t bit, uint64_t j) const {
    if (bit) {
        if (j == 0 || j > ones) return nbits;
        return value_at(j - 1);
    }
    const uint64_t zeros = nbits - ones;
    if (j == 0 || j > zeros) return nbits;
    uint64_t lo = 0, hi = nbits;
    while (lo < hi) {
        const uint64_t mid  = lo + (hi - lo) / 2;
        if (rank(false, mid + 1) >= j) hi = mid;
        else                           lo = mid + 1;
    }
    return lo;
}

uint64_t SDArray::size() const { return nbits; }

uint64_t SDArray::size_in_bytes() const {
    return sizeof(SDArray)
         + high_bits.size()        * sizeof(uint64_t)
         + low_bits.size()         * sizeof(uint64_t)
         + high_superblocks.size() * sizeof(uint32_t)
         + high_blocks.size()      * sizeof(uint16_t);
}