#pragma once
#include "wavelet/wavelet_binary.hpp"

void WaveletTreePointerless::build(const sdsl::int_vector<>& values)
{
    levels_.clear();

    n_ = values.size();

    if (n_ == 0)
        return;

    // Compactar alfabeto

    std::vector<int64_t> symbols(values.begin(), values.end());

    std::sort(symbols.begin(), symbols.end());
    symbols.erase(std::unique(symbols.begin(), symbols.end()), symbols.end());

    sigma_ = static_cast<uint32_t>(symbols.size());

    id_to_symbol_ = symbols;

    symbol_to_id_.clear();

    for (uint32_t i = 0; i < sigma_; ++i)
        symbol_to_id_[symbols[i]] = i;

    // Altura

    height_ = 0;

    while ((1ULL << height_) < sigma_)
        ++height_;

    if (height_ == 0)
        height_ = 1;

    // Datos compactados

    std::vector<uint32_t> data(n_);

    for (uint64_t i = 0; i < n_; ++i)
        data[i] = symbol_to_id_[values[i]];

    // Construcción levelwise

    levels_.resize(height_);

    std::vector<uint32_t> cur = data;
    std::vector<uint32_t> next(n_);

    for (uint32_t level = 0; level < height_; ++level)
    {
        std::vector<bool> bits;
        bits.reserve(n_);

        uint32_t shift = height_ - level - 1;

        uint64_t zero_pos = 0;
        uint64_t one_pos = 0;

        // contar ceros

        for (uint32_t x : cur)
        {
            bool bit = (x >> shift) & 1ULL;

            if (!bit)
                ++zero_pos;
        }

        levels_[level].zero_count = zero_pos;

        one_pos = zero_pos;

        // llenar bitvector y preparar orden para el siguiente nivel

        uint64_t z = 0;
        uint64_t o = zero_pos;

        for (uint32_t x : cur) {
            bool bit = (x >> shift) & 1ULL;

            bits.push_back(bit);

            if (!bit) next[z++] = x;
            else next[o++] = x;
        }

        levels_[level].bv = ConstruirBitVectorAuto(bits);

        cur.swap(next);
    }
}

uint64_t WaveletTreePointerless::rank(int64_t c, uint64_t i) const {
    auto it = symbol_to_id_.find(c);

    if (it == symbol_to_id_.end())
        return 0;

    uint32_t symbol = it->second;

    uint64_t left = 0;
    uint64_t right = i;

    for (uint32_t level = 0; level < height_; ++level) {
        const auto& lvl = levels_[level];

        uint32_t shift = height_ - level - 1;

        bool bit = (symbol >> shift) & 1ULL;

        if (bit) {
            left =
                lvl.zero_count +
                lvl.bv->rank(true, left);

            right =
                lvl.zero_count +
                lvl.bv->rank(true, right);
        }
        else {
            left =
                left -
                lvl.bv->rank(true, left);

            right =
                right -
                lvl.bv->rank(true, right);
        }
    }

    return right - left;
}

int64_t WaveletTreePointerless::access(uint64_t i) const
{
    if (i >= n_) return 1;

    uint32_t symbol = 0;

    for (uint32_t level = 0; level < height_; ++level)
    {
        const auto& lvl = levels_[level];

        bool bit = lvl.bv->access(i);

        symbol <<= 1;

        if (bit) {
            symbol |= 1ULL;
            i = lvl.zero_count + lvl.bv->rank(true, i);
        }
        else {
            i = i - lvl.bv->rank(true, i);
        }
    }

    return id_to_symbol_[symbol];
}

uint64_t WaveletTreePointerless::size_bytes() const
{
    uint64_t size = sizeof(*this);
    for (const auto& lvl : levels_) {
        size += sizeof(Level) + lvl.bv->size_bytes();
    }
    size += id_to_symbol_.size() * sizeof(int64_t);
    size += symbol_to_id_.size() * 24; // Estimación para unordered_map
    return size;
}