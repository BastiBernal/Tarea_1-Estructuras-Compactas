#pragma once

#include "bitvector/bitvector.hpp"
#include "sdsl/int_vector.hpp"

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <memory>

// Implementación de un wavelet tree sin punteros usando el vector de bitvectors
class WaveletTreePointerless
{
private:

    struct Level
    {
        std::unique_ptr<BitVector> bv;
        uint64_t zero_count;
    };

    std::vector<Level> levels_;

    uint64_t n_ = 0;

    uint32_t sigma_ = 0;
    uint32_t height_ = 0;

    std::vector<int64_t> id_to_symbol_;
    std::unordered_map<int64_t,uint32_t> symbol_to_id_;

public:

    void build(const sdsl::int_vector<>& values);

    uint64_t rank(int64_t i, uint64_t c) const;

    int64_t access(uint64_t i) const;

    uint64_t size_in_bytes() const;
};