#pragma once
#include "wavelet/wavelet_binary.hpp"
#include "utils-p/utils.hpp"


void WaveletTreePointerless::build(const sdsl::int_vector<>& values) {
    levels_.clear();
    n_ = values.size();
    if (n_ == 0) return;

    // 1. Compactar alfabeto
    std::vector<int64_t> symbols(values.begin(), values.end());
    std::sort(symbols.begin(), symbols.end());
    symbols.erase(std::unique(symbols.begin(), symbols.end()), symbols.end());

    sigma_ = static_cast<uint32_t>(symbols.size());
    id_to_symbol_ = symbols;
    symbol_to_id_.clear();
    for (uint32_t i = 0; i < sigma_; ++i)
        symbol_to_id_[symbols[i]] = i;

    // 2. Calcular altura (mínimo 1 si hay símbolos)
    height_ = 0;
    while ((1ULL << height_) < sigma_) ++height_;
    if (height_ == 0 && sigma_ > 0) height_ = 1;

    // 3. Construcción Levelwise
    std::vector<uint32_t> cur(n_);
    for (uint64_t i = 0; i < n_; ++i)
        cur[i] = symbol_to_id_[values[i]];

    levels_.resize(height_);
    std::vector<uint32_t> next(n_);

    for (uint32_t level = 0; level < height_; ++level) {
        std::vector<uint8_t> bits(n_);
        uint32_t shift = height_ - level - 1;
        uint64_t zero_pos = 0;

        // Contar ceros para el offset del siguiente nivel
        for (uint32_t x : cur) {
            if (!((x >> shift) & 1ULL)) ++zero_pos;
        }
        levels_[level].zero_count = zero_pos;

        uint64_t z = 0;
        uint64_t o = zero_pos;
        for (uint64_t i = 0; i < n_; ++i) {
            bool bit = (cur[i] >> shift) & 1ULL;
            bits[i] = bit;
            if (!bit) next[z++] = cur[i];
            else      next[o++] = cur[i];
        }

        levels_[level].bv = ConstruirBitVectorAuto(bits);
        cur.swap(next);
    }
}

uint64_t WaveletTreePointerless::rank(int64_t c, uint64_t i) const {
    // Si i es 0, no hay elementos que contar
    if (i == 0) return 0;
    // Si i se pasa del tamaño, lo limitamos a n_
    if (i > n_) i = n_;

    auto it = symbol_to_id_.find(c);
    if (it == symbol_to_id_.end()) return 0;

    uint32_t symbol_id = it->second;
    uint64_t pos = i;
    uint64_t start = 0; // Para calcular el offset del símbolo

    for (uint32_t level = 0; level < height_; ++level) {
        //
        const auto& lvl = levels_[level];
        bool bit = (symbol_id >> (height_ - level - 1)) & 1ULL;
        // Calculamos cuántos 1s hay antes de 'pos' y 'start' para ajustar la posición
        uint64_t ones_before_pos   = lvl.bv->rank(true, pos);
        uint64_t ones_before_start = lvl.bv->rank(true, start);

        if (bit) { // Si el bit es 1, nos movemos a la parte de unos
            pos   = lvl.zero_count + ones_before_pos;
            start = lvl.zero_count + ones_before_start;
        } else { // Si el bit es 0, nos movemos a la parte de ceros
            pos   = pos   - ones_before_pos;
            start = start - ones_before_start;
        }
    }

    return pos - start;
}
    
    // Para obtener el rank real, restamos la posición donde empezaría el símbolo con rank(c, 0)
    uint64_t start_offset = 0;
    uint64_t temp_pos = 0; // Calculamos el 'left' original
    for (uint32_t level = 0; level < height_; ++level) {
        const auto& lvl = levels_[level];
        bool bit = (symbol_id >> (height_ - level - 1)) & 1ULL;
        if (bit) temp_pos = lvl.zero_count + lvl.bv->rank(true, temp_pos);
        else     temp_pos = temp_pos - lvl.bv->rank(true, temp_pos);
    }

    return pos - temp_pos;
}

int64_t WaveletTreePointerless::access(uint64_t i) const {
    if (i >= n_) return -1; // -1 indica fuera de rango

    uint32_t symbol_id = 0;
    uint64_t pos = i;

    for (uint32_t level = 0; level < height_; ++level) {
        const auto& lvl = levels_[level];
        bool bit = lvl.bv->access(pos);

        symbol_id <<= 1;
        if (bit) {
            symbol_id |= 1ULL;
            pos = lvl.zero_count + lvl.bv->rank(true, pos);
        } else {
            pos = pos - lvl.bv->rank(true, pos);
        }
    }

    return id_to_symbol_[symbol_id];
}

uint64_t WaveletTreePointerless::size_in_bytes() const {
    uint64_t total = 0;

    // 1. Variables escalares
    total += sizeof(n_);
    total += sizeof(sigma_);
    total += sizeof(height_);

    // 2. Bitvectores por nivel
    for (const auto& lvl : levels_) {
        total += sizeof(lvl.zero_count);
        if (lvl.bv) {
            total += sdsl::size_in_bytes(*lvl.bv);
        }
    }
    // Overhead del vector levels_ en sí (sin contar el contenido ya contado)
    total += sizeof(levels_);

    // 3. Mapas de alfabeto
    // id_to_symbol_: vector<int64_t>
    total += sizeof(id_to_symbol_);
    total += id_to_symbol_.size() * sizeof(int64_t);

    // symbol_to_id_: unordered_map<int64_t, uint32_t>
    // size_in_bytes exacto no es trivial; aproximación conservadora:
    total += sizeof(symbol_to_id_);
    total += symbol_to_id_.size() * (sizeof(int64_t) + sizeof(uint32_t));

    return total;
}