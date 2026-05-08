#ifndef RRR_ARRAY_HPP
#define RRR_ARRAY_HPP

#include "bitvector.hpp"
#include <cstdint>
#include <vector>

class RRRArray : public BitVector {
private:
    static constexpr uint8_t kBlockSize = 15;

    struct __attribute__((packed)) RRRBlock {
        uint8_t  class_c;   // 4 bits útiles (0–15)
        uint16_t offset;    // hasta 13 bits útiles
    };                      // 3 bytes con packed, vs 4 sin él

    std::vector<RRRBlock> blocks;
    std::vector<uint32_t> partial_sums;
    uint64_t nbits = 0;
    uint64_t ones  = 0;

    uint16_t block_mask(uint64_t block_index) const;
    uint64_t rank1(uint64_t i) const;

    // Build interno templatizado para aceptar bool y uint8_t
    template<typename BitSeq>
    void build_impl(const BitSeq& bits);

public:
    void build(const std::vector<bool>&    bits);
    void build(const std::vector<uint8_t>& bits);

    bool     access(uint64_t i)          const override;
    uint64_t rank  (bool bit, uint64_t i) const override;
    uint64_t select(bool bit, uint64_t j) const override;
    uint64_t size  ()                    const override;
    uint64_t size_in_bytes()             const override; // renombrado desde bytes()
};

#endif // RRR_ARRAY_HPP