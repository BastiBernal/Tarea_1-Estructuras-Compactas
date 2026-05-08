#ifndef SD_ARRAY_HPP
#define SD_ARRAY_HPP

#include "bitvector.hpp"
#include <cstdint>
#include <vector>

class SDArray : public BitVector {
private:
    std::vector<uint64_t> high_bits;
    std::vector<uint64_t> low_bits;
    uint8_t  low_bit_length = 0;
    uint64_t nbits          = 0;
    uint64_t ones           = 0;
    uint64_t high_nbits     = 0;
    std::vector<uint32_t> high_superblocks;
    std::vector<uint16_t> high_blocks;

    uint64_t high_rank1  (uint64_t i) const;
    uint64_t high_select1(uint64_t j) const;
    uint64_t get_low     (uint64_t idx) const;
    uint64_t value_at    (uint64_t k)   const;
    uint64_t rank1       (uint64_t i)   const;

    template<typename BitSeq>
    void build_impl(const BitSeq& bits);

public:
    //void build(const std::vector<uint8_t>&    bits);
    void build(const std::vector<uint8_t>& bits);

    uint8_t     access(uint64_t i)           const override;
    uint64_t rank  (uint8_t bit, uint64_t i) const override;
    uint64_t select(uint8_t bit, uint64_t j) const override;
    uint64_t size  ()                     const override;
    uint64_t size_in_bytes()              const override;
};

#endif // SD_ARRAY_HPP