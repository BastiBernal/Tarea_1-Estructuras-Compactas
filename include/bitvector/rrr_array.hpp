#ifndef RRR_ARRAY_HPP
#define RRR_ARRAY_HPP

#include "bitvector.hpp"

#include <cstdint>
#include <vector>

class RRRArray : public BitVector {
private:
    static constexpr uint8_t kBlockSize = 15;

    struct RRRBlock {
        uint8_t class_c;
        uint16_t offset;
    };
    std::vector<RRRBlock> blocks;
    std::vector<uint32_t> partial_sums;
    uint64_t nbits = 0;
    uint64_t ones = 0;

    uint16_t block_mask(uint64_t block_index) const;
    uint64_t rank1(uint64_t i) const;
public:
    void build(const std::vector<bool>& bits);
    bool access(uint64_t i) const override;
    uint64_t rank(bool bit, uint64_t i) const override;
    uint64_t select(bool bit, uint64_t j) const override;
    uint64_t size() const override;
    uint64_t bytes() const override;
};

#endif // RRR_ARRAY_HPP