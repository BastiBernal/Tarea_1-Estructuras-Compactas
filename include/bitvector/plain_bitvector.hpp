#ifndef PLAIN_BITVECTOR_HPP
#define PLAIN_BITVECTOR_HPP

#include "bitvector.hpp"

#include <cstdint>
#include <vector>

class PlainBitVector : public BitVector {
private:
    std::vector<uint64_t> data;
    std::vector<uint32_t> superblocks;
    std::vector<uint16_t> blocks;
    uint64_t nbits = 0;
    uint64_t ones = 0;

    uint64_t rank1(uint64_t i) const;
    uint64_t select1(uint64_t j) const;
    uint64_t select0(uint64_t j) const;
public:
    void build(const std::vector<bool>& bits);
    bool access(uint64_t i) const override;
    uint64_t rank(bool bit, uint64_t i) const override;
    uint64_t select(bool bit, uint64_t j) const override;
    uint64_t size() const override;
    uint64_t bytes() const override;
};

#endif // PLAIN_BITVECTOR_HPP