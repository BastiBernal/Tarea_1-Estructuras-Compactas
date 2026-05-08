#ifndef BITVECTOR_HPP
#define BITVECTOR_HPP

#include <cstdint>

class BitVector {
public:
    virtual ~BitVector() {}
    virtual uint8_t access(uint64_t i) const = 0;
    virtual uint64_t rank(uint8_t bit, uint64_t i) const = 0;
    virtual uint64_t select(uint8_t bit, uint64_t j) const = 0;
    virtual uint64_t size() const = 0;
    // Aproximación del uso de memoria en RAM (en bytes), incluyendo buffers internos.
    // Nota: no incluye overhead del allocator.
    virtual uint64_t size_in_bytes() const = 0;
};

#endif // BITVECTOR_HPP