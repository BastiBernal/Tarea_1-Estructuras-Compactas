#include <cassert>
#include <cstdint>
#include <string>
#include <vector>

#include "wavelet/wavelet_tree.hpp"

namespace {
	template <typename T>
	uint64_t naive_rank(const std::vector<T>& values, T value, uint64_t i) {
		if (i > values.size()) i = values.size();
		uint64_t count = 0;
		for (uint64_t idx = 0; idx < i; ++idx) {
			if (values[static_cast<std::size_t>(idx)] == value) {
				++count;
			}
		}
		return count;
	}

	template <typename T>
	uint64_t naive_select(const std::vector<T>& values, T value, uint64_t j) {
		if (j == 0) return values.size();
		uint64_t count = 0;
		for (uint64_t idx = 0; idx < values.size(); ++idx) {
			if (values[static_cast<std::size_t>(idx)] == value) {
				++count;
				if (count == j) {
					return idx;
				}
			}
		}
		return values.size();
	}

	void test_wavelet_from_string() {
		const std::string text = "banana_bandana";
		WaveletTree tree(text);

		assert(tree.size() == text.size());
		assert(!tree.empty());
		assert(tree.access(0) == 'b');
		assert(tree.access(1) == 'a');
		assert(tree.access(5) == 'a');
		assert(tree.access(text.size() - 1) == 'a');
		assert(tree.access(text.size()) == '\0');

		assert(tree.rank('a', 0) == 0);
		assert(tree.rank('a', 1) == 0);
		assert(tree.rank('a', 2) == 1);
		assert(tree.rank('a', text.size()) == 6);
		assert(tree.rank('n', text.size()) == 4);
		assert(tree.rank('z', text.size()) == 0);

		assert(tree.select('a', 1) == 1);
		assert(tree.select('a', 2) == 3);
		assert(tree.select('n', 1) == 2);
		assert(tree.select('n', 4) == 12);
		assert(tree.select('z', 1) == text.size());
		assert(tree.select('a', 0) == text.size());

		WaveletTreeBinary binary(text);
		assert(binary.size() == text.size());
		assert(binary.access(0) == 'b');
		assert(binary.access(1) == 'a');
		assert(binary.access(text.size()) == '\0');
		assert(binary.rank('a', text.size()) == 6);
		assert(binary.select('n', 4) == 12);
	}

	void test_wavelet_from_int_vector() {
		const std::vector<int> values = {3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5};
		WaveletTree tree(values);
		WaveletTreeBinary binary(values);

		assert(tree.size() == values.size());
		assert(binary.size() == values.size());
		assert(tree.empty() == false);
		assert(binary.empty() == false);

		for (uint64_t i = 0; i < values.size(); ++i) {
			assert(tree.access_int(i) == values[static_cast<std::size_t>(i)]);
			assert(binary.access_int(i) == values[static_cast<std::size_t>(i)]);
		}
		assert(tree.access_int(values.size()) == 0);
		assert(binary.access_int(values.size()) == 0);

		for (int value : {1, 2, 3, 4, 5, 6, 9, 7}) {
			for (uint64_t i = 0; i <= values.size(); ++i) {
				assert(tree.rank(value, i) == naive_rank(values, value, i));
				assert(binary.rank(value, i) == naive_rank(values, value, i));
			}
			const uint64_t total = naive_rank(values, value, values.size());
			for (uint64_t j = 0; j <= total + 1; ++j) {
				assert(tree.select(value, j) == (j == 0 ? values.size() : naive_select(values, value, j)));
				assert(binary.select(value, j) == (j == 0 ? values.size() : naive_select(values, value, j)));
			}
		}
	}

} // namespace

int main() {
	test_wavelet_from_string();
	test_wavelet_from_int_vector();
	return 0;
}