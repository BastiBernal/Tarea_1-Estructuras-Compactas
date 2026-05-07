
#include "wavelet/wavelet_tree.hpp"
#include "utils-p/utils.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace {
	template <typename T>
	static std::vector<int64_t> ToSignedValues(const std::vector<T>& values) {
		std::vector<int64_t> result;
		result.reserve(values.size());
		for (const auto& value : values) {
			result.push_back(static_cast<int64_t>(value));
		}
		return result;
	}

	static int64_t ToSignedValue(char value) {
		return static_cast<int64_t>(static_cast<unsigned char>(value));
	}

	static int64_t ToSignedValue(int value) {
		return static_cast<int64_t>(value);
	}
}

WaveletTree::WaveletTree(const std::string& text) {
	build(text);
}

WaveletTree::WaveletTree(const std::vector<int>& values) {
	build(values);
}

std::unique_ptr<WaveletTree::Node> WaveletTree::build_rec(
	const std::vector<int64_t>& data,
	int64_t lo,
	int64_t hi) {
	if (data.empty()) return nullptr;

	auto node = std::make_unique<Node>();
	node->lo = lo;
	node->hi = hi;
	node->n = static_cast<uint64_t>(data.size());

	if (lo == hi) {
		return node;
	}

	const int64_t mid = lo + (hi - lo) / 2;

	std::vector<bool> bits;
	bits.resize(data.size());
	std::vector<int64_t> left_data;
	std::vector<int64_t> right_data;
	left_data.reserve(data.size());
	right_data.reserve(data.size());

	for (std::size_t i = 0; i < data.size(); ++i) {
		const int64_t x = data[i];
		const bool go_right = (x > mid);
		bits[i] = go_right;
		if (go_right) right_data.push_back(x);
		else left_data.push_back(x);
	}

	node->bv = ConstruirBitVectorAuto(bits);

	if (!left_data.empty()) {
		node->left = build_rec(left_data, lo, mid);
	}
	if (!right_data.empty()) {
		const unsigned char rlo = static_cast<unsigned char>(mid + 1);
		node->right = build_rec(right_data, rlo, hi);
	}

	return node;
}

void WaveletTree::build(const std::string& text) {
	root_.reset();
	n_ = static_cast<uint64_t>(text.size());
	if (n_ == 0) return;

	std::vector<int64_t> data;
	data.reserve(text.size());
	int64_t min_c = std::numeric_limits<int64_t>::max();
	int64_t max_c = std::numeric_limits<int64_t>::min();

	for (char ch : text) {
		const int64_t c = ToSignedValue(ch);
		data.push_back(c);
		if (c < min_c) min_c = c;
		if (c > max_c) max_c = c;
	}

	root_ = build_rec(data, min_c, max_c);
}

void WaveletTree::build(const std::vector<int>& values) {
	root_.reset();
	n_ = static_cast<uint64_t>(values.size());
	if (n_ == 0) return;

	const std::vector<int64_t> data = ToSignedValues(values);
	int64_t min_c = std::numeric_limits<int64_t>::max();
	int64_t max_c = std::numeric_limits<int64_t>::min();
	for (int64_t value : data) {
		if (value < min_c) min_c = value;
		if (value > max_c) max_c = value;
	}

	root_ = build_rec(data, min_c, max_c);
}

char WaveletTree::access(uint64_t i) const {
	if (!root_ || i >= n_) return '\0';
	const Node* node = root_.get();
	uint64_t pos = i;
	while (node && node->lo != node->hi) {
		const bool b = node->bv->access(pos);
		if (!b) {
			pos = node->bv->rank(false, pos);
			node = node->left.get();
		} else {
			pos = node->bv->rank(true, pos);
			node = node->right.get();
		}
	}
	return node ? static_cast<char>(node->lo) : '\0';
}

int64_t WaveletTree::access_int(uint64_t i) const {
	if (!root_ || i >= n_) return 0;
	const Node* node = root_.get();
	uint64_t pos = i;
	while (node && node->lo != node->hi) {
		const bool b = node->bv->access(pos);
		if (!b) {
			pos = node->bv->rank(false, pos);
			node = node->left.get();
		} else {
			pos = node->bv->rank(true, pos);
			node = node->right.get();
		}
	}
	return node ? node->lo : 0;
}

uint64_t WaveletTree::rank(char c, uint64_t i) const {
	return rank(static_cast<int>(static_cast<unsigned char>(c)), i);
}

uint64_t WaveletTree::rank(int value, uint64_t i) const {
	if (!root_) return 0;
	if (i > n_) i = n_;

	const int64_t uc = ToSignedValue(value);
	const Node* node = root_.get();
	uint64_t pos = i;

	while (node) {
		if (uc < node->lo || uc > node->hi) return 0;
		if (node->lo == node->hi) return pos;

		const unsigned char mid = static_cast<unsigned char>(node->lo + (node->hi - node->lo) / 2);
		if (uc <= mid) {
			pos = node->bv->rank(false, pos);
			node = node->left.get();
		} else {
			pos = node->bv->rank(true, pos);
			node = node->right.get();
		}
	}
	return 0;
}

uint64_t WaveletTree::select(char c, uint64_t j) const {
	return select(static_cast<int>(static_cast<unsigned char>(c)), j);
}

uint64_t WaveletTree::select(int value, uint64_t j) const {
	if (!root_) return 0;
	if (j == 0) return n_;

	const int64_t uc = ToSignedValue(value);
	const Node* node = root_.get();
	std::vector<const Node*> path;
	std::vector<bool> dirs;
	path.reserve(16);
	dirs.reserve(16);

	while (node) {
		if (uc < node->lo || uc > node->hi) return n_;
		if (node->lo == node->hi) break;

		const unsigned char mid = static_cast<unsigned char>(node->lo + (node->hi - node->lo) / 2);
		const bool go_right = (uc > mid);
		path.push_back(node);
		dirs.push_back(go_right);
		node = go_right ? node->right.get() : node->left.get();
	}

	if (!node) return n_;
	if (j > node->n) return n_;

	uint64_t pos = j; // 1-based posición en el nodo actual
	for (std::size_t idx = path.size(); idx-- > 0;) {
		const Node* cur = path[idx];
		const bool dir = dirs[idx];
		const uint64_t sel = cur->bv->select(dir, pos);
		if (sel == cur->bv->size()) return n_;
		pos = sel + 1;
	}

	const uint64_t ans = pos - 1;
	return (ans < n_) ? ans : n_;
}

uint64_t WaveletTree::size_bytes() const {
	uint64_t total = sizeof(WaveletTree);
	auto node_bytes = [&](const Node* node, const auto& self) -> uint64_t {
		if (!node) return 0;
		uint64_t acc = sizeof(Node);
		if (node->bv) acc += node->bv->bytes();
		acc += self(node->left.get(), self);
		acc += self(node->right.get(), self);
		return acc;
	};
	total += node_bytes(root_.get(), node_bytes);
	return total;
}

// ---------------- WaveletTreeBinary (levelwise / pointerless) ----------------

WaveletTreeBinary::WaveletTreeBinary(const std::string& text) {
	build(text);
}

WaveletTreeBinary::WaveletTreeBinary(const std::vector<int>& values) {
	build(values);
}

int32_t WaveletTreeBinary::build_rec(
	std::vector<std::vector<bool>>& level_bits,
	std::vector<std::vector<NodeInfo>>& level_nodes,
	const std::vector<int64_t>& data,
	int64_t lo,
	int64_t hi,
	uint32_t depth) {
	if (data.empty()) return -1;
	if (lo >= hi) return -1; // hoja (un solo símbolo)

	if (level_bits.size() <= depth) {
		level_bits.resize(static_cast<std::size_t>(depth) + 1);
		level_nodes.resize(static_cast<std::size_t>(depth) + 1);
	}

	NodeInfo info;
	info.lo = lo;
	info.hi = hi;
	info.mid = lo + (hi - lo) / 2;
	info.start = static_cast<uint64_t>(level_bits[static_cast<std::size_t>(depth)].size());
	info.end = info.start + static_cast<uint64_t>(data.size());
	info.left = -1;
	info.right = -1;
	info.rank1_before = 0;
	info.rank0_before = 0;

	const int32_t node_idx = static_cast<int32_t>(level_nodes[static_cast<std::size_t>(depth)].size());
	level_nodes[static_cast<std::size_t>(depth)].push_back(info);

	std::vector<int64_t> left_data;
	std::vector<int64_t> right_data;
	left_data.reserve(data.size());
	right_data.reserve(data.size());

	auto& bits_out = level_bits[static_cast<std::size_t>(depth)];
	bits_out.reserve(bits_out.size() + data.size());

	for (int64_t x : data) {
		const bool go_left = (x <= info.mid);
		// Convención pedida: 1 => izquierdo, 0 => derecho
		bits_out.push_back(go_left);
		if (go_left) left_data.push_back(x);
		else right_data.push_back(x);
	}

	// Hijo izquierdo: rango [lo, mid]
	if (!left_data.empty()) {
		if (lo != info.mid) {
			level_nodes[static_cast<std::size_t>(depth)][static_cast<std::size_t>(node_idx)].left =
				build_rec(level_bits, level_nodes, left_data, lo, info.mid, depth + 1);
		}
	}
	// Hijo derecho: rango [mid+1, hi]
	if (!right_data.empty()) {
		const unsigned char rlo = static_cast<unsigned char>(info.mid + 1);
		if (rlo != hi) {
			level_nodes[static_cast<std::size_t>(depth)][static_cast<std::size_t>(node_idx)].right =
				build_rec(level_bits, level_nodes, right_data, rlo, hi, depth + 1);
		}
	}

	return node_idx;
}

void WaveletTreeBinary::build(const std::string& text) {
	levels_.clear();
	n_ = static_cast<uint64_t>(text.size());
	if (n_ == 0) {
		min_c_ = 0;
		max_c_ = 0;
		return;
	}

	std::vector<int64_t> data;
	data.reserve(text.size());
	int64_t min_c = std::numeric_limits<int64_t>::max();
	int64_t max_c = std::numeric_limits<int64_t>::min();
	for (char ch : text) {
		const int64_t c = ToSignedValue(ch);
		data.push_back(c);
		if (c < min_c) min_c = c;
		if (c > max_c) max_c = c;
	}
	min_c_ = min_c;
	max_c_ = max_c;

	// Caso constante: no hay niveles (no hay bits que guardar)
	if (min_c_ == max_c_) {
		return;
	}

	std::vector<std::vector<bool>> level_bits;
	std::vector<std::vector<NodeInfo>> level_nodes;
	build_rec(level_bits, level_nodes, data, min_c_, max_c_, 0);

	levels_.resize(level_bits.size());
	for (std::size_t d = 0; d < level_bits.size(); ++d) {
		levels_[d].nodes = std::move(level_nodes[d]);
		levels_[d].bv = ConstruirBitVectorAuto(level_bits[d]);
		// Precomputar rank en start para selects rápidos por segmento.
		for (auto& node : levels_[d].nodes) {
			node.rank1_before = levels_[d].bv->rank(true, node.start);
			node.rank0_before = node.start - node.rank1_before;
		}
	}
}

char WaveletTreeBinary::access(uint64_t i) const {
	if (n_ == 0 || i >= n_) return '\0';
	// Caso constante
	if (levels_.empty()) return static_cast<char>(min_c_);

	int32_t node_idx = 0;
	uint64_t pos = i;
	for (std::size_t depth = 0; depth < levels_.size(); ++depth) {
		const auto& lvl = levels_[depth];
		if (node_idx < 0 || static_cast<std::size_t>(node_idx) >= lvl.nodes.size()) return '\0';
		const NodeInfo& node = lvl.nodes[static_cast<std::size_t>(node_idx)];
		const uint64_t len = node.end - node.start;
		if (pos >= len) return '\0';

		const uint64_t global_pos = node.start + pos;
		const bool bit = lvl.bv->access(global_pos);
		const uint64_t ones_before = lvl.bv->rank(true, global_pos) - node.rank1_before;
		if (bit) {
			// ir a hijo izquierdo (1)
			pos = ones_before;
			if (node.lo == node.mid) return static_cast<char>(node.lo);
			node_idx = node.left;
		} else {
			// ir a hijo derecho (0)
			pos = pos - ones_before;
			const unsigned char rlo = static_cast<unsigned char>(node.mid + 1);
			if (rlo == node.hi) return static_cast<char>(node.hi);
			node_idx = node.right;
		}
	}

	return '\0';
}

void WaveletTreeBinary::build(const std::vector<int>& values) {
	levels_.clear();
	n_ = static_cast<uint64_t>(values.size());
	if (n_ == 0) {
		min_c_ = 0;
		max_c_ = 0;
		return;
	}

	const std::vector<int64_t> data = ToSignedValues(values);
	int64_t min_c = std::numeric_limits<int64_t>::max();
	int64_t max_c = std::numeric_limits<int64_t>::min();
	for (int64_t value : data) {
		if (value < min_c) min_c = value;
		if (value > max_c) max_c = value;
	}
	min_c_ = min_c;
	max_c_ = max_c;

	if (min_c_ == max_c_) {
		return;
	}

	std::vector<std::vector<bool>> level_bits;
	std::vector<std::vector<NodeInfo>> level_nodes;
	build_rec(level_bits, level_nodes, data, min_c_, max_c_, 0);

	levels_.resize(level_bits.size());
	for (std::size_t d = 0; d < level_bits.size(); ++d) {
		levels_[d].nodes = std::move(level_nodes[d]);
		levels_[d].bv = ConstruirBitVectorAuto(level_bits[d]);
		for (auto& node : levels_[d].nodes) {
			node.rank1_before = levels_[d].bv->rank(true, node.start);
			node.rank0_before = node.start - node.rank1_before;
		}
	}
}

int64_t WaveletTreeBinary::access_int(uint64_t i) const {
	if (n_ == 0 || i >= n_) return 0;
	if (levels_.empty()) return min_c_;

	int32_t node_idx = 0;
	uint64_t pos = i;
	for (std::size_t depth = 0; depth < levels_.size(); ++depth) {
		const auto& lvl = levels_[depth];
		if (node_idx < 0 || static_cast<std::size_t>(node_idx) >= lvl.nodes.size()) return 0;
		const NodeInfo& node = lvl.nodes[static_cast<std::size_t>(node_idx)];
		const uint64_t len = node.end - node.start;
		if (pos >= len) return 0;

		const uint64_t global_pos = node.start + pos;
		const bool bit = lvl.bv->access(global_pos);
		const uint64_t ones_before = lvl.bv->rank(true, global_pos) - node.rank1_before;
		if (bit) {
			pos = ones_before;
			if (node.lo == node.mid) return node.lo;
			node_idx = node.left;
		} else {
			pos = pos - ones_before;
			if (node.mid + 1 == node.hi) return node.hi;
			node_idx = node.right;
		}
	}

	return 0;
}

uint64_t WaveletTreeBinary::rank(char c, uint64_t i) const {
	return rank(static_cast<int>(static_cast<unsigned char>(c)), i);
}

uint64_t WaveletTreeBinary::rank(int value, uint64_t i) const {
	if (n_ == 0) return 0;
	if (i > n_) i = n_;

	const int64_t uc = ToSignedValue(value);
	// Caso constante
	if (levels_.empty()) {
		return (uc == min_c_) ? i : 0;
	}

	int32_t node_idx = 0;
	uint64_t pos = i;
	for (std::size_t depth = 0; depth < levels_.size(); ++depth) {
		const auto& lvl = levels_[depth];
		if (node_idx < 0 || static_cast<std::size_t>(node_idx) >= lvl.nodes.size()) return 0;
		const NodeInfo& node = lvl.nodes[static_cast<std::size_t>(node_idx)];
		if (uc < node.lo || uc > node.hi) return 0;

		const uint64_t len = node.end - node.start;
		if (pos > len) pos = len;

		const uint64_t ones = lvl.bv->rank(true, node.start + pos) - node.rank1_before;
		if (uc <= node.mid) {
			pos = ones;
			if (node.lo == node.mid) return pos; // hoja
			node_idx = node.left;
		} else {
			pos = pos - ones;
			const unsigned char rlo = static_cast<unsigned char>(node.mid + 1);
			if (rlo == node.hi) return pos; // hoja
			node_idx = node.right;
		}
	}
	return 0;
}

uint64_t WaveletTreeBinary::select(char c, uint64_t j) const {
	return select(static_cast<int>(static_cast<unsigned char>(c)), j);
}

uint64_t WaveletTreeBinary::select(int value, uint64_t j) const {
	if (j == 0) return n_;
	if (n_ == 0) return 0;

	const int64_t uc = ToSignedValue(value);
	// Caso constante
	if (levels_.empty()) {
		if (uc != min_c_) return n_;
		return (j <= n_) ? (j - 1) : n_;
	}

	const uint64_t total = rank(static_cast<int>(uc), n_);
	if (j > total) return n_;

	struct Step {
		uint32_t depth;
		int32_t node_idx;
		bool dir; // 1=izq, 0=der
	};

	std::vector<Step> path;
	path.reserve(levels_.size());

	int32_t node_idx = 0;
	for (std::size_t depth = 0; depth < levels_.size(); ++depth) {
		const auto& lvl = levels_[depth];
		if (node_idx < 0 || static_cast<std::size_t>(node_idx) >= lvl.nodes.size()) return n_;
		const NodeInfo& node = lvl.nodes[static_cast<std::size_t>(node_idx)];
		if (uc < node.lo || uc > node.hi) return n_;

		if (uc <= node.mid) {
			path.push_back(Step{static_cast<uint32_t>(depth), node_idx, true});
			if (node.lo == node.mid) break; // hoja
			node_idx = node.left;
		} else {
			path.push_back(Step{static_cast<uint32_t>(depth), node_idx, false});
			const unsigned char rlo = static_cast<unsigned char>(node.mid + 1);
			if (rlo == node.hi) break; // hoja
			node_idx = node.right;
		}
	}

	uint64_t pos = j; // 1-based dentro de la hoja
	for (std::size_t k = path.size(); k-- > 0;) {
		const Step& step = path[k];
		const auto& lvl = levels_[static_cast<std::size_t>(step.depth)];
		const NodeInfo& node = lvl.nodes[static_cast<std::size_t>(step.node_idx)];
		const uint64_t global_j = (step.dir ? node.rank1_before : node.rank0_before) + pos;
		const uint64_t global_pos = lvl.bv->select(step.dir, global_j);
		if (global_pos >= node.end) return n_;
		pos = (global_pos - node.start) + 1;
	}

	const uint64_t ans = pos - 1;
	return (ans < n_) ? ans : n_;
}

uint64_t WaveletTreeBinary::size_bytes() const {
	uint64_t total = sizeof(WaveletTreeBinary);
	// Storage del vector de niveles (objetos Level), más el contenido de cada nivel.
	total += static_cast<uint64_t>(levels_.capacity()) * sizeof(Level);
	for (const auto& lvl : levels_) {
		if (lvl.bv) total += lvl.bv->bytes();
		total += static_cast<uint64_t>(lvl.nodes.capacity()) * sizeof(NodeInfo);
	}
	return total;
}

