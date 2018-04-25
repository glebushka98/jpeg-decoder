#include "huffman.h"

HuffmanCode::HuffmanCode() = default;

HuffmanCode::HuffmanCode(vector<uint32_t> counters, vector<uint32_t> values) {
    size_t max_len = 0;
    for (size_t i = 0; i < counters.size(); ++i) {
        if (counters[i] != 0) {
            max_len = i + 1;
        }
    }
    binary_trie_ = vector<pair<bool, bool>>((2 << (max_len + 1)) + 1, {false, false});
    values_ = vector<uint32_t>((2 << (max_len + 1)) + 1, 0);
    uint32_t values_it = 0;
    vector<uint32_t> current_height_vertexes{0};
    for (size_t i = 0; i < max_len; ++i) {
        if (counters[i] > 2 * current_height_vertexes.size()) {
            throw runtime_error("HuffmanCode error!");
        }

        vector<uint32_t> new_vers;
        for (uint32_t j = 0; j < current_height_vertexes.size(); j++) {
            uint32_t son = current_height_vertexes[j] * 2;

            if (j * 2 < counters[i]) {
                binary_trie_[current_height_vertexes[j]].first = true;
                values_[son + 1] = values[values_it++];
            } else {
                if (values_it < values.size()) {
                    new_vers.push_back(son + 1);
                    binary_trie_[current_height_vertexes[j]].first = true;
                }
            }

            if (j * 2 + 1 < counters[i]) {
                binary_trie_[current_height_vertexes[j]].second = true;
                values_[son + 2] = values[values_it++];
            } else {
                if (values_it < values.size()) {
                    new_vers.push_back(son + 2);
                    binary_trie_[current_height_vertexes[j]].second = true;
                }
            }
        }
        swap(new_vers, current_height_vertexes);
    }
}
