#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using std::array;
using std::pair;
using std::runtime_error;
using std::vector;

class HuffmanCode {
public:
    HuffmanCode();

    HuffmanCode(vector<uint32_t> counters, vector<uint32_t> values);

    template <typename Stream>
    uint32_t FindCode(Stream&& stream) const {
        size_t cur_ver = 0;

        while (binary_trie_.size() > cur_ver &&
               (binary_trie_[cur_ver].first || binary_trie_[cur_ver].second)) {
            if (cur_ver >= binary_trie_.size()) {
                return false;
            }
            int cur_bit = stream.GetBit();
            if (cur_bit == 0 && binary_trie_[cur_ver].first) {
                cur_ver = cur_ver * 2 + 1;
                continue;
            }

            if (cur_bit == 1 && binary_trie_[cur_ver].second) {
                cur_ver = cur_ver * 2 + 2;
                continue;
            }
            throw runtime_error("wrong code");
        }
        if (cur_ver >= binary_trie_.size()) {
            throw runtime_error("wrong code");
        }
        if (binary_trie_[cur_ver].first || binary_trie_[cur_ver].second) {
            throw runtime_error("wrong code");
        }
        return values_[cur_ver];
    }

private:
    vector<pair<bool, bool>> binary_trie_;
    vector<uint32_t> values_;
};