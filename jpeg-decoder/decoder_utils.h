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
#include "image.h"

using std::array;
using std::pair;
using std::runtime_error;
using std::vector;

class BitsReader {
public:
    BitsReader();

    explicit BitsReader(vector<char>& bytes);

    uint32_t GetBits(size_t count, bool peek = false);

    uint32_t GetByte();

    size_t GetBitsPointer();

    void ConsumeBits(size_t count);

    void SquashFF00();

    uint32_t GetBit();

    bool CheckTwoBytes(uint8_t byte_first, uint8_t byte_second);

    uint32_t GetInt(size_t bytes_count);

    bool Empty();

private:
    uint8_t GetUIntLR(uint8_t num, size_t l, size_t r);

    static constexpr int BYTE_SIZE = 8;

    size_t bits_pointer = 0;
    size_t bytes_pointer_ = 0;

    vector<uint8_t> bytes_;
};

class SpiralTravel {
public:
    SpiralTravel(int size);

    pair<int, int> GoNext();

    bool HasNext();

private:
    int size_;
    int cur_i;
    int cur_j;
    int cur_len_;
    bool line_len_inc_;
    bool inc_i_;
    int to_go_;
};