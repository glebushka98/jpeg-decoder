#include "decoder_utils.h"

SpiralTravel::SpiralTravel(int size)
    : size_(size), cur_i(0), cur_j(0), cur_len_(1), line_len_inc_(true), inc_i_(false), to_go_(1) {}

pair<int, int> SpiralTravel::GoNext() {
    int old_i = cur_i, old_j = cur_j;
    if (to_go_ != cur_len_) {
        to_go_++;
        int val = inc_i_ ? 1 : -1;
        cur_i = cur_i + val;
        cur_j = cur_j - val;
        return {old_i, old_j};
    }
    inc_i_ ^= 1;

    if (cur_len_ != size_) {
        cur_len_ += line_len_inc_ ? 1 : -1;
    } else {
        cur_len_ = size_ - 1;
        line_len_inc_ = false;
    }
    if ((line_len_inc_ && cur_i == 0) || (!line_len_inc_ && cur_j == size_ - 1)) {
        if (line_len_inc_) {
            cur_j++;
        } else {
            cur_i++;
        }
    } else {
        if (line_len_inc_) {
            cur_i++;
        } else {
            cur_j++;
        }
    }
    to_go_ = 1;

    return {old_i, old_j};
}

bool SpiralTravel::HasNext() { return (cur_i != size_ && cur_j != size_); }

BitsReader::BitsReader() = default;

BitsReader::BitsReader(vector<char> &bytes) { bytes_.assign(bytes.begin(), bytes.end()); }

uint32_t BitsReader::GetBits(size_t count, bool peek) {
    if (count > 8) {
        return (GetBits(8, peek) << (count - 8)) + GetBits(count - 8, peek);
    }
    if (count + bits_pointer > BYTE_SIZE) {
        if (bytes_pointer_ + 1 >= bytes_.size()) {
            throw runtime_error("Error while reading next byte");
        }
    }

    if (count + bits_pointer <= BYTE_SIZE) {
        uint32_t res = GetUIntLR(bytes_[bytes_pointer_], bits_pointer, count + bits_pointer - 1);
        if (!peek) {
            bits_pointer += count;
            if (bits_pointer == 8) {
                bits_pointer = 0;
                bytes_pointer_ += 1;
            }
        }
        return res;
    }

    auto left = GetUIntLR(bytes_[bytes_pointer_], bits_pointer, 7);
    auto right = GetUIntLR(bytes_[bytes_pointer_ + 1], 0, count - (8 - bits_pointer) - 1);
    if (!peek) {
        bytes_pointer_ += 1;
        bits_pointer = count - (8 - bits_pointer);
        return (left << bits_pointer) + right;
    } else {
        return (left << (count - (8 - bits_pointer))) + right;
    }
}

uint32_t BitsReader::GetByte() {
    if (bits_pointer == 0) {
        if (bytes_pointer_ >= bytes_.size()) {
            throw runtime_error("Error while reading next byte");
        }
        return bytes_[bytes_pointer_++];
    } else {
        return GetBits(8);
    }
}

size_t BitsReader::GetBitsPointer() { return bits_pointer; }

void BitsReader::ConsumeBits(size_t count) {
    bits_pointer += count;
    if (bits_pointer >= 8) {
        bytes_pointer_ += bits_pointer / 8;
        bits_pointer %= 8;
        if ((bytes_pointer_ >= bytes_.size() && bits_pointer != 0) ||
            bytes_pointer_ >= bytes_.size() + 1) {
            throw runtime_error("Error while reading next byte!");
        }
    }
}

void BitsReader::SquashFF00() {
    if (bits_pointer != 0) {
        throw runtime_error("Can't squash with bits_pointer != 0");
    }

    auto cur_pos = bytes_pointer_;
    auto cur = bytes_pointer_;
    while (cur < bytes_.size()) {
        if (cur + 1 < bytes_.size() && static_cast<uint32_t>(bytes_[cur]) == 255 &&
            static_cast<uint32_t>(bytes_[cur + 1]) == 0) {
            bytes_[cur_pos] = bytes_[cur];
            cur_pos++;
            cur += 2;
        } else {
            bytes_[cur_pos] = bytes_[cur];
            cur_pos++;
            cur++;
        }
    }
    bytes_.resize(cur_pos);
}

uint32_t BitsReader::GetBit() { return GetBits(1); }

bool BitsReader::CheckTwoBytes(uint8_t byte_first, uint8_t byte_second) {
    uint32_t byte1 = GetByte();
    uint32_t byte2 = GetByte();
    bytes_pointer_ -= 2;
    return byte1 == byte_first && byte2 == byte_second;
}

uint32_t BitsReader::GetInt(size_t bytes_count) {
    uint32_t res = 0;
    for (size_t i = 0; i < bytes_count; i++) {
        res <<= 8;
        res += GetByte();
    }
    return res;
}

bool BitsReader::Empty() { return bytes_pointer_ >= bytes_.size(); }

uint8_t BitsReader::GetUIntLR(uint8_t num, size_t l, size_t r) {
    return static_cast<uint8_t>(num << l) >> (BYTE_SIZE - (r - l + 1));
}
