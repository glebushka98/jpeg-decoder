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

namespace Jpeg {
static constexpr int CELL_SIZE = 8;
using JpegImageCell = array<array<int, CELL_SIZE>, CELL_SIZE>;
Image Decode(const std::string &filename);
}  // namespace Jpeg
