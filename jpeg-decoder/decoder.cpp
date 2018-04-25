#include "decoder.h"
#include "decoder_utils.h"
#include "huffman.h"

#include <fftw3.h>
#include <array>
#include <cassert>
#include <cctype>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using std::array;
using std::cerr;
using std::cout;
using std::endl;
using std::istream;
using std::max;
using std::min;
using std::pair;
using std::runtime_error;
using std::set;
using std::vector;

namespace Jpeg {

class JpegMarker {
public:
    JpegMarker() = default;

    JpegMarker(uint8_t b_first, uint8_t b_second) : first_byte_(b_first), second_byte_(b_second) {}

    uint8_t GetFirstByte() const { return first_byte_; }

    uint8_t GetSecondByte() const { return second_byte_; }

    bool operator==(const JpegMarker &jpeg_marker) const {
        return jpeg_marker.first_byte_ == first_byte_ && jpeg_marker.second_byte_ == second_byte_;
    }

    bool operator<(const JpegMarker &jpeg_marker) const {
        if (first_byte_ < jpeg_marker.first_byte_) {
            return true;
        }

        return (second_byte_ < jpeg_marker.second_byte_);
    }

    const static JpegMarker SOI;
    const static JpegMarker EOI;
    const static JpegMarker DQT;
    const static JpegMarker SOF0;
    const static JpegMarker DHT;
    const static JpegMarker COM;
    const static JpegMarker SOS;

private:
    uint8_t first_byte_;
    uint8_t second_byte_;
};

const JpegMarker JpegMarker::SOI('\xff', '\xd8');
const JpegMarker JpegMarker::EOI('\xff', '\xd9');
const JpegMarker JpegMarker::DQT('\xff', '\xdb');
const JpegMarker JpegMarker::SOF0('\xff', '\xc0');
const JpegMarker JpegMarker::DHT('\xff', '\xc4');
const JpegMarker JpegMarker::COM('\xff', '\xfe');
const JpegMarker JpegMarker::SOS('\xff', '\xda');

class JpegImageInfoHolder {
public:
    void AddProcessedSection(const JpegMarker &jm) { processed_sections_.insert(jm); }

    bool WasSectionProcessed(const JpegMarker &jm) { return processed_sections_.count(jm) > 0; }

    void SetImageSize(size_t width, size_t height) {
        height_ = height;
        width_ = width;
        image_.SetSize(width, height);
    }

    void SetCompsCount(size_t comps) { comps_count_ = comps; }

    size_t GetCompsCount() const { return comps_count_; }

    void SetMcuSize(size_t mcu_w, size_t mcu_h) {
        mcu_w_ = mcu_w;
        mcu_h_ = mcu_h;
    }

    void SetYQuantizationTableId(size_t y_qid) { y_qid_ = y_qid; }

    size_t GetYQuantizationTableId() const { return y_qid_; }

    void SetCbQuantizationTableId(size_t cb_qid) { cb_qid_ = cb_qid; }

    size_t GetCbQuantizationTableId() const { return cb_qid_; }

    void SetCrQuantizationTableId(size_t cr_qid) { cr_qid_ = cr_qid; }

    size_t GetCrQuantizationTableId() const { return cr_qid_; }

    void SetImageComment(const std::string &s) { image_.SetComment(s); }

    void AddHuffmanTable(size_t id, const HuffmanCode &code) {
        assert(huffmans_.count(id) == 0);
        huffmans_[id] = code;
    }

    HuffmanCode &GetHuffmanTable(size_t id) {
        assert(huffmans_.count(id) > 0);
        return huffmans_[id];
    }

    bool HasHuffmanTable(size_t id) { return huffmans_.count(id) > 0; }

    void SetQuantizationTable(size_t id, const JpegImageCell &cell) {
        assert(dqts_.count(id) == 0);
        dqts_[id] = cell;
    }

    JpegImageCell &GetQuantizationTable(size_t id) { return dqts_[id]; }

    Image GetImage() { return image_; }

    void SetImagePixel(int y, int x, const RGB &rgb) { image_.SetPixel(y, x, rgb); }

    size_t GetImageHeight() const { return height_; }

    size_t GetImageWidth() const { return width_; }

    size_t GetMcuHeight() const { return mcu_h_; }

    size_t GetMcuWidth() const { return mcu_w_; }

private:
    set<JpegMarker> processed_sections_;
    std::unordered_map<size_t, HuffmanCode> huffmans_;
    std::unordered_map<size_t, JpegImageCell> dqts_;

    size_t height_ = 0;
    size_t width_ = 0;

    size_t comps_count_ = 0;
    size_t mcu_w_ = 0;
    size_t mcu_h_ = 0;
    size_t y_qid_ = 0;
    size_t cb_qid_ = 0;
    size_t cr_qid_ = 0;

    Image image_;
};

class JpegDecoder {
public:
    explicit JpegDecoder(std::istream &in) : in_(in) {}

    void Initialize() {
        in_.seekg(0, in_.end);
        auto length_ = static_cast<uint32_t>(in_.tellg());
        in_.seekg(0, in_.beg);
        vector<char> bytes_(max(length_, 1u));
        in_.read(bytes_.data(), length_);
        if (!in_) {
            throw std::runtime_error("error: only " + std::to_string(in_.gcount()) +
                                     "bytes could be read!");
        }

        bits_reader_ = BitsReader(bytes_);
    }

    void Validate(const std::string &error_descr, bool is_ok) {
        if (!is_ok) {
            throw runtime_error("Error = " + error_descr);
        }
    }

    bool ProcessSOI() {
        const JpegMarker &marker = JpegMarker::SOI;
        if (!bits_reader_.CheckTwoBytes(marker.GetFirstByte(), marker.GetSecondByte())) {
            return false;
        }

        Validate("Jpeg must have only one SOI", !info_holder_.WasSectionProcessed(marker));
        bits_reader_.ConsumeBits(16);
        info_holder_.AddProcessedSection(marker);
        return true;
    }

    bool ProcessEOI() {
        const JpegMarker &marker = JpegMarker::EOI;
        if (!bits_reader_.CheckTwoBytes(marker.GetFirstByte(), marker.GetSecondByte())) {
            return false;
        }
        Validate("Jpeg must have only one EOI", !info_holder_.WasSectionProcessed(marker));
        bits_reader_.ConsumeBits(16);

        info_holder_.AddProcessedSection(marker);
        return true;
    }

    bool ProcessDQT() {
        const JpegMarker &marker = JpegMarker::DQT;
        if (!bits_reader_.CheckTwoBytes(marker.GetFirstByte(), marker.GetSecondByte())) {
            return false;
        }
        bits_reader_.ConsumeBits(16);

        size_t length = bits_reader_.GetInt(2);
        Validate("length >= 2", length >= 2);
        length -= 2;
        while (length) {
            size_t dig_size = bits_reader_.GetBits(4) + 1;
            size_t table_id = bits_reader_.GetBits(4);

            Validate("Invalid DQT section size", dig_size * 64 + 1 <= length);

            SpiralTravel sp(8);

            JpegImageCell jic{};

            for (int i = 0; i < 64; i++) {
                auto indexes = sp.GoNext();

                jic[indexes.first][indexes.second] = bits_reader_.GetInt(dig_size);
            }
            info_holder_.SetQuantizationTable(table_id, jic);
            length -= 64 * dig_size + 1;
        }

        info_holder_.AddProcessedSection(marker);
        return true;
    }

    bool ProcessSOF0() {
        const JpegMarker &marker = JpegMarker::SOF0;
        if (!bits_reader_.CheckTwoBytes(marker.GetFirstByte(), marker.GetSecondByte())) {
            return false;
        }
        bits_reader_.ConsumeBits(16);

        Validate("Jpeg must have only one SOF0", !info_holder_.WasSectionProcessed(marker));

        uint32_t len = bits_reader_.GetInt(2);
        Validate("SOF0 invalid size_of byte", bits_reader_.GetByte() == 8);
        size_t height_ = bits_reader_.GetInt(2);
        size_t width_ = bits_reader_.GetInt(2);
        size_t comps_count_ = bits_reader_.GetInt(1);

        info_holder_.SetImageSize(width_, height_);
        info_holder_.SetCompsCount(comps_count_);

        Validate("Invalid SOF0 components count", comps_count_ == 1 || comps_count_ == 3);

        Validate("Invalid SOF0 section length", len == 8 + 3 * comps_count_);

        // Y component
        uint32_t y_id = bits_reader_.GetInt(1);
        uint32_t y_h = bits_reader_.GetBits(4);
        uint32_t y_v = bits_reader_.GetBits(4);
        info_holder_.SetYQuantizationTableId(bits_reader_.GetInt(1));

        if (comps_count_ == 3) {
            // Cb component
            uint32_t cb_id = bits_reader_.GetInt(1);
            uint32_t cb_h = bits_reader_.GetBits(4);
            uint32_t cb_v = bits_reader_.GetBits(4);
            info_holder_.SetCbQuantizationTableId(bits_reader_.GetInt(1));

            // Cr component
            uint32_t cr_id = bits_reader_.GetInt(1);
            uint32_t cr_h = bits_reader_.GetBits(4);
            uint32_t cr_v = bits_reader_.GetBits(4);
            info_holder_.SetCrQuantizationTableId(bits_reader_.GetInt(1));

            auto h_max = max<uint32_t>({y_h, cb_h, cr_h});
            auto v_max = max<uint32_t>({y_v, cb_v, cr_v});

            Validate("Invalid IDS order", y_id == 1 && cb_id == 2 && cr_id == 3);
            Validate("No downsampling for y", y_h == h_max && y_v == v_max);
            Validate("Invalid equality of downsampling", cb_h == cr_h && cr_v == cb_v);

            info_holder_.SetMcuSize(h_max / cb_h, v_max / cb_v);
        } else {
            info_holder_.SetMcuSize(1, 1);
        }

        info_holder_.AddProcessedSection(marker);
        return true;
    }

    bool ProcessDHT() {
        const JpegMarker &marker = JpegMarker::DHT;
        if (!bits_reader_.CheckTwoBytes(marker.GetFirstByte(), marker.GetSecondByte())) {
            return false;
        }

        bits_reader_.ConsumeBits(16);

        uint32_t len = bits_reader_.GetInt(2);
        len -= 2;
        while (len) {
            uint32_t id = bits_reader_.GetInt(1);
            len -= 1;
            vector<uint32_t> counters(16, 0);
            vector<uint32_t> values;
            size_t cnt = 0;
            Validate("Invalid DHT section len", len >= 16);
            len -= 16;

            for (size_t i = 0; i < 16; i++) {
                counters[i] = bits_reader_.GetInt(1);
                cnt += counters[i];
            }
            Validate("Invalid DHT section len", len >= cnt);

            for (size_t i = 0; i < cnt; i++) {
                values.push_back(bits_reader_.GetInt(1));
            }

            info_holder_.AddHuffmanTable(id, HuffmanCode(counters, values));
            len -= cnt;
        }

        info_holder_.AddProcessedSection(marker);

        return true;
    }

    bool ProcessAPPn() {
        for (int i = 0; i < 15; i++) {
            if (bits_reader_.CheckTwoBytes('\xff', static_cast<uint8_t>(14 * 16 + i))) {
                bits_reader_.ConsumeBits(16);
                size_t len = bits_reader_.GetInt(2);
                bits_reader_.ConsumeBits(8 * (len - 2));
                return true;
            }
        }
        return false;
    }

    bool ProcessCOM() {
        const JpegMarker &marker = JpegMarker::COM;
        if (!bits_reader_.CheckTwoBytes(marker.GetFirstByte(), marker.GetSecondByte())) {
            return false;
        }
        bits_reader_.ConsumeBits(16);
        uint32_t len = bits_reader_.GetInt(2);
        std::string comment;
        for (size_t i = 0; i < len - 2; i++) {
            comment.push_back(static_cast<char>(bits_reader_.GetByte()));
        }
        info_holder_.SetImageComment(comment);
        info_holder_.AddProcessedSection(marker);
        return true;
    }

    void ReadImageCell(JpegImageCell &res, const HuffmanCode &hf_dc, const HuffmanCode &hf_ac) {
        auto set_value = [&](int i, int j, size_t len) mutable {
            uint32_t val = bits_reader_.GetBits(len);
            if (val & (1u << (len - 1))) {
                res[i][j] = val;
            } else {
                res[i][j] = static_cast<int>(val) - (1 << len) + 1;
            }
        };

        {
            uint32_t dc = hf_dc.FindCode(bits_reader_);

            if (dc == 0) {
                res[0][0] = 0;
            } else {
                set_value(0, 0, dc);
            }
        }

        SpiralTravel st(8);
        st.GoNext();

        bool zeros = false;

        while (st.HasNext()) {
            if (zeros) {
                auto index = st.GoNext();
                res[index.first][index.second] = 0;
                continue;
            }
            uint32_t ac = hf_ac.FindCode(bits_reader_);
            if (ac == 0) {
                auto index = st.GoNext();
                zeros = true;
                res[index.first][index.second] = 0;
                continue;
            }
            int zeros_count = ac >> 4;
            size_t len = ac & (0xF);
            while (zeros_count && st.HasNext()) {
                auto index = st.GoNext();
                res[index.first][index.second] = 0;
                --zeros_count;
            }

            Validate("While reading matrix in sos", !zeros_count);
            if (len) {
                Validate("While reading matrix in sos", st.HasNext());

                auto index = st.GoNext();
                set_value(index.first, index.second, len);
            } else {
                Validate("While reading matrix in sos", st.HasNext());

                auto index = st.GoNext();
                res[index.first][index.second] = 0;
            }
        }
    }

    static void DotCells(JpegImageCell &lhs, JpegImageCell &rhs) {
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 8; j++) {
                lhs[i][j] *= rhs[i][j];
            }
        }
    }

    void TransformWithFastFft(JpegImageCell &ar) {
        using fftw_plan_t = std::unique_ptr<fftw_plan_s, decltype(&fftw_destroy_plan)>;
        using fftw_double_t = std::unique_ptr<double, decltype(&fftw_free)>;

        auto buffer_size = CELL_SIZE * CELL_SIZE * sizeof(double);

        fftw_double_t m_in = {(double *)fftw_malloc(buffer_size), &fftw_free};
        fftw_double_t m_out = {(double *)fftw_malloc(buffer_size), &fftw_free};

        for (size_t i = 0; i < CELL_SIZE; i++) {
            for (size_t j = 0; j < CELL_SIZE; j++) {
                m_in.get()[i * CELL_SIZE + j] = ar[i][j];

                if (i == 0) {
                    m_in.get()[i * CELL_SIZE + j] *= sqrt(2);
                }

                if (j == 0) {
                    m_in.get()[i * CELL_SIZE + j] *= sqrt(2);
                }
            }
        }

        fftw_plan_t m_plan = {(fftw_plan_r2r_2d(CELL_SIZE, CELL_SIZE, m_in.get(), m_out.get(),
                                                FFTW_REDFT01, FFTW_REDFT01, 0)),
                              &fftw_destroy_plan};

        fftw_execute(m_plan.get());

        for (size_t i = 0; i < CELL_SIZE; i++) {
            for (size_t j = 0; j < CELL_SIZE; j++) {
                ar[i][j] = static_cast<int>(m_out.get()[i * CELL_SIZE + j] / (2 * CELL_SIZE));
            }
        }
    }

    static RGB YCbCrToRGB(int y, int cb, int cr) {
        RGB rgb{};

        auto cut_bounds = [&](int val) {
            if (val > 255) {
                return 255;
            }
            if (val < 0) {
                return 0;
            }
            return val;
        };

        rgb.r = static_cast<int>(y + 1.402 * cr + 128);
        rgb.g = static_cast<int>(y - 0.34414 * cb - 0.71414 * cr + 128);
        rgb.b = static_cast<int>(y + 1.772 * cb + 128);

        rgb.r = cut_bounds(rgb.r);
        rgb.g = cut_bounds(rgb.g);
        rgb.b = cut_bounds(rgb.b);
        return rgb;
    }

    void TransofrmSosWithFft(size_t ds_h, size_t ds_w, size_t y_qti, size_t cb_qti, size_t cr_qti,
                             std::vector<std::vector<JpegImageCell>> &y_mats, JpegImageCell &cb_mat,
                             JpegImageCell &cr_mat) {
        for (size_t ini = 0; ini < ds_h; ini++) {
            for (size_t inj = 0; inj < ds_w; inj++) {
                DotCells(y_mats[ini][inj], info_holder_.GetQuantizationTable(y_qti));
                TransformWithFastFft(y_mats[ini][inj]);
            }
        }
        {
            DotCells(cb_mat, info_holder_.GetQuantizationTable(cb_qti));
            TransformWithFastFft(cb_mat);
        }
        {
            DotCells(cr_mat, info_holder_.GetQuantizationTable(cr_qti));
            TransformWithFastFft(cr_mat);
        }
    }

    void SetFinalPixelsInSos(size_t ds_h, size_t ds_w,
                             std::vector<std::vector<JpegImageCell>> &y_mats, JpegImageCell &cb_mat,
                             JpegImageCell &cr_mat, size_t actual_height, size_t actual_width,
                             size_t big_block_i, size_t big_block_j) {
        for (size_t ini = 0; ini < ds_h; ini++) {
            for (size_t inj = 0; inj < ds_w; inj++) {
                for (int block_i = 0; block_i < 8; block_i++) {
                    for (int block_j = 0; block_j < 8; block_j++) {
                        auto cur_i = ini * 8 + block_i;
                        auto cur_j = inj * 8 + block_j;

                        auto c_i = cur_i / ds_h;
                        auto c_j = cur_j / ds_w;

                        RGB rgb = YCbCrToRGB(y_mats[ini][inj][block_i][block_j], cb_mat[c_i][c_j],
                                             cr_mat[c_i][c_j]);
                        auto si = big_block_i * 8 * ds_h + cur_i;
                        auto sj = big_block_j * 8 * ds_w + cur_j;
                        if (si < actual_height && sj < actual_width) {
                            info_holder_.SetImagePixel(si, sj, rgb);
                        }
                    }
                }
            }
        }
    }

    void ReadSosMcu(uint32_t comps, uint32_t y_dc, uint32_t y_ac, uint32_t cb_dc, uint32_t cb_ac,
                    uint32_t cr_dc, uint32_t cr_ac, int &last_dc_y, int &last_dc_cb,
                    int &last_dc_cr, size_t ds_h_, size_t ds_w_,
                    vector<vector<JpegImageCell>> &y_mats, JpegImageCell &cb_mat,
                    JpegImageCell &cr_mat) {
        {
            {
                // Y reading
                size_t cur_y_i = 0;
                size_t cur_y_j = 0;

                Validate("Invalid huffman table id",
                         info_holder_.HasHuffmanTable((0 << 4) + y_dc) &&
                             info_holder_.HasHuffmanTable((1 << 4) + y_ac));
                auto &hf_dc = info_holder_.GetHuffmanTable((0 << 4) + y_dc);
                auto &hf_ac = info_holder_.GetHuffmanTable((1 << 4) + y_ac);

                for (size_t ds_size = 0; ds_size < ds_h_ * ds_w_; ++ds_size) {
                    ReadImageCell(y_mats[cur_y_i][cur_y_j], hf_dc, hf_ac);

                    y_mats[cur_y_i][cur_y_j][0][0] = last_dc_y + y_mats[cur_y_i][cur_y_j][0][0];
                    last_dc_y = y_mats[cur_y_i][cur_y_j][0][0];

                    if (cur_y_i == 0 && cur_y_j == ds_w_ - 1) {
                        cur_y_i = 1;
                        cur_y_j = 0;
                    } else {
                        cur_y_j++;
                    }
                }
            }
        }
        if (comps >= 3) {
            // Cb reading

            Validate("Invalid huffman table id",
                     info_holder_.HasHuffmanTable((0 << 4) + cb_dc) &&
                         info_holder_.HasHuffmanTable((1 << 4) + cb_ac));
            auto &hf_dc = info_holder_.GetHuffmanTable((0 << 4) + cb_dc);
            auto &hf_ac = info_holder_.GetHuffmanTable((1 << 4) + cb_ac);
            ReadImageCell(cb_mat, hf_dc, hf_ac);
            cb_mat[0][0] = last_dc_cb + cb_mat[0][0];
            last_dc_cb = cb_mat[0][0];
        }
        if (comps >= 3) {
            // Cr reading
            Validate("Invalid huffman table id",
                     info_holder_.HasHuffmanTable((0 << 4) + cr_dc) &&
                         info_holder_.HasHuffmanTable((1 << 4) + cr_ac));

            auto &hf_dc = info_holder_.GetHuffmanTable((0 << 4) + cr_dc);
            auto &hf_ac = info_holder_.GetHuffmanTable((1 << 4) + cr_ac);
            ReadImageCell(cr_mat, hf_dc, hf_ac);
            cr_mat[0][0] = last_dc_cr + cr_mat[0][0];
            last_dc_cr = cr_mat[0][0];
        }
    }
    bool ProcessSOS() {
        const JpegMarker &marker = JpegMarker::SOS;
        if (!bits_reader_.CheckTwoBytes(marker.GetFirstByte(), marker.GetSecondByte())) {
            return false;
        }

        Validate("Jpeg must have only one SOS", !info_holder_.WasSectionProcessed(marker));

        {
            Validate("Was SOI section", info_holder_.WasSectionProcessed(JpegMarker::SOI));
            Validate("Was SOF0 section", info_holder_.WasSectionProcessed(JpegMarker::SOF0));
        }

        auto actual_width = info_holder_.GetImageWidth();
        auto actual_height = info_holder_.GetImageHeight();

        auto block_h = (8 * info_holder_.GetMcuHeight());
        auto block_w = (8 * info_holder_.GetMcuWidth());

        auto width_ = ((actual_width + block_w - 1) / block_w) * block_w;
        auto height_ = ((actual_height + block_h - 1) / block_h) * block_h;

        bits_reader_.ConsumeBits(16);

        auto head_len = bits_reader_.GetInt(2);
        auto comps = bits_reader_.GetInt(1);

        Validate("Equals comps", comps == info_holder_.GetCompsCount());
        Validate("Head len", head_len == 6 + comps * 2);
        auto y_id = bits_reader_.GetInt(1);
        Validate("Y id", y_id == 1);
        auto y_dc = bits_reader_.GetBits(4);
        auto y_ac = bits_reader_.GetBits(4);

        uint32_t cb_id = 0;
        uint32_t cb_dc = 0;
        uint32_t cb_ac = 0;
        uint32_t cr_id = 0;
        uint32_t cr_dc = 0;
        uint32_t cr_ac = 0;

        if (comps == 3) {
            cb_id = bits_reader_.GetInt(1);
            Validate("Cb id", cb_id == 2);
            cb_dc = bits_reader_.GetBits(4);
            cb_ac = bits_reader_.GetBits(4);

            cr_id = bits_reader_.GetInt(1);
            Validate("Cr id", cr_id == 3);
            cr_dc = bits_reader_.GetBits(4);
            cr_ac = bits_reader_.GetBits(4);
        }

        Validate("Poebota", bits_reader_.GetInt(1) == 0 && bits_reader_.GetInt(1) == '\x3f' &&
                                bits_reader_.GetInt(1) == 0);

        Validate("% 8 *", height_ % block_h == 0 && width_ % block_w == 0);

        std::vector<std::vector<JpegImageCell>> y_mats(
            info_holder_.GetMcuHeight(), vector<JpegImageCell>(info_holder_.GetMcuWidth()));
        JpegImageCell cb_mat{};
        JpegImageCell cr_mat{};

        int last_dc_y = 0;
        int last_dc_cb = 0;
        int last_dc_cr = 0;
        bits_reader_.SquashFF00();

        auto ds_h_ = info_holder_.GetMcuHeight();
        auto ds_w_ = info_holder_.GetMcuWidth();

        auto y_qti = info_holder_.GetYQuantizationTableId();
        auto cb_qti = info_holder_.GetCbQuantizationTableId();
        auto cr_qti = info_holder_.GetCrQuantizationTableId();

        for (size_t big_block_i = 0; big_block_i < height_ / block_h; ++big_block_i) {
            for (size_t big_block_j = 0; big_block_j < width_ / block_w; ++big_block_j) {
                ReadSosMcu(comps, y_dc, y_ac, cb_dc, cb_ac, cr_dc, cr_ac, last_dc_y, last_dc_cb,
                           last_dc_cr, ds_h_, ds_w_, y_mats, cb_mat, cr_mat);

                TransofrmSosWithFft(ds_h_, ds_w_, y_qti, cb_qti, cr_qti, y_mats, cb_mat, cr_mat);

                SetFinalPixelsInSos(ds_h_, ds_w_, y_mats, cb_mat, cr_mat, actual_height,
                                    actual_width, big_block_i, big_block_j);
            }
        }

        if (bits_reader_.GetBitsPointer() != 0) {
            bits_reader_.ConsumeBits(8 - bits_reader_.GetBitsPointer());
        }

        info_holder_.AddProcessedSection(marker);
        return true;
    }

    Image Decode() {
        vector<bool (JpegDecoder::*)(void)> section_processors{
            &JpegDecoder::ProcessDQT, &JpegDecoder::ProcessAPPn, &JpegDecoder::ProcessCOM,
            &JpegDecoder::ProcessDHT, &JpegDecoder::ProcessEOI,  &JpegDecoder::ProcessSOF0,
            &JpegDecoder::ProcessSOI, &JpegDecoder::ProcessSOS};

        while (!bits_reader_.Empty()) {
            bool was_one_section = false;
            for (auto processor : section_processors) {
                if ((this->*processor)()) {
                    was_one_section = true;
                    break;
                }
            }
            if (info_holder_.WasSectionProcessed(JpegMarker::EOI)) {
                break;
            }
            Validate("There was no marker here!", was_one_section);
        }

        Validate("There was no SOS marker!", info_holder_.WasSectionProcessed(JpegMarker::SOS));

        return info_holder_.GetImage();
    }

    JpegImageInfoHolder info_holder_;
    BitsReader bits_reader_;
    istream &in_;
};

Image Decode(const std::string &filename) {
    std::ifstream in_file(filename, std::ios_base::binary);
    JpegDecoder jpeg_decoder(in_file);
    jpeg_decoder.Initialize();
    auto res_im = jpeg_decoder.Decode();
    in_file.close();
    return res_im;
}
}  // namespace Jpeg
