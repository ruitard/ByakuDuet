#include <vector>
#include <cstdint>
#include <stdexcept>
#include <cstring>
#include <cassert>
#include <iostream>

#include "ico_image.hpp"

namespace {

uint8_t read_byte(const std::vector<uint8_t> &bytes, size_t offset) { return bytes.at(offset); }

uint16_t read_uint16(const std::vector<uint8_t> &bytes, size_t offset) {
    uint16_t u0 = read_byte(bytes, offset);
    uint16_t u1 = read_byte(bytes, offset + 1);
    return u0 | (u1 << 8);
}

uint16_t read_uint32(const std::vector<uint8_t> &bytes, size_t offset) {
    uint32_t u0 = read_byte(bytes, offset);
    uint32_t u1 = read_byte(bytes, offset + 1);
    uint32_t u2 = read_byte(bytes, offset + 2);
    uint32_t u3 = read_byte(bytes, offset + 3);

    return u0 | (u1 << 8) | (u2 << 16) | (u3 << 24);
}

template <typename T> void require_equal(T a, T b) {
    if (a != b) {
        throw std::logic_error("not supported format");
    }
}

} // namespace

namespace ico {

static constexpr size_t BMP_INFO_HEADER_SIZE = 40;

struct DirHeader {
    uint16_t reserved;
    uint16_t type;
    uint16_t count;
};

struct ImageHeader {
    uint8_t width;
    uint8_t height;
    uint8_t color_palette;
    uint8_t reserved;
    uint16_t planes;
    uint16_t bit_count;
    uint32_t size;
    uint16_t offset;
};

GrayImage GrayImage::from_icon(const std::vector<uint8_t> &buffer) {
    if (buffer.size() < sizeof(DirHeader) + sizeof(ImageHeader)) {
        throw std::invalid_argument("");
    }
    if (std::memcmp(buffer.data(), "\x00\x00\x01\x00", 4) != 0) {
        throw std::invalid_argument("not ico format");
    }

    uint8_t width = read_byte(buffer, sizeof(DirHeader) + offsetof(ImageHeader, width));
    uint8_t height = read_byte(buffer, sizeof(DirHeader) + offsetof(ImageHeader, height));
    GrayImage image(width, height);

    require_equal(read_byte(buffer, sizeof(DirHeader) + offsetof(ImageHeader, color_palette)), uint8_t{});
    require_equal(read_uint16(buffer, sizeof(DirHeader) + offsetof(ImageHeader, bit_count)), uint16_t{32});

    uint32_t bmp_size = read_uint32(buffer, sizeof(DirHeader) + offsetof(ImageHeader, size));
    uint16_t bmp_offset = read_uint16(buffer, sizeof(DirHeader) + offsetof(ImageHeader, offset));
    require_equal((bmp_offset + bmp_size) <= buffer.size(), true);

    require_equal((BMP_INFO_HEADER_SIZE + width * height * 4) <= bmp_size, true);

    require_equal(std::memcmp(buffer.data() + bmp_offset, "\x28\x00\x00\x00", 4), 0);

    std::vector<float> pixels;
    pixels.reserve(width * height);
    for (uint32_t i = 1; i <= height; i++) {
        for (uint32_t j = 0; j < width; j++) {
            uint32_t anchor = bmp_offset + BMP_INFO_HEADER_SIZE + ((height - i) * width + j) * 4;
            pixels.push_back(static_cast<float>(buffer[anchor] + buffer[anchor + 1] + buffer[anchor + 2]) / 3);
        }
    }

    image.pixels(pixels);
    return image;
}

} // namespace ico