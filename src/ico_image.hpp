#pragma once

#include <vector>
#include <cstdint>

namespace ico {

class GrayImage {
public:
    uint8_t width() const { return m_width; }
    uint8_t height() const { return m_height; }
    std::vector<float> pixels() const { return m_pixels; }

private:
    GrayImage(uint8_t width, uint8_t height) : m_width{width}, m_height{height} {}
    uint8_t m_width;
    uint8_t m_height;
    std::vector<float> m_pixels;
    void pixels(const std::vector<float> &pixels) { m_pixels = pixels; }

public:
    static GrayImage from_icon(const std::vector<uint8_t> &buffer);
};

} // namespace ico