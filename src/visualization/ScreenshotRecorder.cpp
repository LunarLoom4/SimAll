// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ScreenshotRecorder.cpp
// =============================================================================
#include "visualization/ScreenshotRecorder.hpp"

#include <cstdint>
#include <fstream>
#include <sstream>
#include <vector>

namespace simall::visualization {

bool ScreenshotRecorder::save_ppm(const Image& image,
                                   const std::filesystem::path& path)
{
    if (image.empty()) return false;
    std::ofstream out(path, std::ios::binary);
    if (!out) return false;
    out << "P6\n" << image.width << ' ' << image.height << "\n255\n";
    std::vector<std::uint8_t> rgb;
    rgb.reserve(static_cast<std::size_t>(image.width) * image.height * 3u);
    const std::size_t n = static_cast<std::size_t>(image.width) * image.height;
    for (std::size_t i = 0; i < n; ++i) {
        rgb.push_back(image.pixels[i*4 + 0]);
        rgb.push_back(image.pixels[i*4 + 1]);
        rgb.push_back(image.pixels[i*4 + 2]);
    }
    out.write(reinterpret_cast<const char*>(rgb.data()),
              static_cast<std::streamsize>(rgb.size()));
    return out.good();
}

bool ScreenshotRecorder::save_rgba_raw(const Image& image,
                                        const std::filesystem::path& path)
{
    if (image.empty()) return false;
    std::ofstream out(path, std::ios::binary);
    if (!out) return false;
    const std::uint32_t w = image.width, h = image.height;
    out.write(reinterpret_cast<const char*>(&w), sizeof(w));
    out.write(reinterpret_cast<const char*>(&h), sizeof(h));
    out.write(reinterpret_cast<const char*>(image.pixels.data()),
              static_cast<std::streamsize>(image.pixels.size()));
    return out.good();
}

Image ScreenshotRecorder::load_ppm(const std::filesystem::path& path)
{
    Image img;
    std::ifstream in(path, std::ios::binary);
    if (!in) return img;
    std::string magic;
    int w = 0, h = 0, maxVal = 0;
    in >> magic;
    if (magic != "P6") return img;
    // Skip comments and whitespace.
    auto skip_ws = [&]() {
        while (in.good()) {
            int c = in.peek();
            if (c == '#') { std::string s; std::getline(in, s); }
            else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') in.get();
            else return;
        }
    };
    skip_ws();  in >> w;
    skip_ws();  in >> h;
    skip_ws();  in >> maxVal;
    in.get();   // consume single whitespace
    if (w <= 0 || h <= 0 || maxVal != 255) return img;
    img.resize(static_cast<std::uint32_t>(w), static_cast<std::uint32_t>(h));
    const std::size_t n = static_cast<std::size_t>(w) * h;
    std::vector<std::uint8_t> rgb(n * 3);
    in.read(reinterpret_cast<char*>(rgb.data()),
            static_cast<std::streamsize>(rgb.size()));
    if (!in) { img.pixels.clear(); img.width = img.height = 0; return img; }
    for (std::size_t i = 0; i < n; ++i) {
        img.pixels[i*4 + 0] = rgb[i*3 + 0];
        img.pixels[i*4 + 1] = rgb[i*3 + 1];
        img.pixels[i*4 + 2] = rgb[i*3 + 2];
        img.pixels[i*4 + 3] = 255;
    }
    return img;
}

}  // namespace simall::visualization
