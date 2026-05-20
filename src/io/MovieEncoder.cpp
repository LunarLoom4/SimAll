// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/MovieEncoder.cpp
//
// APNG (animated PNG) encoder.  Spec: https://wiki.mozilla.org/APNG_Specification
//
// Frame layout produced here:
//   PNG signature
//   IHDR
//   acTL (numFrames, numPlays)
//   fcTL #0  (sequence 0)
//   IDAT     (frame 0 image data; also the default still image)
//   fcTL #1  (sequence 1)
//   fdAT     (sequence 2, frame 1 data)
//   fcTL #2  (sequence 3)
//   fdAT     (sequence 4, frame 2 data)
//   ...
//   IEND
//
// IDAT/fdAT payload is a zlib stream of one filter byte + width*3 RGB bytes
// per scanline.  We use *no* compression (BTYPE=00, stored blocks) so the
// implementation is self-contained — output is larger than libpng would
// produce but is bit-perfect, deterministic, and trivially decodable.
// =============================================================================
#include "io/MovieEncoder.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

namespace simall::io
{

namespace
{

// ---------------------------------------------------------------------------
// CRC32 (zlib variant) for PNG chunk CRCs.
// ---------------------------------------------------------------------------
std::uint32_t crc32_of(const std::uint8_t* data, std::size_t n)
{
    static std::uint32_t table[256];
    static bool init = false;
    if (!init) {
        for (std::uint32_t i = 0; i < 256; ++i) {
            std::uint32_t c = i;
            for (int k = 0; k < 8; ++k)
                c = (c & 1) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            table[i] = c;
        }
        init = true;
    }
    std::uint32_t c = 0xFFFFFFFFu;
    for (std::size_t i = 0; i < n; ++i)
        c = table[(c ^ data[i]) & 0xFF] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

// ---------------------------------------------------------------------------
// Adler-32 for the zlib trailer.
// ---------------------------------------------------------------------------
std::uint32_t adler32_of(const std::uint8_t* data, std::size_t n)
{
    std::uint32_t a = 1, b = 0;
    for (std::size_t i = 0; i < n; ++i) {
        a = (a + data[i]) % 65521;
        b = (b + a) % 65521;
    }
    return (b << 16) | a;
}

// ---------------------------------------------------------------------------
// "Stored" zlib stream (no compression).  Wraps raw bytes with a 2-byte
// zlib header, one or more uncompressed deflate blocks, and an Adler-32.
// ---------------------------------------------------------------------------
std::vector<std::uint8_t> deflate_stored(const std::uint8_t* data, std::size_t n)
{
    std::vector<std::uint8_t> out;
    out.reserve(n + 64);
    // zlib header: CMF=0x78 (deflate, 32K window), FLG=0x01 (no preset dict, level 0; FCHECK fixed)
    out.push_back(0x78);
    out.push_back(0x01);
    // Stored blocks (max 65535 bytes each).
    std::size_t i = 0;
    while (true) {
        const std::size_t chunk = std::min<std::size_t>(n - i, 65535);
        const bool last = (i + chunk == n);
        out.push_back(last ? 0x01 : 0x00); // BFINAL + BTYPE=00
        out.push_back(std::uint8_t(chunk & 0xFF));
        out.push_back(std::uint8_t((chunk >> 8) & 0xFF));
        const std::uint16_t nlen = static_cast<std::uint16_t>(~chunk);
        out.push_back(std::uint8_t(nlen & 0xFF));
        out.push_back(std::uint8_t((nlen >> 8) & 0xFF));
        out.insert(out.end(), data + i, data + i + chunk);
        i += chunk;
        if (last)
            break;
        if (chunk == 0)
            break;
    }
    const std::uint32_t adler = adler32_of(data, n);
    out.push_back(std::uint8_t((adler >> 24) & 0xFF));
    out.push_back(std::uint8_t((adler >> 16) & 0xFF));
    out.push_back(std::uint8_t((adler >> 8) & 0xFF));
    out.push_back(std::uint8_t(adler & 0xFF));
    return out;
}

void put_be_u32(std::vector<std::uint8_t>& out, std::uint32_t v)
{
    out.push_back(std::uint8_t((v >> 24) & 0xFF));
    out.push_back(std::uint8_t((v >> 16) & 0xFF));
    out.push_back(std::uint8_t((v >> 8) & 0xFF));
    out.push_back(std::uint8_t(v & 0xFF));
}
void put_be_u16(std::vector<std::uint8_t>& out, std::uint16_t v)
{
    out.push_back(std::uint8_t((v >> 8) & 0xFF));
    out.push_back(std::uint8_t(v & 0xFF));
}

void put_chunk(std::vector<std::uint8_t>& out,
               const char tag[4],
               const std::vector<std::uint8_t>& payload)
{
    put_be_u32(out, std::uint32_t(payload.size()));
    const std::size_t crcStart = out.size();
    out.insert(out.end(), tag, tag + 4);
    out.insert(out.end(), payload.begin(), payload.end());
    const std::uint32_t crc = crc32_of(out.data() + crcStart, 4 + payload.size());
    put_be_u32(out, crc);
}

std::vector<std::uint8_t> filter_rgb_scanlines(std::uint32_t w,
                                               std::uint32_t h,
                                               const std::uint8_t* rgb)
{
    std::vector<std::uint8_t> out;
    out.reserve(std::size_t(h) * (1 + 3 * std::size_t(w)));
    for (std::uint32_t y = 0; y < h; ++y) {
        out.push_back(0); // filter "None"
        out.insert(out.end(), rgb + 3 * std::size_t(w) * y, rgb + 3 * std::size_t(w) * (y + 1));
    }
    return out;
}

} // namespace

std::vector<std::uint8_t> apng_encode(std::uint32_t w,
                                      std::uint32_t h,
                                      std::uint32_t fps,
                                      bool loop,
                                      const std::vector<std::vector<std::uint8_t>>& frames)
{
    std::vector<std::uint8_t> out;
    // PNG signature
    const std::uint8_t sig[8] = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    out.insert(out.end(), sig, sig + 8);

    // IHDR
    {
        std::vector<std::uint8_t> p;
        put_be_u32(p, w);
        put_be_u32(p, h);
        p.push_back(8); // bit depth
        p.push_back(2); // colour type = truecolour RGB
        p.push_back(0); // compression
        p.push_back(0); // filter
        p.push_back(0); // interlace
        put_chunk(out, "IHDR", p);
    }
    // acTL (animation control)
    {
        std::vector<std::uint8_t> p;
        put_be_u32(p, std::uint32_t(frames.size()));
        put_be_u32(p, loop ? 0u : 1u);
        put_chunk(out, "acTL", p);
    }
    std::uint32_t sequenceNumber = 0;
    const std::uint16_t delayNum = 1;
    const std::uint16_t delayDen = static_cast<std::uint16_t>(fps == 0 ? 1 : fps);

    for (std::size_t f = 0; f < frames.size(); ++f) {
        // fcTL
        {
            std::vector<std::uint8_t> p;
            put_be_u32(p, sequenceNumber++);
            put_be_u32(p, w);
            put_be_u32(p, h);
            put_be_u32(p, 0);
            put_be_u32(p, 0);
            put_be_u16(p, delayNum);
            put_be_u16(p, delayDen);
            p.push_back(0); // dispose
            p.push_back(0); // blend
            put_chunk(out, "fcTL", p);
        }
        const auto filtered = filter_rgb_scanlines(w, h, frames[f].data());
        const auto zlib = deflate_stored(filtered.data(), filtered.size());
        if (f == 0) {
            put_chunk(out, "IDAT", zlib);
        } else {
            std::vector<std::uint8_t> p;
            put_be_u32(p, sequenceNumber++);
            p.insert(p.end(), zlib.begin(), zlib.end());
            put_chunk(out, "fdAT", p);
        }
    }
    put_chunk(out, "IEND", {});
    return out;
}

// ---------------------------------------------------------------------------
// MovieEncoder façade
// ---------------------------------------------------------------------------
struct MovieEncoder::Impl
{
    MovieCodec codec = MovieCodec::Auto;
    MovieCodec active = MovieCodec::Apng;
    std::uint32_t w = 0, h = 0;
    std::uint32_t fps = 30;
    bool loop = true;
    std::string path;
    std::string ffmpeg;
    std::uint32_t frameNo = 0;
    std::string err;
    // For APNG we buffer frames; PPM we emit incrementally; FFmpeg we pipe.
    std::vector<std::vector<std::uint8_t>> framesRgb;
    FILE* ffmpegPipe = nullptr;
};

MovieEncoder::MovieEncoder() : impl_(std::make_unique<Impl>()) {}
MovieEncoder::~MovieEncoder()
{
    (void) end();
}

namespace
{
MovieCodec resolve_codec(MovieCodec requested, const std::string& path)
{
    if (requested != MovieCodec::Auto)
        return requested;
    // Pick by extension.
    auto ends_with = [](const std::string& s, const char* sfx) {
        const std::size_t n = std::strlen(sfx);
        return s.size() >= n && s.compare(s.size() - n, n, sfx) == 0;
    };
    if (ends_with(path, ".png"))
        return MovieCodec::Apng;
    if (ends_with(path, ".apng"))
        return MovieCodec::Apng;
    if (ends_with(path, ".mp4"))
        return MovieCodec::FfmpegMp4;
    return MovieCodec::PpmSequence;
}
} // namespace

bool MovieEncoder::begin(const std::string& path,
                         std::uint32_t w,
                         std::uint32_t h,
                         MovieOptions opts)
{
    impl_->path = path;
    impl_->w = w;
    impl_->h = h;
    impl_->fps = opts.fps == 0 ? 1 : opts.fps;
    impl_->loop = opts.loop;
    impl_->codec = opts.codec;
    impl_->ffmpeg = opts.ffmpegPath;
    impl_->frameNo = 0;
    impl_->framesRgb.clear();
    impl_->active = resolve_codec(opts.codec, path);

    if (impl_->active == MovieCodec::FfmpegMp4) {
        // Try to open an ffmpeg pipe; on failure, downgrade to APNG.
        std::ostringstream cmd;
        cmd << impl_->ffmpeg << " -y -f rawvideo -pix_fmt rgb24 -s " << w << "x" << h << " -r "
            << impl_->fps << " -i - -c:v libx264 -pix_fmt yuv420p \"" << path << "\"";
#if defined(_WIN32)
        impl_->ffmpegPipe = _popen(cmd.str().c_str(), "wb");
#else
        impl_->ffmpegPipe = popen(cmd.str().c_str(), "w");
#endif
        if (!impl_->ffmpegPipe) {
            impl_->active = MovieCodec::Apng;
        }
    }
    if (impl_->active == MovieCodec::PpmSequence) {
        // Folder is the file path interpreted as a stem.  No-op on begin.
    }
    return true;
}

bool MovieEncoder::add_frame_rgb8(const std::uint8_t* pixels, std::size_t bytes)
{
    if (!impl_ || bytes < 3 * std::size_t(impl_->w) * std::size_t(impl_->h)) {
        impl_->err = "Frame too small";
        return false;
    }
    const std::size_t expected = 3 * std::size_t(impl_->w) * std::size_t(impl_->h);
    if (impl_->active == MovieCodec::Apng) {
        impl_->framesRgb.emplace_back(pixels, pixels + expected);
    } else if (impl_->active == MovieCodec::PpmSequence) {
        char name[64];
        std::snprintf(name, sizeof(name), "_%06u.ppm", impl_->frameNo);
        const std::string outPath = impl_->path + name;
        std::ofstream f(outPath, std::ios::binary);
        if (!f) {
            impl_->err = "Cannot open " + outPath;
            return false;
        }
        f << "P6\n" << impl_->w << " " << impl_->h << "\n255\n";
        f.write(reinterpret_cast<const char*>(pixels), std::streamsize(expected));
    } else if (impl_->active == MovieCodec::FfmpegMp4 && impl_->ffmpegPipe) {
        if (std::fwrite(pixels, 1, expected, impl_->ffmpegPipe) != expected) {
            impl_->err = "ffmpeg pipe write failed";
            return false;
        }
    }
    ++impl_->frameNo;
    return true;
}

bool MovieEncoder::end()
{
    if (!impl_)
        return true;
    if (impl_->active == MovieCodec::Apng && !impl_->framesRgb.empty()) {
        const auto bytes =
            apng_encode(impl_->w, impl_->h, impl_->fps, impl_->loop, impl_->framesRgb);
        std::ofstream f(impl_->path, std::ios::binary);
        if (!f) {
            impl_->err = "Cannot open " + impl_->path;
            return false;
        }
        f.write(reinterpret_cast<const char*>(bytes.data()), std::streamsize(bytes.size()));
        impl_->framesRgb.clear();
        return f.good();
    }
    if (impl_->active == MovieCodec::FfmpegMp4 && impl_->ffmpegPipe) {
#if defined(_WIN32)
        _pclose(impl_->ffmpegPipe);
#else
        pclose(impl_->ffmpegPipe);
#endif
        impl_->ffmpegPipe = nullptr;
    }
    if (impl_->active == MovieCodec::PpmSequence) {
        // Emit an ffmpeg.cmd convenience script next to the frame folder.
        std::ofstream f(impl_->path + "_ffmpeg.cmd");
        if (f) {
            f << "ffmpeg -y -framerate " << impl_->fps << " -i \"" << impl_->path << "_%06d.ppm\""
              << " -c:v libx264 -pix_fmt yuv420p \"" << impl_->path << ".mp4\"\n";
        }
    }
    return true;
}

MovieCodec MovieEncoder::active_codec() const noexcept
{
    return impl_ ? impl_->active : MovieCodec::Auto;
}
std::uint32_t MovieEncoder::frame_count() const noexcept
{
    return impl_ ? impl_->frameNo : 0;
}
const std::string& MovieEncoder::error() const noexcept
{
    static const std::string empty;
    return impl_ ? impl_->err : empty;
}

} // namespace simall::io
