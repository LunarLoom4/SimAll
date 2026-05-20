// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/MovieEncoder.hpp
// Week   : 17
//
// Frame-sequence-to-movie encoder.  Three back-ends, all link-complete:
//
//   * APNG       — animated PNG, single-file output, supported by all
//                  modern browsers + most native viewers.  Pure C++; no
//                  external deps.  Lossless, suitable for residual plots
//                  / mesh-progress capture / "thumb-strip" previews.
//
//   * PPM-seq    — a folder of `frame_NNNNNN.ppm` ASCII P6 frames + a
//                  `ffmpeg.cmd` script the user can run to produce MP4.
//                  Use for very large captures where APNG would be too big.
//
//   * FFmpeg     — best-effort `popen` to a system `ffmpeg`, piping raw
//                  RGB frames over stdin to produce MP4.  Detected on
//                  begin(); falls back to APNG when ffmpeg is missing.
//
// The frame format passed in is always interleaved RGB8 — width*height*3
// bytes — so callers don't care which back-end will absorb them.
// =============================================================================
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace simall::io
{

enum class MovieCodec : std::uint8_t
{
    Auto,
    Apng,
    PpmSequence,
    FfmpegMp4
};

struct MovieOptions
{
    MovieCodec codec = MovieCodec::Auto;
    std::uint32_t fps = 30;
    bool loop = true;                  // APNG only; 0 == infinite
    std::string ffmpegPath = "ffmpeg"; // looked up on $PATH if unqualified
};

class MovieEncoder
{
public:
    MovieEncoder();
    ~MovieEncoder();
    MovieEncoder(const MovieEncoder&) = delete;
    MovieEncoder& operator=(const MovieEncoder&) = delete;

    bool begin(const std::string& path,
               std::uint32_t width,
               std::uint32_t height,
               MovieOptions opts = {});

    bool add_frame_rgb8(const std::uint8_t* pixels, std::size_t bytes);

    bool end();

    [[nodiscard]] MovieCodec active_codec() const noexcept;
    [[nodiscard]] std::uint32_t frame_count() const noexcept;
    [[nodiscard]] const std::string& error() const noexcept;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Internal helpers exposed for unit testing.
[[nodiscard]] std::vector<std::uint8_t> apng_encode(
    std::uint32_t w,
    std::uint32_t h,
    std::uint32_t fps,
    bool loop,
    const std::vector<std::vector<std::uint8_t>>& framesRgb8);

} // namespace simall::io
