/**
    @file main.cpp
    @author DaanSystems
    @version 0.1
    @copyright
    @parblock
        Groovy.Video
        Copyright (c) 2022, DaanSystems
        All rights reserved.

        This file is licensed under the GPLv2 license.
    @endparblock
*/
#include <av/Resample.hpp>
#include <av/StreamReader.hpp>
#include <av/StreamWriter.hpp>
#include <chrono>
#include <iostream>

extern "C" {
#include <libavutil/audio_fifo.h>
#include <libavutil/timestamp.h>
}

#include <csignal>
#include <libprojectM/ProjectM.hpp>
#include <libprojectM/Renderer/SOIL2/stb_image.h>

#define CRCPP_USE_CPP11
#include "CRC.h"

#include "argparse.h"
#include "initgl.hpp"
using namespace argparse;

namespace av {
void writeLog(LogLevel level, internal::SourceLocation &&loc,
              std::string msg) noexcept {
  std::cerr << loc.toString() << ": " << msg << std::endl;
}
} // namespace av

template <typename... Args>
void println(std::string_view fmt, Args &&...args) noexcept {
  std::cerr << av::internal::format(fmt, std::forward<Args>(args)...)
            << std::endl;
}

template <typename Return>
Return assertExpected(av::Expected<Return> &&expected) noexcept {
  if (!expected) {
    std::cerr << " === Expected failure == \n"
              << expected.errorString() << std::endl;
    exit(EXIT_FAILURE);
  }

  if constexpr (std::is_same_v<Return, void>)
    return;
  else
    return expected.value();
}

void flipPixels(uint32_t *flippedPixels, uint32_t *imgPixels, int imgw,
                int imgh) {
  if (imgPixels != NULL) {
    for (int y = 0; y < imgh; y++) {
      for (int x = 0; x < imgw; x++) {
        flippedPixels[((imgh - y - 1) * imgw) + x] = imgPixels[(y * imgw) + x];
      }
    }
  }
}

inline uint32_t raw_blend(const uint32_t src, const uint32_t dest) {
  // setup and calculate α

  uint32_t src_a = src >> 24;
  uint32_t src_a_neg = 255 - src_a;
  uint32_t dest_a = dest >> 24;

  uint32_t res = src_a + ((dest_a * src_a_neg) >> 8);

  // setup and calculate R

  uint32_t src_r = (src >> 16) & 255;
  uint32_t dest_r = (dest >> 16) & 255;

  res = (res << 8) | (((src_r * src_a) + (dest_r * src_a_neg)) >> 8);

  // setup and calculate G

  uint32_t src_g = (src >> 8) & 255;
  uint32_t dest_g = (dest >> 8) & 255;

  res = (res << 8) | (((src_g * src_a) + (dest_g * src_a_neg)) >> 8);

  // setup and calculate B

  uint32_t src_b = src & 255;
  uint32_t dest_b = dest & 255;

  return (res << 8) | (((src_b * src_a) + (dest_b * src_a_neg)) >> 8);
}

void draw(const uint32_t *source, uint32_t *dest, const uint32_t canvasWidth,
          const uint32_t canvasHeight, const uint32_t width,
          const uint32_t height, const uint32_t targetx,
          const uint32_t targety) {
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      long tileIndex = (x + y * width);
      long canvasIndex = ((x + targetx) + (y + targety) * canvasWidth);

      dest[canvasIndex] = raw_blend(source[tileIndex], dest[canvasIndex]);
    }
  }
}

// trim from end of string (right)
inline std::string &rtrim(std::string &s, const char *t) {
  s.erase(s.find_last_not_of(t) + 1);
  return s;
}

// trim from beginning of string (left)
inline std::string &ltrim(std::string &s, const char *t) {
  s.erase(0, s.find_first_not_of(t));
  return s;
}

std::vector<uint32_t> split(std::string s, std::string delimiter) {
  std::vector<uint32_t> elems;
  char *token = std::strtok((char *)s.c_str(), delimiter.c_str());
  while (token) {
    elems.push_back(std::stoul(std::string(token)));
    token = std::strtok(NULL, ",");
  }
  return elems;
}

bool find(std::vector<uint32_t> haystack, uint32_t needle) {
  for (uint32_t item : haystack) {
    if (item == needle) {
      return true;
    }
  }
  return false;
}

std::string formatTime(int seconds) {

  // convret seconds into hr, min, sec
  int min = ((int)(seconds / 60)) % 60;
  int sec = (int)(seconds % 60);
  char buf[6];

  std::snprintf(buf, 6, "%02d:%02d", min, sec);

  return buf;
}

int main(int argc, const char **argv) {
  ArgumentParser parser(argv[0], "Generate visuals video from audio file");
  parser.add_argument("-i", "--input", "input audio file", true);
  parser.add_argument("-o", "--output", "output video file", true);
  parser.add_argument("-w", "--width", "video width (default 1280)", false);
  parser.add_argument("-h", "--height", "video height (default 720)", false);
  parser.add_argument("-d", "--duration", "preset duration in seconds (default 15)", false);
  parser.add_argument("-p", "--presets", "presets (comma separated crc32 of filename)", false);
  parser.add_argument("-l", "--layer", "filename of overlay image", false);
  parser.add_argument("-m", "--watermark", "filename of watermark image", false);

  // parser.enable_help();
  auto parseErr = parser.parse(argc, argv);
  if (parseErr) {
    std::cerr << parseErr << std::endl;
    parser.print_help();
    return -1;
  }

  std::vector<uint32_t> crcs;
  if (parser.exists("p")) {
    crcs = split(parser.get<std::string>("p"), ",");
  }

  auto startTime = std::chrono::steady_clock::now();
  double totalFrames;
  uint32_t frameCounter = 0;

  auto settings = new class ProjectM::Settings();
  settings->windowWidth = parser.get<unsigned int>("w") ?: 1280;
  settings->windowHeight = parser.get<unsigned int>("h") ?: 720;
  settings->textureSize = 2048;
  settings->meshX = 96;
  settings->meshY = 72;
  // settings->presetURL = "./linemilk";
  settings->presetURL = "./presets-cream-of-the-crop";
  settings->presetDuration = parser.get<unsigned int>("d") ?: 15;
  settings->softCutDuration = 4;
  settings->fps = 35;
  // settings->titleFontURL = "./fonts2/Vera.ttf";
  // settings->menuFontURL = "./fonts2/VeraMono.ttf";
  settings->datadir = ".";

  uint32_t bufsize =
      sizeof(uint32_t) * settings->windowWidth * settings->windowHeight;
  GLubyte *frameBuffer = (GLubyte *)malloc(bufsize);
  if (!frameBuffer) {
    fprintf(stderr, "Alloc image frameBuffer failed!\n");
    return -1;
  }

  EGLInternalData2 *context =
      initGL(settings->windowWidth, settings->windowHeight);

  ProjectM *p = new ProjectM(*settings);

  auto &pcm = p->Pcm();
  if (crcs.size()) {
    for (int i = 0; i < p->PlaylistSize(); i++) {
      auto preset = p->PresetURL(i);
      preset = rtrim(preset, ".milk");
      preset = ltrim(preset, settings->presetURL.c_str());
      std::uint32_t crc =
          CRC::Calculate(preset.data(), preset.size(), CRC::CRC_32());
      if (!find(crcs, crc)) {
        p->RemovePreset(i);
        i--;
      }
    }
  }
  std::cerr << "presets read..." << std::endl;
  fprintf(stderr, "Playlist size: %d\n", p->PlaylistSize());
  av_log_set_level(AV_LOG_VERBOSE);

  p->SelectRandom(true);
  if (p->PlaylistSize() == 1) {
    p->SetPresetLocked(true);
  }

  {
    auto reader = assertExpected(
        av::StreamReader::create(parser.get<std::string>("i"), true));
    auto writer =
        assertExpected(av::StreamWriter::create(parser.get<std::string>("o")));
    {
      auto framerate = av_make_q(settings->fps, 1);
      framerate = av_inv_q(framerate);
      // av::OptValueMap codecOpts = {{"preset", "13"}, {"crf", 29}};
      // av::OptValueMap codecOpts = {{"preset", 3}};
      av::OptValueMap codecOpts = {{"global_quality", "70"}, {"bf", "0"}};

      // assertExpected(writer->addVideoStream("libsvtav1",
      // settings->windowWidth, settings->windowHeight, AV_PIX_FMT_RGBA,
      // framerate, std::move(codecOpts)));

      // assertExpected(writer->addVideoStream(
      //     "libvpx-vp9", settings->windowWidth, settings->windowHeight,
      //     AV_PIX_FMT_RGBA, framerate, std::move(codecOpts)));

      assertExpected(writer->addVideoStream(
          "vp9_vaapi", settings->windowWidth, settings->windowHeight,
          AV_PIX_FMT_RGBA, framerate, std::move(codecOpts)));
    }

    {
      auto bitrate = 128 * 1000;
      assertExpected(writer->addAudioStream(
          "libopus", reader->channels(), AV_SAMPLE_FMT_FLT, 48000,
          reader->channels(), AV_SAMPLE_FMT_FLT, 48000, bitrate));
    }
    assertExpected(writer->open());

    auto audioFrame = av::Frame();
    double totalTime = 0;
    AVAudioFifo *fifo =
        av_audio_fifo_alloc(reader->sampleFormat(), reader->channels(), 2048);

    while (true) {
      if (!assertExpected(reader->readFrame(audioFrame))) {
        fprintf(stderr, "No more frames.\n");
        break;
      }

      if (audioFrame.type() == AVMEDIA_TYPE_AUDIO) {
        auto audioF = audioFrame.native();
        totalTime +=
            double(audioF->nb_samples) / double(reader->audio()->time_base.den);
        auto err = av_audio_fifo_write(fifo, (void **)audioF->data,
                                       audioF->nb_samples);
        if (err < 0) {
          LOG_AV_ERROR("Error allocating split frames: {}",
                       av::avErrorStr(err));
          return false;
        }
      } else {
        std::cerr << "TYPE:" << audioFrame.type() << std::endl;
      }
    }

    auto totalSamples = av_audio_fifo_size(fifo);

    fprintf(stderr, "Total samples: %d\n", totalSamples);
    auto superFrame = av_frame_alloc();
    superFrame->nb_samples = totalSamples;
    superFrame->format = reader->sampleFormat();
    superFrame->channel_layout = reader->channel_layout();
    superFrame->channels = reader->channels();
    auto err = av_frame_get_buffer(superFrame, 0);
    if (err < 0) {
      LOG_AV_ERROR("Error allocating superFrame: {}", av::avErrorStr(err));
      return false;
    }

    int totalSamplesRead =
        av_audio_fifo_read(fifo, (void **)superFrame->data, totalSamples);
    fprintf(stderr, "Total samples read: %jd duration: %jd\n", totalSamplesRead,
            superFrame->pkt_duration);

    unsigned int lastPreset = -1;

    auto swrExp = assertExpected(av::Resample::create(
        reader->channels(), reader->sampleFormat(), reader->sampleRate(),
        reader->channels(), AV_SAMPLE_FMT_FLT, 48000));

    uint64_t outSamples = swr_get_out_samples(swrExp->native(), totalSamples);

    fprintf(stderr, "insamples: %jd outsamples: %jd\n", totalSamples,
            outSamples);

    auto superOpusFrame = av_frame_alloc();
    superOpusFrame->nb_samples = outSamples;
    superOpusFrame->format = AV_SAMPLE_FMT_FLT;
    superOpusFrame->channel_layout = reader->channel_layout();
    superOpusFrame->channels = reader->channels();
    err = av_frame_get_buffer(superOpusFrame, 0);
    if (err < 0) {
      LOG_AV_ERROR("Error allocating superOpusFrame: {}", av::avErrorStr(err));
      return false;
    }

    swrExp->convertData((uint8_t **)superOpusFrame->data, outSamples,
                        (const uint8_t **)superFrame->data, totalSamples);

    av_frame_free(&superFrame);
    av_audio_fifo_free(fifo);

    totalFrames = totalTime * settings->fps;

    uint64_t floatsPerVideoFrame = double(outSamples) / totalFrames;

    fprintf(stderr, "Total time: %f\n", totalTime);
    fprintf(stderr, "Total Frames needed: %f\n", totalFrames);
    fprintf(stderr, "Floats per frame needed: %d\n", floatsPerVideoFrame);

    uint64_t frameSize = writer->stream(1)->encoder->native()->frame_size;
    fprintf(stderr, "Audio frame size: %jd\n", frameSize);
    int audioPos = 0;

    int overlayWidth;
    int overlayHeight;
    int overlayChannels;
    unsigned char *overlay = NULL;
    auto overlayerPath = parser.get<std::string>("l");

    if (!overlayerPath.empty()) {
      overlay = stbi_load(overlayerPath.c_str(), &overlayWidth, &overlayHeight,
                          &overlayChannels, 4);
      if (overlay == NULL) {
        fprintf(stderr, "Error in loading the overlay image\n");
        exit(1);
      }
      fprintf(stderr,
              "Loaded overlay image with a width of %dpx, a height of %dpx and "
              "%d channels\n",
              overlayWidth, overlayHeight, overlayChannels);

      if (overlayWidth > settings->windowWidth ||
          overlayHeight > settings->windowHeight) {
        fprintf(stderr, "Overlayer is larger than window size %d x %d\n",
                settings->windowWidth, settings->windowHeight);
        exit(1);
      }
    }

    int watermarkWidth;
    int watermarkHeight;
    int watermarkChannels;
    unsigned char *watermark = NULL;

    auto watermarkPath = parser.get<std::string>("m");

    if (!watermarkPath.empty()) {
      watermark = stbi_load(watermarkPath.c_str(), &watermarkWidth,
                            &watermarkHeight, &watermarkChannels, 4);
      if (watermark == NULL) {
        fprintf(stderr, "Error in loading the watermark image\n");
        exit(1);
      }
      fprintf(stderr,
              "Loaded image with a width of %dpx, a height of %dpx and %d "
              "channels\n",
              watermarkWidth, watermarkHeight, watermarkChannels);
    }

    auto frame = assertExpected(av::Frame::create(
        settings->windowWidth, settings->windowHeight, AV_PIX_FMT_RGBA));
    auto f = frame.get();
    f->type(AVMEDIA_TYPE_VIDEO);

    for (uint64_t i = 0; i < outSamples; i += floatsPerVideoFrame) {
      if (av_frame_make_writable(f->native()) < 0) {
        LOG_AV_ERROR("Error making frame writable: {}", av::avErrorStr(err));
        return false;
      }
      unsigned int currentPreset = 0;
      if (p->SelectedPresetIndex(currentPreset) &&
          currentPreset != lastPreset) {
        std::cerr << "Preset at " << formatTime(frameCounter / settings->fps)
                  << " " << p->PresetURL(currentPreset) << std::endl;
        lastPreset = currentPreset;
      }

      size_t maxSamples = std::min(outSamples - i, floatsPerVideoFrame);

      pcm.AddStereo((float *)(superOpusFrame->data[0]) + i * 2, maxSamples);
      p->RenderFrame();

      glReadPixels(0, 0, settings->windowWidth, settings->windowHeight, GL_RGBA,
                   GL_UNSIGNED_BYTE, frameBuffer);
      flipPixels((uint32_t *)f->native()->data[0], (uint32_t *)frameBuffer,
                 settings->windowWidth, settings->windowHeight);

      if (overlay && frameCounter < settings->fps * 5) {
        draw((uint32_t *)overlay, (uint32_t *)f->native()->data[0],
             settings->windowWidth, settings->windowHeight, overlayWidth,
             overlayHeight, (settings->windowWidth - overlayWidth) / 2,
             (settings->windowHeight - overlayHeight) / 2);
      }

      if (watermark) {
        draw((uint32_t *)watermark, (uint32_t *)f->native()->data[0],
             settings->windowWidth, settings->windowHeight, watermarkWidth,
             watermarkHeight, settings->windowWidth - watermarkWidth,
             settings->windowHeight - watermarkHeight);
      }
      assertExpected(writer->write(*f, 0));

      for (int k = 0;
           k <= floatsPerVideoFrame / frameSize && audioPos < outSamples; k++) {
        auto opusFrame = assertExpected(
            writer->stream(1)->encoder->newWriteableAudioFrame());
        auto opusF = opusFrame->native();

        size_t maxSamples = std::min(outSamples - audioPos, frameSize);
        opusF->nb_samples = maxSamples;
        auto err = av_frame_get_buffer(opusF, 0);
        if (err < 0) {
          LOG_AV_ERROR("Error allocating split frames: {}",
                       av::avErrorStr(err));
          return false;
        }
        err = av_samples_copy(opusF->data, superOpusFrame->data, 0, audioPos,
                              opusF->nb_samples, opusF->channels,
                              writer->stream(1)->encoder->native()->sample_fmt);
        if (err < 0) {
          LOG_AV_ERROR("Error copying samples: {}", av::avErrorStr(err));
          return false;
        }
        assertExpected(writer->write(*(opusFrame.get()), 1));
        audioPos += opusF->nb_samples;
      }
      /* if (frameCounter % settings->fps == 0) {
        fprintf(stderr, "frame rendered %d / %d = %d %%\n", frameCounter,
                uint32_t(totalFrames),
                int(100.0 * (double(frameCounter) / totalFrames)));
      } */
      frameCounter++;
    }
    av_frame_free(&superOpusFrame);
    if (watermark) {
      stbi_image_free(watermark);
    }
    if (overlay) {
      stbi_image_free(overlay);
    }
  }

  free(frameBuffer);
  delete p;
  delete settings;

  auto duration = std::chrono::duration_cast<std::chrono::seconds>(
                      std::chrono::steady_clock::now() - startTime)
                      .count();
  std::cerr << "Done elapsed time in seconds: " << duration << "."
            << " Total frames rendered: " << frameCounter << "."
            << " Avg frames/sec: " << int(totalFrames / duration) << "."
            << std::endl;
  return 0;
}
