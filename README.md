# Groovy.Video

![Groovy.Video](groovy.video.jpg)

Create cool videos from your music.

## Description

Groovy.Video is a commandline tool to render a visuals video from an input
audio file. It uses offscreen OpenGL rendering and hardware accelerated encoding
of a MP4 video in VP9/Opus format.

## Generate Video Online

I've created an online service to generate a video at
[https://groovy.video](https://groovy.video)

If you'd like to support the project you can pay a small
amount to generate a Full HD video online.

## Building yourself

### Dependencies

* This has only been tested on [Arch Linux](https://archlinux.org).
* The visuals are generated by [projectM](https://github.com/projectM-visualizer/projectm) using the [Cream of the Crop](https://github.com/projectM-visualizer/presets-cream-of-the-crop) presets.
* It uses the VAAPI interface of [FFmpeg](https://ffmpeg.org/) to do hardware accelerated video encoding of VP9. It only works on Intel Graphics.
* It uses the [Mesa](https://www.mesa3d.org/) OpenGL drivers and [EGL](https://www.khronos.org/egl) libary to create an offscreen hardware accelerated rendering context.
* It uses [libva-cpp](https://github.com/GregoryIstratov/libav-cpp) to wrap the FFmpeg API.
* It uses [argparse](https://github.com/jamolnng/argparse) to parse commandline arguments.
* it Uses [CRC++](https://github.com/d-bahr/CRCpp) to encode the CRC32 of presets.

### Installing (Arch Linux)

```
pacman -S ffmpeg mesa libva-intel-driver intel-media-driver
git clone --recurse-submodules https://github.com/daansystems/groovy.video
make
```

### Executing program

* Please note on Intel architectures earlier than Icelake you need to set the LIBVA_DRIVER_NAME=i965 environment variable. Otherwise VP9 encoding will not be available.
```
Usage: ./groovy.video [options...]
Options:
    -i, --input            input audio file        (Required)
    -o, --output           output video file       (Required)
    -w, --width            video width (default 1280)
    -h, --height           video height (default 720)
    -d, --duration         preset duration in seconds (default 15)
    -p, --presets          presets (comma separated crc32 of filename)
    -l, --layer            filename of overlay image
    -m, --watermark        filename of watermark image
```

Example:

```
./groovy.video -i test/stereo-test.mp3 -o stereo-test.mp4 -p 2948676508 -m test/groovy.video.png
```
## Support

If you have an issue please create one on the Github project. If you'd like to contribute please create a Pull Request.

## Authors

[DaanSystems](https://www.daansystems.com)

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the GPLv2 License - see the LICENSE file for details
