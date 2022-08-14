#include <stdio.h>

extern "C" {
#include <EGL/egl.h>
#include <EGL/eglext.h>
#include <GL/gl.h>
}

struct EGLInternalData2 {
  bool m_isInitialized;

  int m_windowWidth;
  int m_windowHeight;
  int m_renderDevice{-1};

  EGLBoolean success;
  EGLint num_configs;
  EGLConfig egl_config;
  EGLSurface egl_surface;
  EGLContext egl_context;
  EGLDisplay egl_display;

  EGLInternalData2()
      : m_isInitialized(false), m_windowWidth(0), m_windowHeight(0) {}
};

EGLInternalData2 *initGL(int m_windowWidth, int m_windowHeight) {
  EGLBoolean success;
  EGLint num_configs;
  EGLConfig egl_config;
  EGLSurface egl_surface;
  EGLContext egl_context;
  EGLDisplay egl_display;

  EGLint egl_config_attribs[] = {EGL_RED_SIZE,
                                 8,
                                 EGL_GREEN_SIZE,
                                 8,
                                 EGL_BLUE_SIZE,
                                 8,
                                 EGL_DEPTH_SIZE,
                                 8,
                                 EGL_ALPHA_SIZE,
                                 8,
                                 EGL_SURFACE_TYPE,
                                 EGL_PBUFFER_BIT,
                                 EGL_RENDERABLE_TYPE,
                                 EGL_OPENGL_BIT,
                                 EGL_NONE};

  EGLint egl_pbuffer_attribs[] = {
      EGL_WIDTH, m_windowWidth, EGL_HEIGHT, m_windowHeight, EGL_NONE,
  };

  EGLInternalData2 *m_data = new EGLInternalData2();

  PFNEGLQUERYDEVICESEXTPROC eglQueryDevicesEXT =
      (PFNEGLQUERYDEVICESEXTPROC)eglGetProcAddress("eglQueryDevicesEXT");
  PFNEGLGETPLATFORMDISPLAYEXTPROC eglGetPlatformDisplayEXT =
      (PFNEGLGETPLATFORMDISPLAYEXTPROC)eglGetProcAddress(
          "eglGetPlatformDisplayEXT");
  PFNEGLQUERYDEVICESTRINGEXTPROC eglQueryDeviceStringEXT =
      (PFNEGLQUERYDEVICESTRINGEXTPROC)eglGetProcAddress(
          "eglQueryDeviceStringEXT");

  // Query EGL Devices
  const int max_devices = 32;
  EGLDeviceEXT egl_devices[max_devices];
  EGLint num_devices = 0;
  EGLint egl_error = eglGetError();
  if (!eglQueryDevicesEXT(max_devices, egl_devices, &num_devices) ||
      egl_error != EGL_SUCCESS) {
    fprintf(stderr, "eglQueryDevicesEXT Failed.\n");
    m_data->egl_display = EGL_NO_DISPLAY;
  }

  // Query EGL Screens
  if (m_data->m_renderDevice == -1) {
    // Chose default screen, by trying all
    for (EGLint i = 0; i < num_devices; ++i) {

      std::string deviceString =
          eglQueryDeviceStringEXT(egl_devices[i], EGL_EXTENSIONS);

      fprintf(stderr, "DEVICESTRING: %s\n", deviceString.c_str());

      // skip software renderers.
      if (deviceString.find("software") != std::string::npos) {
        continue;
      }
      // Set display
      EGLDisplay display = eglGetPlatformDisplayEXT(EGL_PLATFORM_DEVICE_EXT,
                                                    egl_devices[i], NULL);

      if (eglGetError() == EGL_SUCCESS && display != EGL_NO_DISPLAY) {
        int major, minor;
        EGLBoolean initialized = eglInitialize(display, &major, &minor);
        if (eglGetError() == EGL_SUCCESS && initialized == EGL_TRUE) {
          m_data->egl_display = display;
        }
      }
    }
  } else {
    // Chose specific screen, by using m_renderDevice
    if (m_data->m_renderDevice < 0 || m_data->m_renderDevice >= num_devices) {
      fprintf(stderr, "Invalid render_device choice: %d < %d.\n",
              m_data->m_renderDevice, num_devices);
      exit(EXIT_FAILURE);
    }

    // Set display
    EGLDisplay display = eglGetPlatformDisplayEXT(
        EGL_PLATFORM_DEVICE_EXT, egl_devices[m_data->m_renderDevice], NULL);
    if (eglGetError() == EGL_SUCCESS && display != EGL_NO_DISPLAY) {
      int major, minor;
      EGLBoolean initialized = eglInitialize(display, &major, &minor);
      if (eglGetError() == EGL_SUCCESS && initialized == EGL_TRUE) {
        m_data->egl_display = display;
      }
    }
  }

  if (!eglInitialize(m_data->egl_display, NULL, NULL)) {
    fprintf(stderr, "Unable to initialize EGL\n");
    exit(EXIT_FAILURE);
  }

  m_data->success = eglBindAPI(EGL_OPENGL_API);
  if (!m_data->success) {
    // TODO: Properly handle this error (requires change to default window
    // API to change return on all window types to bool).
    fprintf(stderr, "Failed to bind OpenGL API.\n");
    exit(EXIT_FAILURE);
  }

  m_data->success =
      eglChooseConfig(m_data->egl_display, egl_config_attribs,
                      &m_data->egl_config, 1, &m_data->num_configs);
  if (!m_data->success) {
    // TODO: Properly handle this error (requires change to default window
    // API to change return on all window types to bool).
    fprintf(stderr, "Failed to choose config (eglError: %d)\n", eglGetError());
    exit(EXIT_FAILURE);
  }
  if (m_data->num_configs != 1) {
    fprintf(stderr, "Didn't get exactly one config, but %d\n",
            m_data->num_configs);
    exit(EXIT_FAILURE);
  }

  m_data->egl_surface = eglCreatePbufferSurface(
      m_data->egl_display, m_data->egl_config, egl_pbuffer_attribs);
  if (m_data->egl_surface == EGL_NO_SURFACE) {
    fprintf(stderr, "Unable to create EGL surface (eglError: %d)\n",
            eglGetError());
    exit(EXIT_FAILURE);
  }

  m_data->egl_context = eglCreateContext(
      m_data->egl_display, m_data->egl_config, EGL_NO_CONTEXT, NULL);
  if (!m_data->egl_context) {
    fprintf(stderr, "Unable to create EGL context (eglError: %d)\n",
            eglGetError());
    exit(EXIT_FAILURE);
  }

  m_data->success = eglMakeCurrent(m_data->egl_display, m_data->egl_surface,
                                   m_data->egl_surface, m_data->egl_context);
  if (!m_data->success) {
    fprintf(stderr, "Failed to make context current (eglError: %d)\n",
            eglGetError());
    exit(EXIT_FAILURE);
  }

  fprintf(stderr, "GL_RENDERER = %s\n", (char *)glGetString(GL_RENDERER));
  fprintf(stderr, "GL_VERSION = %s\n", (char *)glGetString(GL_VERSION));
  fprintf(stderr, "GL_VENDOR = %s\n", (char *)glGetString(GL_VENDOR));
  fprintf(stderr, "GL_SHADING_LANGUAGE_VERSION = %s\n",
          (char *)glGetString(GL_SHADING_LANGUAGE_VERSION));
  fprintf(stderr, "GL_EXTENSIONS = %s\n", (char *)glGetString(GL_EXTENSIONS));

  return m_data;
}
