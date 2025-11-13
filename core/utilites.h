#ifndef GONS_CORE_UTILITES_H_
#define GONS_CORE_UTILITES_H_

#include "config.h"
#include <algorithm>
#include <iostream>
#ifdef _WIN32
#include <stdio.h>
#include <windows.h>
void SetColor(int textColor, int bgColor = 0) {
  HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
  SetConsoleTextAttribute(hConsole, (bgColor << 4) | textColor);
}
#endif
namespace gons {
namespace utilites {
// LOG METHOD   need to be adapted according to platform
namespace LOG_MSG {

const char *const ERROR_PREFIX = "ERROR: ";
const char *const WARNING_PREFIX = "WARNING: ";

#define SHW(x) std::cout << x
#define ERR(x) std::cerr << x
#define WRN(x) std::cerr << x

#if (LOG_LVL == 0)

#define LOG(msg)
#define LOG_ERROR(msg)
#define LOG_WARNING(msg)

#elif (LOG_LVL == 1)
#ifdef __linux__
#define LOG(msg) SHW("\033[1;32m" << msg << "\033[0m\n");
#define LOG_ERROR(msg) ERR("\033[1;31m" << msg << "\033[0m\n");
#define LOG_WARNING(msg) WRN("\033[1;33m" << msg << "\033[0m\n");
#endif
#ifdef _WIN32
#define LOG(msg) SHW(msg);
#define LOG_ERROR(msg)                                                         \
  SetColor(4);                                                                 \
  SHW(ERROR_PREFIX << msg << "\n");                                            \
  SetColor(7);
#define LOG_WARNING(msg)                                                       \
  SetColor(6);                                                                 \
  SHW(WARNING_PREFIX << msg << "\n");                                          \
  SetColor(7);
#endif
#elif (LOG_LVL == 2)
#ifdef __linux__
#define LOG(msg)                                                               \
  SHW("\033[1;32m[" << __FILE__ << ": L. " << __LINE__ << "] " << msg          \
                    << "\033[0m\n");
#define LOG_ERROR(msg)                                                         \
  ERR("\033[1;31m[" << __FILE__ << ": L. " << __LINE__ << "] " << msg          \
                    << "\033[0m\n");
#define LOG_WARNING(msg)                                                       \
  WRN("\033[1;33m[" << __FILE__ << ": L. " << __LINE__ << "] " << msg          \
                    << "\033[0m\n");
#endif
#ifdef _WIN32
#define LOG(msg)                                                               \
  SetColor(4);                                                                 \
  SHW("[" << __FILE__ << ": L. " << __LINE__ << "] " << msg << "\n");          \
  SetColor(7);
#define LOG_ERROR(msg)                                                         \
  SetColor(6);                                                                 \
  SHW("[" << __FILE__ << ": L. " << __LINE__ << "] " << ERROR_PREFIX << msg    \
          << "\n");                                                            \
  SetColor(7);
#define LOG_WARNING(msg)                                                       \
  SetColor(9);                                                                 \
  SHW("[" << __FILE__ << ": L. " << __LINE__ << "] " << WARNING_PREFIX << msg  \
                    << "\n");                                                  \
  SetColor(7);
#endif
#else
#error "Log Level Setting Invalid"
#define LOG(msg)
#define LOG_ERROR(msg)
#define LOG_WARNING(msg)
#endif
#define CHECK(condition, msg)                                                  \
  if ((condition)) {                                                           \
    throw std::runtime_error(std::string(msg));                                \
  }

#define ASSERT(condition)                                                      \
  if ((condition)) {                                                           \
    LOG_ERROR("Assertion failed: " #condition);                                \
    std::abort();                                                              \
  }
#define UNUSED(x) (void)(x);

#define F_CHECK(condition, msg)                                                \
  if ((condition)) {                                                           \
    F_LOG_ERROR(msg);                                                          \
    std::abort();                                                              \
  }
} // namespace LOG_MSG

// Get Version
inline const char *GetVersion() { return GONS_VERSION; }

// Float Equal
#define FLT_EQUAL(a, b) (std::abs(a - b) < GONS_FLT_EPSILON)

// Float Not Equal
#define FLT_NOT_EQUAL(a, b) (std::abs(a - b) >= GONS_FLT_EPSILON)

// Float Equal to Zero
#define FLT_EQUAL_ZERO(a) (std::abs(a) < GONS_FLT_EPSILON)

} // namespace utilites
} // namespace gons

#endif