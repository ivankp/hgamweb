#ifndef IVANP_ERROR_HH
#define IVANP_ERROR_HH

#include <stdexcept>
#include <cerrno>

#include "macros.hh"
#include "string.hh"

#define IVANP_ERROR_PREF __FILE__ ":" STR(__LINE__) ": "
#define THROW(...) throw ivanp::error( \
  IVANP_ERROR_PREF, __VA_ARGS__ )
#define THROW_ERRNO(F,...) throw ivanp::error( \
  IVANP_ERROR_PREF F "(", __VA_ARGS__, "): ", std::strerror(errno) )

namespace ivanp {

struct error: std::runtime_error {
  using std::runtime_error::runtime_error;
  template <typename... T> [[ gnu::always_inline ]]
  error(T&&... x): std::runtime_error(cat(std::forward<T>(x)...)) { };
  [[ gnu::always_inline ]]
  error(const char* str): std::runtime_error(str) { };
};

}

#endif
