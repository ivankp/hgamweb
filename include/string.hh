#ifndef IVANP_STRING_HH
#define IVANP_STRING_HH

#include <string>
#include <cstring>
#include <string_view>
#include <type_traits>

namespace ivanp {

[[nodiscard]]
inline const char* cstr(const char* x) noexcept { return x; }
[[nodiscard]]
inline const char* cstr(const std::string& x) noexcept { return x.c_str(); }

namespace impl {

template <typename T>
[[nodiscard]]
[[gnu::always_inline]]
inline std::string_view to_string_view(const T& x) noexcept {
  return x;
}
[[nodiscard]]
[[gnu::always_inline]]
inline std::string_view to_string_view(const char& x) noexcept {
  return { &x, 1 };
}

template <typename T>
inline constexpr bool is_sv = std::is_same_v<T,std::string_view>;

}

[[nodiscard]]
inline std::string cat() noexcept { return { }; }
[[nodiscard]]
inline std::string cat(std::string x) noexcept { return x; }
[[nodiscard]]
inline std::string cat(const char* x) noexcept { return x; }
[[nodiscard]]
inline std::string cat(char x) noexcept { return std::string(1,x); }
[[nodiscard]]
inline std::string cat(std::string_view x) noexcept { return std::string(x); }

template <typename... T>
[[nodiscard]]
[[gnu::always_inline]]
inline auto cat(T... x) -> std::enable_if_t<
  (sizeof...(x) > 1) && (impl::is_sv<T> && ...),
  std::string
> {
  std::string s((x.size() + ...),{});
  char* p = s.data();
  ( ( memcpy(p, x.data(), x.size()), p += x.size() ), ...);
  return s;
}

template <typename... T>
[[nodiscard]]
[[gnu::always_inline]]
inline auto cat(const T&... x) -> std::enable_if_t<
  (sizeof...(x) > 1) && !(impl::is_sv<T> && ...),
  std::string
> {
  return cat(impl::to_string_view(x)...);
}

struct chars_less {
  using is_transparent = void;
  bool operator()(const char* a, const char* b) const noexcept {
    return strcmp(a,b) < 0;
  }
  template <typename T> requires std::is_class_v<T>
  bool operator()(const T& a, const char* b) const noexcept {
    return a < b;
  }
  template <typename T> requires std::is_class_v<T>
  bool operator()(const char* a, const T& b) const noexcept {
    return a < b;
  }
};

} // end namespace ivanp

#endif
