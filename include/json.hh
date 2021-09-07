#include <variant>
#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <iostream>

namespace ivanp {

struct json;

template <typename T>
struct json_converter {
  static void from_json(json&,T&);
  static void to_json(json&,T&);
};

struct json {
  struct  null_t { };
  using   bool_t = bool;
  using  float_t = double;
  using    int_t = int64_t;
  using string_t = std::string;
  using  array_t = std::vector<json>;
  using object_t = std::map<std::string,json,std::less<>>;

  using node_t = std::variant<
    null_t,
    bool_t,
    int_t,
    float_t,
    string_t,
    array_t,
    object_t
  >;
  node_t node;

  template <typename T>
  static constexpr bool is_null_v = std::is_same_v<T,null_t>;
  template <typename T>
  static constexpr bool is_bool_v = std::is_same_v<T,bool_t>;
  template <typename T>
  static constexpr bool is_int_v = std::is_same_v<T,int_t>;
  template <typename T>
  static constexpr bool is_float_v = std::is_same_v<T,float_t>;
  template <typename T>
  static constexpr bool is_string_v = std::is_same_v<T,string_t>;
  template <typename T>
  static constexpr bool is_array_v = std::is_same_v<T,array_t>;
  template <typename T>
  static constexpr bool is_object_v = std::is_same_v<T,object_t>;
  template <typename T>
  static constexpr bool is_value_v =
    is_null_v  <T> ||
    is_bool_v  <T> ||
    is_int_v   <T> ||
    is_float_v <T> ||
    is_string_v<T> ||
    is_array_v <T> ||
    is_object_v<T>;

  struct error: std::runtime_error {
    using std::runtime_error::runtime_error;
  };

  template <typename T>
  requires(std::is_same_v<T,node_t>)
  json(T&& node): node(std::move(node)) { }

  json() = default;
  json(std::string_view);

  template <typename F>
  decltype(auto) operator()(F&& f) {
    return std::visit(f,node);
  }
  template <typename F>
  decltype(auto) operator()(F&& f) const {
    return std::visit(f,node);
  }

  json& operator[](size_t key) {
    return std::get<array_t>(node).at(key);
  }
  const json& operator[](size_t key) const {
    return std::get<array_t>(node).at(key);
  }
  json& operator[](std::string_view key) {
    auto& obj = std::get<object_t>(node);
    auto it = obj.find(key);
    if (it!=obj.end()) return it->second;
    throw error("non-existent key");
  }
  const json& operator[](std::string_view key) const {
    auto& obj = std::get<object_t>(node);
    auto it = obj.find(key);
    if (it!=obj.end()) return it->second;
    throw error("non-existent key");
  }

  const char* type_name() const;

  constexpr bool is_null  () const noexcept { return node.index()==0; }
  constexpr bool is_bool  () const noexcept { return node.index()==1; }
  constexpr bool is_int   () const noexcept { return node.index()==2; }
  constexpr bool is_float () const noexcept { return node.index()==3; }
  constexpr bool is_string() const noexcept { return node.index()==4; }
  constexpr bool is_array () const noexcept { return node.index()==5; }
  constexpr bool is_object() const noexcept { return node.index()==6; }

  // get ------------------------------------------------------------
  template <typename T>
  requires(
    is_value_v<std::remove_cvref_t<T>> &&
    std::is_reference_v<T> &&
    !std::is_const_v<std::remove_reference_t<T>>
  )
  constexpr T get() {
    if constexpr (std::is_rvalue_reference_v<T>)
      return std::get<std::remove_cvref_t<T>>(std::move(node));
    else
      return std::get<std::remove_cvref_t<T>>(node);
  }
  template <typename T>
  requires(
    is_value_v<std::remove_cvref_t<T>> &&
    std::is_reference_v<T> &&
    std::is_const_v<std::remove_reference_t<T>>
  )
  constexpr T get() const {
    return std::get<std::remove_cvref_t<T>>(node);
  }

  template <typename T>
  requires(!std::is_reference_v<T>)
  constexpr T get() const {
    if constexpr (std::is_arithmetic_v<T>) {
      return std::visit([&]<typename X>(const X& x) -> T {
        if constexpr (std::is_arithmetic_v<X>) {
          return x;
        } else if constexpr (is_bool_v<T>) {
          if constexpr (json::is_null_v<X>)
            return false;
          else
            return !x.empty();
        } else
          throw error("cannot cast non-arithmetic json value to arithmetic type");
      },node);
    } else if constexpr (is_value_v<std::remove_cv_t<T>>) {
      return get<const T&>();
    } else {
      using type = std::remove_cvref_t<T>;
      type x;
      json_converter<type>::from_json(*this,x);
      return x;
    }
  }

  template <typename T>
  void get(T& x) const {
    if constexpr (std::is_arithmetic_v<T>)
      x = get<T>();
    else
      json_converter<std::remove_cv_t<T>>::from_json(*this,x);
  }

  // cast operators -------------------------------------------------
  template <typename T>
  operator T() const { return get<T>(); }
  // ----------------------------------------------------------------
};

std::ostream& operator<<(std::ostream& o, const json& j);

}
