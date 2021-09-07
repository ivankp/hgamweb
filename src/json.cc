#include "json.hh"
#include <cstring>
// #include <charconv>

#define STR1(x) #x
#define STR(x) STR1(x)

#define TEST(var) std::cout << \
  "\033[33m" STR(__LINE__) ": " \
  "\033[36m" #var ":\033[0m " << (var) << std::endl;

namespace ivanp {
namespace {

bool is_blank(char c) noexcept {
  switch (c) {
    case ' ':
    case '\t':
    case '\r':
    case '\n':
      return true;
    default:
      return false;
  }
}
bool is_end(char c) noexcept {
  switch (c) {
    case ' ':
    case '\t':
    case '\r':
    case '\n':
    case ',':
    case ']':
    case '}':
    case '\0':
      return true;
    default:
      return false;
  }
}

const char* skip_blank(const char* s, const char* end) {
  for (; s<end; ++s)
    if (!is_blank(*s)) break;
  if (s>=end)
    throw json::error("json source ended unexpectidly");
  return s;
}

struct parser_ret_t {
  json::node_t node;
  const char* ptr;
};

parser_ret_t parse_number(const char* s, const char* end) {
  /*
  { json::int_t x;
    auto [ptr,ec] = std::from_chars(s,end,x);
    if (ptr!=s) {
      const char c = *ptr;
      if (!(c=='.' || c=='e' || c=='E'))
        return {x,ptr};
    }
    if (ec==std::errc::invalid_argument)
      throw json::error("invalid json number");
  }
  { json::float_t x;
    auto [ptr,ec] = std::from_chars(s,end,x);
    if (ptr==s)
      throw json::error("invalid json number");
    return {x,ptr};
  }
  */
  char* num_end;
  { const json::int_t x = std::strtoll(s,&num_end,10);
    if (num_end==s)
      throw json::error("invalid json number");
    if (end<num_end) // TODO: avoid this
      throw json::error("*** number continued past the end");
    const char c = *num_end;
    if (!(c=='.' || c=='e' || c=='E'))
      return {x,num_end};
  }
  { const json::float_t x = std::strtod(s,&num_end);
    if (num_end==s)
      throw json::error("invalid json number");
    if (end<num_end) // TODO: avoid this
      throw json::error("*** number continued past the end");
    return {x,num_end};
  }
}

parser_ret_t parse_string(const char* s, const char* end) {
  unsigned esc = 0;
  json::string_t str;
  for (char c; s!=end; ++s) {
    c = *s;
    if (c=='\\') ++esc;
    else {
      while (esc > 1) {
        str += '\\';
        esc -= 2;
      }
      if (esc) {
        // TODO: interpret escaped chars
        str += '\\';
        str += c;
        esc = 0;
      } else if (c=='"') {
        return {std::move(str),s+1};
      } else {
        str += c;
      }
    }
  }
  throw json::error("unterminated json string");
}

void match(const char* seq, const char* s, const char* end) {
  const size_t len = end - s;
  const size_t seq_len = strlen(++seq);
  if (len<seq_len ||
      strncmp(s,seq,std::min(seq_len,len)) ||
      ( len>seq_len && !is_end(*(s+seq_len)) )
  ) throw json::error("expected "+std::string(seq-1));
}

parser_ret_t parse_node(const char* s, const char* end) {
  switch (*s++) {
    case '"':
      return parse_string(s,end);
    case '[': {
      json::array_t arr;
      s = skip_blank(s,end);
      for (char p='[';;) {
        char c = *s;
        if (c==']') {
          if (p==',') throw json::error("unexpected ]");
          ++s;
          break;
        } else if (c==',') {
          if (p) throw json::error("unexpected ,");
          s = skip_blank(s+1,end);
          p = ',';
        } else {
          if (!p) throw json::error("expected , or ]");
          auto [val,s2] = parse_node(s,end);
          s = skip_blank(s2,end);
          arr.emplace_back(std::move(val));
          p = '\0';
        }
      }
      return {std::move(arr),s};
    } break;
    case '{': {
      json::object_t obj;
      s = skip_blank(s,end);
      for (bool not_last=(*s!='}'); not_last; ++s) {
        s = skip_blank(s,end);
        if (*s!='"')
          throw json::error("json key must be a string");
        auto [key,s2] = parse_string(s+1,end);
        s = skip_blank(s2,end);
        if (*s!=':')
          throw json::error("expected : after a key");
        s = skip_blank(s+1,end);
        auto [val,s3] = parse_node(s,end);
        s = skip_blank(s3,end);
        switch (*s) {
          case '}':
            not_last = false;
          case ',':
            if (!obj.try_emplace(
              std::get<json::string_t>(std::move(key)),
              std::move(val)
            ).second)
              throw json::error("failed to insert key "
                + std::get<json::string_t>(key));
            break;
          default:
            throw json::error("expected , or }");
        }
      }
      return {std::move(obj),s};
    } break;
    case '0':
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
    case '9':
    case '-':
      return parse_number(s-1,end);
    case 'n':
      match("null",s,end);
      return {{},s+3};
    case 't':
      match("true",s,end);
      return {true,s+3};
    case 'f':
      match("false",s,end);
      return {false,s+4};
    default :
      throw json::error(std::string("invalid json token ")+*(s-1));
  }
}

parser_ret_t parse_first_node(const char* s, const char* end) {
  return parse_node(skip_blank(s,end),end);
}
parser_ret_t parse_first_node(const char* s, size_t n) {
  return parse_first_node(s,s+n);
}

} // end namespace

json::json(std::string_view src)
: node(parse_first_node(src.data(),src.size()).node)
{ }

const char* json::type_name() const {
  switch (node.index()) {
    case 0: return "null";
    case 1: return "bool";
    case 2: return "int";
    case 3: return "float";
    case 4: return "string";
    case 5: return "array";
    case 6: return "object";
    default: return nullptr;
  }
}

std::ostream& operator<<(std::ostream& o, const json& j) {
  using ivanp::json;
  std::visit([&]<typename T>(const T& x){
    if constexpr (json::is_null_v<T>)
      o << "null";
    else if constexpr (json::is_bool_v<T>)
      o << (x ? "true" : "false");
    else if constexpr (json::is_int_v<T> || json::is_float_v<T>)
      o << x;
    else if constexpr (json::is_string_v<T>)
      o << '"' << x << '"';
    else if constexpr (json::is_array_v<T>) {
      o << '[';
      bool first = true;
      for (const auto& val : x) {
        if (first) first = false;
        else o << ',';
        o << val;
      }
      o << ']';
    } else if constexpr (json::is_object_v<T>) {
      o << '{';
      bool first = true;
      for (const auto& [key,val] : x) {
        if (first) first = false;
        else o << ',';
        o << '"' << key << "\":" << val;
      }
      o << '}';
    }
  },j.node);
  return o;
}

} // end namespace ivanp
