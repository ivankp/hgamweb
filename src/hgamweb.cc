#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <charconv>

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "error.hh"

using std::cout;
using std::cerr;
using std::endl;
// using std::get;
using ivanp::cat;

bool is_mc;
double lumi_data = 0, lumi_mc = 0;
double myy, var, truth, weight;

static_assert(sizeof(float )==4);
static_assert(sizeof(double)==8);

struct file {
  char* m; size_t n;
  file(const char* filename) {
    const int fd = ::open(filename,O_RDONLY);
    if (fd < 0) THROW_ERRNO("open",filename);

    struct stat sb;
    if (::fstat(fd,&sb) < 0) THROW_ERRNO("fstat",filename);
    // if (!S_ISREG(sb.st_mode)) THROW(filename,": not a regular file");

    const size_t len = sb.st_size;
    if (len < 16) THROW(filename,": file is too short");
    m = reinterpret_cast<char*>(::malloc(len));
    const auto nread = ::read(fd,m,len);
    if (nread < 0) THROW_ERRNO("read",filename);
    if (size_t(nread) < len) THROW(filename,": couldn't read whole file");

    float lumi_f;
    memcpy(&lumi_f,m+8,sizeof(lumi_f));
    if (*m == 'd') { // data
      is_mc = false;
      lumi_data += lumi_f;
    } else if (*m == 'm') { // MC
      is_mc = true;
      lumi_mc += lumi_f;
    } else THROW(filename,": unexpected first byte");

    uint32_t nevents;
    memcpy(&nevents,m+12,sizeof(nevents));
    n = nevents;
    n *= (is_mc ? 4 : 2);
    if (n*4 != len-16)
      THROW(filename,": file length doesn't match expected number of events");

    ::close(fd);
  }
  ~file() { ::free(m); }
  struct range { float *begin, *end; };
  range operator*() const {
    float* a = reinterpret_cast<float*>(m+16);
    return { a, a+n };
  }
};

std::vector<double> getnums(std::string_view s) {
  std::vector<double> v;
  double x;
  const char *p = s.data(), *end = p+s.size();
  while (p!=end) {
    char c = *p;
    if (c==' ' || c=='\t' || c==',' || c==';' || c=='\n' || c=='\r') {
      ++p;
      continue;
    }
    auto [ptr,ec] = std::from_chars(p,end,x);
    if (ec == std::errc{}) {
      v.push_back(x);
      p = ptr;
    } else THROW("invalid string: ",s);
  }
  return v;
}

int main(int argc, char* argv[]) {
  if (argc<4) {
    cout << "usage: " << argv[0] <<
      " bin_edges S_params events.dat ...\n";
    return 1;
  }

  using clock = std::chrono::system_clock;
  using tpoint = std::chrono::time_point<clock>;
  const tpoint start = clock::now();

  auto edges = getnums(argv[1]);
  if (edges.size()<2) THROW("need at least two edges");
  std::sort(edges.begin(),edges.end());
  for (auto x : edges) TEST(x)
  auto S = getnums(argv[2]);
  if (S.size()!=6) THROW("DSCB needs 6 parameters");
  for (auto x : S) TEST(x)

  for (int argi=3; argi<argc; ++argi) {
    file events(argv[argi]);
    if (is_mc) {
      for (auto [event,end] = *events; event!=end; event+=4) {
        myy    = event[0];
        var    = event[1];
        truth  = event[2];
        weight = event[3];
      }
    } else {
      for (auto [event,end] = *events; event!=end; event+=2) {
        myy    = event[0];
        var    = event[1];
      }
    }
  }

  const auto time = std::chrono::duration<double>(clock::now()-start).count();
  TEST(time)
}
