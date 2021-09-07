#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <charconv>

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

// #include "linalg.hh"
#include "wls.hh"
#include "error.hh"

using std::cout;
using std::cerr;
using std::endl;
// using std::get;
// using ivanp::cat;
// using ivanp::sq;

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

constexpr unsigned overflow = -1;
unsigned find_bin(double x, double a, double b, unsigned n) noexcept {
  // if (x < a) return 0;
  // if (!(x < b)) return n + 1u;
  // return unsigned(n*(x-a)/(b-a)) + 1u;
  if (x < a || !(x < b)) return overflow;
  return n*(x-a)/(b-a);
}
unsigned find_bin(double x, double* e, unsigned n) noexcept {
  return std::upper_bound(e,e+n,x) - e;
}
double center(unsigned i, double a, double b, unsigned n) noexcept {
  return a + (i*2+1)*(b-a)/(n*2);
}

int main(int argc, char* argv[]) {
  if (argc<2) {
    cout << "usage: " << argv[0] << " files ...\n";
    return 1;
  }

  using clock = std::chrono::system_clock;
  using tpoint = std::chrono::time_point<clock>;
  const tpoint start = clock::now();

  // for (auto x : getnums(argv[1])) TEST(x)

  // auto edges = getnums(argv[2]);
  // if (edges.size()<2) THROW("need at least two edges");
  // std::sort(edges.begin(),edges.end());
  // for (auto x : edges) TEST(x)
  // auto S = getnums(argv[2]);
  // if (S.size()!=6) THROW("DSCB needs 6 parameters");
  // for (auto x : S) TEST(x)

  unsigned nS = 110;
  double* histS = new double[nS];

  unsigned nB = 110;
  double* histB = new double[nB];

  double mS[3] { };

  for (int argi=1; argi<argc; ++argi) {
    file events(argv[argi]);
    if (is_mc) {
      for (auto [event,end] = *events; event!=end; event+=4) {
        myy    = event[0] - 125.;
        var    = event[1];
        truth  = event[2];
        weight = event[3];

        unsigned i = find_bin(myy,-20,35,nS);
        if (i!=overflow) histS[i] += weight;

        double dm = weight;
        for (int i=0;;) {
          mS[i] += dm;
          if (++i > 2) break;
          dm *= myy;
        }
      }
    } else {
      for (auto [event,end] = *events; event!=end; event+=2) {
        myy    = event[0] - 125.;
        var    = event[1];

        unsigned i = find_bin(myy,-20,35,nB);
        if (i!=overflow) ++histB[i];
      }
    }
  }

  mS[1] /= mS[0];
  mS[2] = mS[2]/mS[0] - mS[1]*mS[1];
  TEST(mS[1])
  TEST(mS[2])

  double cs[4];
  const unsigned nc = std::size(cs);
  { const unsigned nb[3] {
      nB*(121-105)/(160-105),
      nB*(129-121)/(160-105),
      nB*(160-129)/(160-105)
    };
    const unsigned nB2 = nb[0]+nb[2];
    double* ys = new double[nB2];
    double* us = new double[nB2];
    double* A = new double[nB2*nc];
    for (unsigned i=0, j=0; i<nB2; ++i, ++j) {
      if (j == nb[0]) j += nb[1];
      double y = histB[j],
            *a = A + i,
             x = center(j,-20,35,nB),
             f = 1;
      ys[i] = y;
      us[i] = y ?: 1.;
      for (unsigned j=0;;) {
        *a = f;
        if (++j >= nc) break;
        a += nB2;
        f *= x;
      }
    }
    ivanp::wls(A, ys, us, nB2, nc, cs);
    delete[] us;
  }

  cout << "hist = [";
  for (unsigned i=0; i<nB; ++i) {
    if (i) cout << ',';
    cout << histB[i];
  }
  cout << ']' << endl;
  cout << "cs = [";
  for (unsigned i=0; i<nc; ++i) {
    if (i) cout << ',';
    cout << cs[i];
  }
  cout << ']' << endl;

  delete[] histB;
  delete[] histS;

  const auto time = std::chrono::duration<double>(clock::now()-start).count();
  TEST(time)
}
