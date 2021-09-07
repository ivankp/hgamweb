#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <charconv>
#include <filesystem>

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

// #include "linalg.hh"
#include "wls.hh"
#include "json.hh"
#include "error.hh"

using std::cout;
using std::cerr;
using std::endl;
// using std::get;
using ivanp::cat;
// using ivanp::sq;

double lumi_data = 0, lumi_mc = 0;

static_assert(sizeof(float )==4);
static_assert(sizeof(double)==8);

struct reader {
  char* m; size_t n; bool is_mc;
  reader(const char* filename) {
    const int fd = ::open(filename,O_RDONLY);
    if (fd < 0) THROW_ERRNO("open",filename);

    struct stat sb;
    if (::fstat(fd,&sb) < 0) THROW_ERRNO("fstat",filename);

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
  ~reader() { ::free(m); }
  struct range { float *begin, *end; };
  range operator*() const {
    float* a = reinterpret_cast<float*>(m+16);
    return { a, a+n };
  }
};

constexpr unsigned overflow = -1;
unsigned find_bin(double x, double a, double b, unsigned n) noexcept {
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
  if (argc!=2) {
    cout << "usage: " << argv[0] << " json_config_string\n";
    return 1;
  }

  using clock = std::chrono::system_clock;
  using tpoint = std::chrono::time_point<clock>;
  const tpoint start = clock::now();

  ivanp::json card(argv[1]);

  const unsigned nS = (unsigned)card["S"]["ndiv"] * (160-105);
  double* histS = new double[nS];
  double mS[3] { };

  const unsigned nB = (unsigned)card["B"]["ndiv"] * (160-105);
  double* histB = new double[nB];
  const unsigned nb[3] {
    nB*(121-105)/(160-105),
    nB*(129-121)/(160-105),
    nB*(160-129)/(160-105)
  };

  const auto prefix = (const std::string&)card["var"]+'-';
  for (const auto& file : std::filesystem::directory_iterator(
    (const std::string&)card["set"]
  )) {
    if (!( file.is_regular_file()
        && file.path().filename().native().starts_with(prefix)
    )) continue;
    reader events(file.path().c_str());
    if (events.is_mc) {
      for (auto [event,end] = *events; event!=end; event+=4) {
        double myy    = event[0] - 125.;
        double var    = event[1];
        double truth  = event[2];
        double weight = event[3];

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
        double myy    = event[0] - 125.;
        double var    = event[1];

        unsigned i = find_bin(myy,-20,35,nB);
        if (i!=overflow) ++histB[i];
      }
    }
  }

  mS[1] /= mS[0];
  mS[2] = mS[2]/mS[0] - mS[1]*mS[1];

  const unsigned nc = (unsigned)card["B"]["deg"]+1;
  double* csB = new double[nc];
  { const unsigned nB2 = nb[0]+nb[2];
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
    ivanp::wls(A, ys, us, nB2, nc, csB);
    delete[] us;
    delete[] ys;
    delete[] A;
  }

  std::stringstream out;
  out.precision(8);
  out << "{\"S\":{\"hist\":["
      << histS[0];
  for (unsigned i=1; i<nS; ++i)
    out << ',' << histS[i];
  out << "],\"mean\":" << mS[1]
      << ",\"stdev\":" << mS[2]
      << "},\"B\":{\"hist\":["
      << histB[0];
  { unsigned i = 1;
    for (; i<nb[0]; ++i)
      out << ',' << histB[i];
    i += nb[1];
    for (; i<nB; ++i)
      out << ',' << histB[i];
  }
  out << "],\"cs\":["
      << csB[0];
  for (unsigned i=1; i<nc; ++i)
    out << ',' << csB[i];
  out << "]},\"time\":";

  const auto time = std::chrono::duration<double>(clock::now()-start).count();
  out << time << '}';

  puts(std::move(out).str().c_str());

  delete[] csB;
  delete[] histB;
  delete[] histS;
}
