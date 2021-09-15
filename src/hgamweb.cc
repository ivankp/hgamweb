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
#include <type_traits>

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <gsl/gsl_multimin.h>

#include "linalg.hh"
#include "wls.hh"
#include "json.hh"
#include "error.hh"

using std::cout;
using std::cerr;
using std::endl;
using ivanp::cat;
using ivanp::sq;

static_assert(sizeof(float )==4);
static_assert(sizeof(double)==8);

struct reader {
  char* m; size_t n; float lumi; bool is_mc;
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

    memcpy(&lumi,m+8,sizeof(lumi));
    if (*m == 'd') { // data
      is_mc = false;
    } else if (*m == 'm') { // MC
      is_mc = true;
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
  ++n;
  auto i = std::upper_bound(e,e+n,x) - e;
  return (i==0 || i==n) ? overflow : i-1;
}
double center(unsigned i, double a, double b, unsigned n) noexcept {
  return a + (i*2+1)*(b-a)/(n*2);
}

template <typename T, size_t N>
struct pool_array {
  T* const m[N];
  ~pool_array() { free(m[0]); }
  template <size_t I>
  friend T* get(const pool_array& x) { return x.m[I]; }
};
namespace std {
  template <typename T, size_t N>
  struct tuple_size<pool_array<T,N>>: std::integral_constant<size_t,N> { };

  template <size_t I, typename T, size_t N>
  struct tuple_element<I,pool_array<T,N>> { using type = T*; };
}
template <typename T, bool zero>
auto pool(auto... n) -> pool_array<T,sizeof...(n)> {
  // T* m = new T[(n + ...)];
  T* m;
  if constexpr (zero)
    m = reinterpret_cast<T*>(calloc((n + ...),sizeof(T)));
  else
    m = reinterpret_cast<T*>(malloc((n + ...)*sizeof(T)));
  return {{ ( m+=n, m-n ) ... }};
}

template <typename S, typename F>
struct stream_decorator {
  S s;
  F f;
  template <typename T>
  friend stream_decorator& operator<<(stream_decorator& d, T&& x) {
    d.f(d.s,std::forward<T>(x));
    return d;
  }
};

template <typename F> // n ≥ 1
double simpson(double a, double b, unsigned n, F&& f) {
  double h = (b-a)/n, h2 = h/2, sum1 = 0, sum2 = f(h2);
  for (unsigned i=1; i<n; ++i) {
    double x = a + h*i;
    sum1 += f(x);
    sum2 += f(x+h2);
  }
  return h/6 * (f(a) + f(b) + sum1*2 + sum2*4);
}

double poly(double x, const double* c, unsigned n) {
  const double *end = c+n;
  double f = *c, y = 1;
  while (++c < end) f += *c*(y*=x);
  return f;
}

struct tail_logL_params {
  double μ, σ;
  const double *hist;
  unsigned n;
};
double tail_logL(const gsl_vector* point, void* params) {
  auto [ μ, σ, ws, nw ] = *reinterpret_cast<const tail_logL_params*>(params);
  const double
    α = gsl_vector_get(point,0),
    p = gsl_vector_get(point,1),
    eα = std::exp(-0.5*α*α),
    pα = p/α;

  auto f = [=](double t) -> double {
    return t > α
    ? eα * std::pow( (pα-α+t)/pα, -p )
    : std::exp(-0.5*t*t);
  };

  // TODO: fit left tail too
  ws += nw*(125-105)/(160-105);
  const double xmax = (160-125), tmax = (xmax-μ)/σ;

  double wtot = 0, Itot = 0, logL = 0;

  const unsigned
    ni = nw*(160-125)/(160-105),
    nj = std::max(4*(160-105)/nw,1u); // at least 4 integration steps per GeV
  const double h=tmax/(ni*nj), h2=h/2;
  double sum1=0, sum2=f(h2), fa=f(0), fb;

  for (unsigned i=0; i<ni; ++i) { // loop over bins
    for (unsigned j=1; j<nj; ++j) { // loop over integration segments
      const double t = h*(j+nj*i);
      sum1 += f(t);
      sum2 += f(t+h2);
    }
    fb = f(tmax);
    const double I = h/6 * (fa + fb + sum1*2 + sum2*4); // simpson
    const double w = ws[i];
    Itot += I;
    logL -= w * std::log(I);
    wtot += w;
    fa = fb;
  }

  return logL + wtot * std::log(Itot);
}
void fit_tails(
  const double* mS, const double* histS,
  unsigned nbins, unsigned nm, unsigned nS
) {
  // https://www.gnu.org/software/gsl/doc/html/multimin.html
  gsl_vector *point = gsl_vector_alloc(2),
             *step  = gsl_vector_alloc(2);

  tail_logL_params params;
  params.n = nS;

  gsl_multimin_function minex_func;
  minex_func.n = 2;
  minex_func.f = tail_logL;
  minex_func.params = &params;

  gsl_multimin_fminimizer *minimizer = gsl_multimin_fminimizer_alloc(
    gsl_multimin_fminimizer_nmsimplex2, 2);
  gsl_multimin_fminimizer_set(minimizer, &minex_func, point, step);

  for (unsigned b=0; b<nbins; ++b) { // loop over var bins
    gsl_vector_set(point, 0, 1.5); gsl_vector_set(point, 1, 5.0);
    gsl_vector_set(step , 0, 0.2); gsl_vector_set(step , 1, 1.0);

    const double* m = mS+nm*b;
    params.μ = m[1];
    params.σ = m[2];
    params.hist = histS+nS*b;

    int status;
    double size;
    for (unsigned iter=0;;) {
      status = gsl_multimin_fminimizer_iterate(minimizer);
      if (status) break;
      size = gsl_multimin_fminimizer_size(minimizer);
      status = gsl_multimin_test_size(size,1e-2);
      printf("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n",
        iter,
        gsl_vector_get(point,0),
        gsl_vector_get(point,1),
        minimizer->fval, size);
      if (status != GSL_CONTINUE || ++iter > 100) break;
    }
  }

  gsl_multimin_fminimizer_free(minimizer);
  gsl_vector_free(step);
  gsl_vector_free(point);
}

int main(int argc, char* argv[]) {
  if (argc!=2) {
    cout << "usage: " << argv[0] << " json_config_string\n";
    return 1;
  }

  // initialize timer to measure execution time ---------------------
  using clock = std::chrono::high_resolution_clock;
  const auto start = clock::now();

  // read configuration string --------------------------------------
  cout << argv[1] << endl;
  ivanp::json card(argv[1]);

  // get numbers of values ------------------------------------------
  const unsigned
    nm = 3,
    nbins = card["edges"].size()-1,
    nS = (unsigned)card["Sdiv"] * (160-105),
    nB = (unsigned)card["Bdiv"] * (160-105),
    nV = card["nV"],
    nc = (unsigned)card["Bdeg"]+1,
    ncS = 6;

  const unsigned nb[] {
    nB*(121-105)/(160-105),
    nB*(129-121)/(160-105),
    nB*(160-129)/(160-105)
  };

  // allocate arrays ------------------------------------------------
  const auto [
    edges  , histS   , mS       , w2S  , histB   , histVfine, histV, cB      ,
    bkg    , rewbkg  , cS
  ] = pool<double,true>(
    nbins+1, nS*nbins, nm*nbins , nbins, nB*nbins, nV       , nbins, nc*nbins,
    nbins  , nbins   , ncS*nbins
  );

  // parse edges ----------------------------------------------------
  { const ivanp::json::array_t& xs = card["edges"];
    for (unsigned i=0; i<=nbins; ++i) {
      const auto& x = xs[i];
      double v;
      if (x.is_string())
        v = std::stod(x.get<ivanp::json::string_t>());
      else
        v = x;
      edges[i] = v;
    }
  }

  // read MC event files --------------------------------------------
  { reader events(cat(card["H"],'/',card["V"],"-mc").c_str());
    const double aV = edges[0], bV = edges[nbins];
    for (auto [event,end] = *events; event!=end; event+=4) {
      double myy    = event[0] - 125;
      double var    = event[1];
      double truth  = event[2];
      double weight = event[3];

      unsigned b = find_bin(var,edges,nbins);
      if (b==overflow) continue;
      histV[b] += weight;

      histVfine[find_bin(var,aV,bV,nV)] += weight;

      w2S[b] += weight*weight;

      unsigned iS = find_bin(myy,-20,35,nS);
      if (iS!=overflow) histS[b*nS+iS] += weight;

      for (unsigned i=0; ; ++i) {
        mS[b*3+i] += weight;
        if (i == 2) break;
        weight *= myy;
      }
    }
  }

  // read data event files ------------------------------------------
  float lumi;
  { reader events(cat(card["H"],'/',card["V"],"-data").c_str());
    lumi = events.lumi;
    for (auto [event,end] = *events; event!=end; event+=2) {
      double myy = event[0] - 125;
      double var = event[1];

      unsigned b = find_bin(var,edges,nbins);
      if (b==overflow) continue;

      unsigned iB = find_bin(myy,-20,35,nB);
      if (iB!=overflow) ++histB[b*nB+iB];
    }
  }

  // background fit (weighted least squares) ------------------------
  const bool fitexp = card["Bf"].get<std::string>() == "exppoly";
  { const unsigned nB2 = nb[0]+nb[2];
    const auto [ ys, us, A ] = pool<double,false>( nB2, nB2, nB2*nc );
    std::fill(A, A+nB2, 1.);
    const std::string& Bf = card["Bf"];
    if (!fitexp) { // fit polynomial
      for (unsigned b=0; b<nbins; ++b) {
        for (unsigned i=0, j=0; i<nB2; ++i, ++j) {
          if (j == nb[0]) j += nb[1];
          double y = histB[b*nB+j],
                *a = A + i,
                 x = center(j,-20,35,nB),
                 f = 1.;
          ys[i] = y;
          us[i] = y>0 ? y : 1.;
          for (unsigned k=1; k<nc; ++k)
            *(a+=nB2) = (f*=x);
        }
        ivanp::wls(A, ys, us, nB2, nc, cB+b*nc);

        // integration
        const double *m = mS + b*nm, *c = cB + b*nc;
        bkg[b] = simpson(-4,4,10,[=](double x){
          return poly(x,c,nc);
        });
        rewbkg[b] = simpson(-20,35,50,[=](double x){
          return std::exp(-0.5*sq((x-m[1])/m[2])) * poly(x,c,nc);
        });
      }
    } else { // fit exp(polynomial)
      for (unsigned b=0; b<nbins; ++b) {
        for (unsigned i=0, j=0; i<nB2; ++i, ++j) {
          if (j == nb[0]) j += nb[1];
          double y = histB[b*nB+j],
                *a = A + i,
                 x = center(j,-20,35,nB),
                 f = 1.;
          if (y>0) {
            ys[i] = std::log(y);
            us[i] = 1./y;
          } else {
            ys[i] = 0.;
            us[i] = 1.;
          }
          for (unsigned k=1; k<nc; ++k)
            *(a+=nB2) = (f*=x);
        }
        double *m = mS + b*nm, *c = cB + b*nc;

        ivanp::wls(A, ys, us, nB2, nc, c);

        // integration
        bkg[b] = simpson(-4,4,10,[=](double x){
          return std::exp(poly(x,c,nc));
        });
        rewbkg[b] = simpson(-20,35,50,[=](double x){
          return std::exp(-0.5*sq((x-m[1])/m[2])) * std::exp(poly(x,c,nc));
        });
      }
    }
  }

  // signal distribution moments ------------------------------------
  for (double* m = mS+nm*nbins; m!=mS; ) {
    m -= nm;
    m[1] /= m[0];
    m[2] = std::sqrt(m[2]/m[0] - m[1]*m[1]);
  }

  // fit signal tails -----------------------------------------------
  fit_tails(mS, histS, nbins, nm, nS);

  // JSON output ----------------------------------------------------
  stream_decorator out {
    std::stringstream{},
    [](auto& s, const auto& x) {
      if constexpr (std::is_floating_point_v<std::remove_cvref_t<decltype(x)>>) {
        const bool ok = std::isfinite(x);
        if (!ok) s << '\"';
        s << x;
        if (!ok) s << '\"';
      } else {
        s << x;
      }
    }
  };
  out.s.precision(8);

  { double* x; unsigned i;
    out << "{\"lumi\":" << lumi
        << ",\"edges\":["
        << edges[0];
    for (i=1; i<=nbins; ++i)
      out << ',' << edges[i];
    out << "],\"bins\":[";
    for (unsigned b=0; b<nbins; ++b) {
      if (b) out << ',';
      out << "{\"S\":{\"hist\":[";

      x = histS + b*nS;
      out << x[0];
      for (i=1; i<nS; ++i)
        out << ',' << x[i];

      x = mS + b*nm;
      out << "],\"mean\":" << x[1]
          << ",\"stdev\":" << x[2]
          << ",\"rootw2\":" << std::sqrt(w2S[b])
          << "},\"B\":{\"hist\":[";

      x = histB + b*nB;
      out << x[0];
      for (i=1; i<nb[0]; ++i)
        out << ',' << x[i];
      for (i+=nb[1]; i<nB; ++i)
        out << ',' << x[i];
      out << "],\"cs\":[";

      x = cB + b*nc;
      out << x[0];
      for (i=1; i<nc; ++i)
        out << ',' << x[i];
      out << "],\"bkg\":[" << bkg[b] << ',' << rewbkg[b] << "]}}";
    }
    out << "],\"V\":[["
        << histVfine[0];
    for (i=1; i<nV; ++i)
      out << ',' << histVfine[i];
    out << "],["
        << histV[0];
    for (i=1; i<nbins; ++i)
      out << ',' << histV[i];
    out << "]],\"time\":";
  }

  const auto time = std::chrono::duration<double>(clock::now()-start).count();
  out << time << '}';

  puts(std::move(out.s).str().c_str());
}
