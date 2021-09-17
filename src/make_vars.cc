#include <iostream>
#include <type_traits>
#include <limits>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cerrno>

#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <ctre.hpp>

#include <TFile.h>
#include <TKey.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
// #include <TLorentzVector.h>

#include "tcnt.hh"
#include "string.hh"

using std::cout;
using std::cerr;
using std::endl;
using ivanp::cat;

static_assert(sizeof(float )==4);
static_assert(sizeof(double)==8);
static_assert(std::is_same_v<float,Float_t>);

constexpr float fnan = std::numeric_limits<float>::quiet_NaN();

bool is_mc;
double mc_factor;
float lumi=0, weight, weight_j;
TTreeReader* reader_ptr;

class output_file {
  int fd = -1;
  unsigned nevents = 0;

  static constexpr unsigned buf_size = (1<<20)/4;
  float *buf, *p;

public:
  output_file(const char* filename, char h1)
  : buf(new float[buf_size]), p(buf)
  {
    fd = ::open(filename, O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) throw std::runtime_error(cat(
      "failed to open output file \"",filename,"\": ",std::strerror(errno) ));

    char head[16];
    std::fill(head,head+16,' ');
    head[0] = h1;
    ::write(fd,head,16);
  }
  ~output_file() {
    if (p!=buf) flush();
    delete[] buf;

    char head[8];
    ::lseek(fd,0,SEEK_SET);
    ::read(fd,head,1);

    float x = head[0]=='d' ? lumi : 1.f;
    memcpy(head,&x,4);
    memcpy(head+4,&nevents,4);

    ::lseek(fd,8,SEEK_SET);
    ::write(fd,head,8);

    ::close(fd);
  }

  void flush() {
    ::write(fd,buf,(p-buf)*4);
    p = buf;
  }

  template <typename... T>
  void write(T... x) requires (
    std::is_same_v<T,float> && ... && (sizeof...(x)==2 || sizeof...(x)==4)
  ) {
    if (p > buf+(buf_size-sizeof...(x))) flush();
    ( (*p++ = x), ... );
    ++nevents;
  }
};

struct vardef_base {
  using V0 = TTreeReaderValue<double>;
  std::aligned_storage_t<sizeof(V0),alignof(V0)> m;
  std::string name;

  template <typename S>
  vardef_base(S&& name): name(std::forward<S>(name)) { }

  virtual void init() = 0;
  virtual void free() = 0;
};
std::vector<vardef_base*> vars;

template <typename T, bool mc=false, T(*f)(T)=nullptr>
struct vardef: vardef_base {
  using V = TTreeReaderValue<T>;
  static_assert( sizeof(V) <= sizeof(V0) );
  static constexpr bool has_f = (f != nullptr);

  template <typename S>
  vardef(S&& name): vardef_base(std::forward<S>(name))
  { vars.push_back(this); }

  void init() {
    if (!mc || is_mc) new(&m) V(*reader_ptr,name.c_str());
  }
  void free() {
    if (!mc || is_mc) reinterpret_cast<V*>(&m)->~V();
  }

  decltype(auto) operator*() {
    if constexpr (has_f) return (*f)(**reinterpret_cast<V*>(&m));
    else return **reinterpret_cast<V*>(&m);
  }
};

template <typename T> T f_1e3(T x) { return x*1e-3; }
template <typename T> T f_abs(T x) { return std::abs(x); }

/*
template <typename T>
struct vec { T pt, eta, phi, m; };

template <typename T>
std::vector<TLorentzVector> select_photons(vec<T>& v, double myy) {
  const auto &pt = *v.pt, &eta = *v.eta, &phi = *v.phi, &m = *v.m;
  const unsigned n = pt.size();
  if (n < 2) throw std::runtime_error("fewer than 2 photons");
  std::vector<TLorentzVector> photons(n);
  for (unsigned i=0; i<n; ++i)
    photons[i].SetPtEtaPhiM( pt[i], eta[i], phi[i], m[i] );
  if (n > 2) {
    unsigned a=0, b=1;
    double dmin = std::numeric_limits<double>::infinity();
    for (unsigned i=0; i<n; ++i) {
      for (unsigned j=i+1; j<n; ++j) {
        const double d = std::abs( (photons[i] + photons[j]).M() - myy );
        if (d < dmin) {
          dmin = d;
          a = i;
          b = j;
        }
      }
    }
    photons = { photons[a], photons[b] };
  }
  return photons;
}

template <typename T>
std::vector<TLorentzVector> select_jets(vec<T>& v) {
  const auto &pt = *v.pt, &eta = *v.eta, &phi = *v.phi, &m = *v.m;
  const unsigned n = pt.size();
  std::vector<TLorentzVector> jets(n);
  for (unsigned i=0; i<n; ++i)
    jets[i].SetPtEtaPhiM( pt[i], eta[i], phi[i], m[i] );
  std::sort(jets.begin(),jets.end(),[](const auto& a, const auto& b){
    return a.Pt() > b.Pt();
  });
  return jets;
}
*/

int main(int argc, char* argv[]) {
  if (argc<3) {
    cout << "usage: " << argv[0] << " output_dir input.root ...\n";
    return 1;
  }

#define VAR_PREF "HGamEventInfoAuxDyn."
#define VAR_PREF_TRUTH "HGamTruthEventInfoAuxDyn."

#define GET_MACRO(_1,_2,_3,_4,NAME,...) NAME

#define VAR(...) GET_MACRO(__VA_ARGS__, VAR_4, VAR_3, VAR_2, VAR_1)(__VA_ARGS__)
#define VAR_4(NAME,F,TYPE,ACTUAL) \
  vardef<TYPE,0,F> v_##NAME(VAR_PREF ACTUAL); \
  vardef<TYPE,1,F> v_##NAME##_truth(VAR_PREF_TRUTH ACTUAL); \
  output_file f_##NAME##_data(cat(argv[1],'/',#NAME"-data").c_str(),'d'); \
  output_file f_##NAME##_mc  (cat(argv[1],'/',#NAME"-mc"  ).c_str(),'m');

#define VAR_3(NAME,F,TYPE) VAR_4(NAME,F,TYPE,#NAME)
#define VAR_2(NAME,F) VAR_3(NAME,F,Float_t)
#define VAR_1(NAME) VAR_2(NAME,nullptr)

#define VARJ(...) GET_MACRO(__VA_ARGS__, _, VARJ_3, VARJ_2, VARJ_1)(__VA_ARGS__)
#define VARJ_3(NAME,F,TYPE) VAR(NAME,F,TYPE,#NAME"_30")
#define VARJ_2(NAME,F) VAR(NAME,F,Float_t,#NAME"_30")
#define VARJ_1(NAME) VAR(NAME,nullptr,Float_t,#NAME"_30")

  vardef<Char_t>    v_isPassed   (VAR_PREF "isPassed");
  vardef<Float_t,1> v_cs_br_fe   (VAR_PREF "crossSectionBRfilterEff");
  vardef<Float_t,1> v_weight     (VAR_PREF "weight");
  vardef<Float_t,1> v_weight_jvt (VAR_PREF "weightJvt_30");

  VARJ(N_j,nullptr,Int_t)

  VAR(m_yy,f_1e3)
  VAR(pT_yy,f_1e3)  VAR(yAbs_yy)      VAR(cosTS_yy,f_abs)
  VAR(pTt_yy,f_1e3) VAR(Dy_y_y,f_abs)

  VARJ(HT,f_1e3)         VARJ(HTall,f_1e3)
  VARJ(pT_j1,f_1e3)      VARJ(pT_j2,f_1e3)      VARJ(pT_j3,f_1e3)
  VARJ(yAbs_j1)          VARJ(yAbs_j2)
  VARJ(Dphi_j_j,f_abs)   VAR(Dphi_j_j_signed,nullptr,Float_t,"Dphi_j_j_30_signed")
  VAR (pT_jj,f_1e3)      VARJ(m_jj,f_1e3)
  VARJ(Dy_j_j,f_abs)

  VARJ(sumTau_yyj,f_1e3) VARJ(maxTau_yyj,f_1e3)
  VARJ(Dphi_yy_jj,f_abs) VAR(cosTS_yyjj,f_abs)

  VARJ(pT_yyj,f_1e3)     VARJ(m_yyj,f_1e3)
  VARJ(pT_yyjj,f_1e3)    VAR (m_yyjj,f_1e3)

/*
#define PHOTONS_PREF "HGamPhotonsAuxDyn."
  vec<vardef<std::vector<float>>> v_photons {
    PHOTONS_PREF "pt",
    PHOTONS_PREF "eta",
    PHOTONS_PREF "phi",
    PHOTONS_PREF "m"
  };

#define JETS_PREF "HGamAntiKt4PFlowCustomVtxHggJetsAuxDyn."
  vec<vardef<std::vector<float>>> v_jets {
    JETS_PREF "pt",
    JETS_PREF "eta",
    JETS_PREF "phi",
    JETS_PREF "m"
  };
*/

  for (int argi=2; argi<argc; ++argi) {
    cout << argv[argi] << endl;
    TFile file(argv[argi]);
    if (file.IsZombie()) return 1;

    if (auto m = ctre::match<
      "(?:^.*/)?mc16([ade])\\."
    >(argv[argi])) {
      is_mc = true;
      for (auto* key : *file.GetListOfKeys()) {
        const std::string_view name = key->GetName();
        if (!(name.starts_with("CutFlow_") &&
              name.ends_with("_noDalitz_weighted"))) continue;
        TH1 *h = dynamic_cast<TH1*>(static_cast<TKey*>(key)->ReadObj());
        if (!h) {
          cerr << name << " is not a TH1\n";
          return 1;
        }
        cout << name << '\n';
        const double n_all = h->GetBinContent(3);
        cout << h->GetXaxis()->GetBinLabel(3) << " = " << n_all << '\n';
        mc_factor = 1e3/n_all;
        cout << "MC factor = " << mc_factor << endl;
        break;
      }
    } else if (auto m = ctre::match<
      R"((?:^.*/)?data(\d{2})_.*\.DS\d+_(\d+)ipb\.)"
    >(argv[argi])) {
    }

    TTreeReader reader("CollectionTree",&file);
    reader_ptr = &reader;
    for (auto* var : vars) var->init();

    for (ivanp::tcnt ent(reader.GetEntries(true)); reader.Next(); ++ent) {
      if (!*v_isPassed) continue;

      const float myy = *v_m_yy;
      if (myy < 105 || 160 <= myy) continue;

      if (is_mc) {
        double w = (*v_weight) * (*v_cs_br_fe) * mc_factor;
        weight = w;
        weight_j = w * (*v_weight_jvt);
      } else {
        // blind data in the signal region
        if (121 < myy && myy < 129) continue;
      }

#define WRITE(NAME) \
      if (!is_mc) f_##NAME##_data.write(myy,(float)*v_##NAME); \
      else f_##NAME##_mc.write(myy,(float)*v_##NAME,(float)*v_##NAME##_truth,weight);

      WRITE(pT_yy)
      WRITE(HT)
      WRITE(HTall)
      WRITE(yAbs_yy)
      WRITE(cosTS_yy)
      WRITE(pTt_yy)
      WRITE(Dy_y_y)

#undef WRITE
#define WRITE(NAME) \
      if (!is_mc) f_##NAME##_data.write(myy,(float)*v_##NAME); \
      else f_##NAME##_mc.write(myy,(float)*v_##NAME,(float)*v_##NAME##_truth,weight_j);

      WRITE(N_j)

      const auto nj = *v_N_j;
      if (nj < 1) continue; // 1 jet --------------------------------

      WRITE(pT_j1)
      WRITE(yAbs_j1)
      WRITE(sumTau_yyj)
      WRITE(maxTau_yyj)
      WRITE(pT_yyj)
      WRITE(m_yyj)

      /*
      const auto photons = select_photons(v_photons,myy);
      const auto jets = select_jets(v_jets);

      const TLorentzVector Hj = photons[0]+photons[1]+jets[0];
      const float m_yyj = Hj.M();
      if (!is_mc) f_m_yyj_data.write(myy,m_yyj);
      else f_m_yyj_mc.write(myy,m_yyj,fnan,weight_j);
      */

      if (nj < 2) continue; // 2 jet --------------------------------

      WRITE(pT_j2)
      WRITE(yAbs_j2)
      WRITE(cosTS_yyjj)
      WRITE(Dphi_j_j_signed)
      WRITE(Dphi_j_j)
      WRITE(Dy_j_j)
      WRITE(pT_jj)
      WRITE(m_jj)
      WRITE(pT_yyjj)
      WRITE(m_yyjj)

      if (nj < 3) continue; // 3 jet --------------------------------

      WRITE(pT_j3)

#undef WRITE
    } // end event loop

    for (auto* var : vars) var->free();
  } // end loop over input files
}
