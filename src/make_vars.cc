#include <iostream>
#include <string>
#include <string_view>
#include <type_traits>

#include <TFile.h>
#include <TKey.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
// #include <TLorentzVector.h>

using std::cout;
using std::cerr;
using std::endl;

static_assert(sizeof(float )==4);
static_assert(sizeof(double)==8);

bool is_mc;
double lumi, mc_factor;
float weight;
TTreeReader* reader_ptr;

struct vardef_base {
  using V0 = TTreeReaderValue<double>;
  std::aligned_storage_t<sizeof(V0),alignof(V0)> m;
  std::string name;
  // void* f;

  // template <typename S>
  // vardef_base(S&& name, void* f=nullptr)
  // : name(std::forward<S>(name)), f(f) { }

  template <typename S>
  vardef_base(S&& name): name(std::forward<S>(name)) { }

  virtual void init() = 0;
  virtual void free() = 0;
};
std::vector<vardef_base*> vars;
template <typename T, T(*f)(T) = nullptr>
struct vardef: vardef_base {
  using V = TTreeReaderValue<T>;

  // template <typename S>
  // vardef(S&& name, F* f=nullptr)
  // : vardef_base(std::forward<S>(name), reinterpret_cast<void*>(f))
  // { vars.push_back(this); }

  template <typename S>
  vardef(S&& name): vardef_base(std::forward<S>(name))
  { vars.push_back(this); }

  void init() {
    static_assert( sizeof(V) <= sizeof(V0) );
    new(&m) V(*reader_ptr,name.c_str());
  }
  void free() { reinterpret_cast<V*>(&m)->~V(); }

  T operator*() {
    if constexpr (f)
      return (*f)(**reinterpret_cast<V*>(&m));
    else
      return **reinterpret_cast<V*>(&m);
  }
};

template <typename T> T f_1e3(T x) { return x*1e-3; }
template <typename T> T f_abs(T x) { return std::abs(x); }

int main(int argc, char* argv[]) {
  if (argc!=3) {
    cout << "usage: " << argv[0] << " output_dir input.root\n";
    return 1;
  }

#define VAR_PREF "HGamEventInfoAuxDyn."
#define VAR_PREF_TRUTH "HGamTruthEventInfoAuxDyn."

  vardef<Char_t> v_isPassed(VAR_PREF "isPassed");
  vardef<Float_t, f_1e3> v_m_yy(VAR_PREF "m_yy");

  for (int argi=2; argi<argc; ++argi) {
    cout << argv[argi] << endl;
    TFile file(argv[argi]);
    if (file.IsZombie()) return 1;

    is_mc = false;
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
      is_mc = true;
      break;
    }
    if (!is_mc) {
      // get lumi
    }

    TTreeReader reader("CollectionTree",&file);
    reader_ptr = &reader;
    for (auto* var : vars) var->init();

    for (;;) {
      if (!*v_isPassed) continue;

      const auto myy = *v_m_yy;
      if (myy < 105 || 160 <= myy) continue;

      // if (is_mc) {
      //   weight = (**b_weight) * (**b_cs_br_fe) * mc_factor;
      // } else { // blind data in the signal region
      //   if (121 < myy && myy < 129) continue;
      // }
    }

    for (auto* var : vars) var->free();
  }
}
