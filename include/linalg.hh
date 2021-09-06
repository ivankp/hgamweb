#ifndef LIN_ALG_HH
#define LIN_ALG_HH

namespace ivanp {

template <typename... T>
[[ gnu::always_inline ]]
constexpr auto sq(const T&... x) noexcept { return (... + (x*x)); }

namespace linalg {

constexpr unsigned utn(unsigned n) noexcept { return n*(n+1) >> 1; }

void cholesky(double* A, unsigned N) noexcept;
void solve_triang(const double* L, double* v, unsigned n) noexcept;
void solve_triang_T(const double* L, double* v, unsigned n) noexcept;
void inv_triang(double* L, unsigned n) noexcept;
double dot(const double* a, const double* b, unsigned n) noexcept;

// void L_LT(const double* L, double* prod, unsigned n) noexcept;
void LT_L(double* L, unsigned n) noexcept;

void change_poly_coords(double* c, unsigned n, double a, double b) noexcept;

} }

#endif
