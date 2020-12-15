#ifndef VECT_HELPERS_HH
#define VECT_HELPERS_HH

// Declare a member operator which takes 0 arguments

#define DECLARE_MEMBER_OPERATOR_0(fn, op)                                      \
                                                                               \
  vect fn() const {                                                            \
    vect r;                                                                    \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = op elt[d];                                                        \
    return r;                                                                  \
  }

// Declare a member operator which takes 0 arguments and returns type
// R

#define DECLARE_MEMBER_OPERATOR_0_RET(fn, op, R)                               \
                                                                               \
  vect<R, D> fn() const {                                                      \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = op elt[d];                                                        \
    return r;                                                                  \
  }

// Declare a member operator which takes 1 argument and returns a
// reference

#define DECLARE_MEMBER_OPERATOR_1_REF(fn, op)                                  \
                                                                               \
  vect &fn(const T &x) {                                                       \
    for (int d = 0; d < D; ++d)                                                \
      elt[d] op x;                                                             \
    return *this;                                                              \
  }                                                                            \
                                                                               \
  vect &fn(const vect &a) {                                                    \
    for (int d = 0; d < D; ++d)                                                \
      elt[d] op a[d];                                                          \
    return *this;                                                              \
  }

// Declare a function which takes 1 argument and returns type R

#define DECLARE_FUNCTION_1_RET(fn, R)                                          \
                                                                               \
  template <typename T, int D> inline vect<R, D> fn(const vect<T, D> &a) {     \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = fn(a[d]);                                                         \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D, int E>                                          \
  inline vect<R, D> fn(const vect<vect<T, D>, E> &a) {                         \
    vect<R, D> r;                                                              \
    for (int e = 0; e < E; ++e)                                                \
      r[e] = fn(a[e]);                                                         \
    return r;                                                                  \
  }

// Declare a function which takes 1 argument

#define DECLARE_FUNCTION_1(fn) DECLARE_FUNCTION_1_RET(fn, T)

// Declare a function which takes 2 arguments and returns type R

#define DECLARE_FUNCTION_2_RET(fn, R)                                          \
                                                                               \
  template <typename T, int D>                                                 \
  inline vect<R, D> fn(const vect<T, D> &a, const vect<T, D> &b) {             \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = fn(a[d], b[d]);                                                   \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D>                                                 \
  inline vect<R, D> fn(const T &a, const vect<T, D> &b) {                      \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = fn(a, b[d]);                                                      \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D>                                                 \
  inline vect<R, D> fn(const vect<T, D> &a, const T &b) {                      \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = fn(a[d], b);                                                      \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D, int E>                                          \
  inline vect<vect<R, D>, E> fn(const vect<vect<T, D>, E> &a,                  \
                                const vect<vect<T, D>, E> &b) {                \
    vect<vect<R, D>, E> r;                                                     \
    for (int e = 0; e < E; ++e)                                                \
      r[e] = fn(a[e], b[e]);                                                   \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D, int E>                                          \
  inline vect<vect<R, D>, E> fn(const T &a, const vect<vect<T, D>, E> &b) {    \
    vect<vect<R, D>, E> r;                                                     \
    for (int e = 0; e < E; ++e)                                                \
      r[e] = fn(a, b[e]);                                                      \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D, int E>                                          \
  inline vect<vect<R, D>, E> fn(const vect<vect<T, D>, E> &a, const T &b) {    \
    vect<vect<R, D>, E> r;                                                     \
    for (int e = 0; e < E; ++e)                                                \
      r[e] = fn(a[e], b);                                                      \
    return r;                                                                  \
  }

// Declare a function which takes 2 arguments

#define DECLARE_FUNCTION_2(fn) DECLARE_FUNCTION_2_RET(fn, T)

// Declare an operator which takes 1 argument and returns type R

#define DECLARE_OPERATOR_1_RET(fn, op, R)                                      \
                                                                               \
  template <typename T, int D> inline vect<R, D> fn(const vect<T, D> &a) {     \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = op a[d];                                                          \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D> inline vect<R, D> fn(const T &a) {              \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = op a;                                                             \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D, int E>                                          \
  inline vect<vect<R, D>, E> fn(const vect<vect<T, D>, E> &a) {                \
    vect<vect<R, D>, E> r;                                                     \
    for (int e = 0; e < E; ++e)                                                \
      r[e] = op a[e];                                                          \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D, int E>                                          \
  inline vect<vect<R, D>, E> fn(const T &a) {                                  \
    vect<vect<R, D>, E> r;                                                     \
    for (int e = 0; e < E; ++e)                                                \
      r[e] = op a;                                                             \
    return r;                                                                  \
  }

// Declare an operator which takes 2 arguments and returns type R

#define DECLARE_OPERATOR_2_RET(fn, op, R)                                      \
                                                                               \
  template <typename T, int D>                                                 \
  inline vect<R, D> fn(const vect<T, D> &a, const vect<T, D> &b) {             \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = a[d] op b[d];                                                     \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D>                                                 \
  inline vect<R, D> fn(const T &a, const vect<T, D> &b) {                      \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = a op b[d];                                                        \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D>                                                 \
  inline vect<R, D> fn(const vect<T, D> &a, const T &b) {                      \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = a[d] op b;                                                        \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D, int E>                                          \
  inline vect<vect<R, D>, E> fn(const vect<vect<T, D>, E> &a,                  \
                                const vect<vect<T, D>, E> &b) {                \
    vect<vect<R, D>, E> r;                                                     \
    for (int e = 0; e < E; ++e)                                                \
      r[e] = a[e] op b[e];                                                     \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D, int E>                                          \
  inline vect<vect<R, D>, E> fn(const T &a, const vect<vect<T, D>, E> &b) {    \
    vect<vect<R, D>, E> r;                                                     \
    for (int e = 0; e < E; ++e)                                                \
      r[e] = a op b[e];                                                        \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <typename T, int D, int E>                                          \
  inline vect<vect<R, D>, E> fn(const vect<vect<T, D>, E> &a, const T &b) {    \
    vect<vect<R, D>, E> r;                                                     \
    for (int e = 0; e < E; ++e)                                                \
      r[e] = a[e] op b;                                                        \
    return r;                                                                  \
  }

// Declare an operator which takes 2 arguments

#define DECLARE_OPERATOR_2(fn, op) DECLARE_OPERATOR_2_RET(fn, op, T)

// Declare a reduction function which takes 1 argument of type T and
// returns type R

#define DECLARE_REDUCTION_OPERATOR_1_T_RET(fn, init, op, final, T, R)          \
                                                                               \
  template <typename U, int D> inline vect<R, D> fn(const vect<U, D> &a) {     \
    vect<R, D> r;                                                              \
    for (int d = 0; d < D; ++d)                                                \
      r[d] = fn(a[d]);                                                         \
    return r;                                                                  \
  }                                                                            \
                                                                               \
  template <int D> inline R fn(const vect<T, D> &a) {                          \
    R r(init);                                                                 \
    for (int d = 0; d < D; ++d)                                                \
      r op a[d];                                                               \
    return final(r);                                                           \
  }

// Declare a reduction function which takes 1 argument

#define DECLARE_REDUCTION_OPERATOR_1(fn, init, op, final)                      \
                                                                               \
  template <typename T, int D> inline T fn(const vect<T, D> &a) {              \
    T r(init);                                                                 \
    for (int d = 0; d < D; ++d)                                                \
      r op a[d];                                                               \
    return final(r);                                                           \
  }

// Declare a reduction function which takes 1 argument

#define DECLARE_REDUCTION_FUNCTION_1(fn, init, op, final)                      \
                                                                               \
  template <typename T, int D> inline T fn(const vect<T, D> &a) {              \
    T r(init);                                                                 \
    for (int d = 0; d < D; ++d)                                                \
      op(r, a[d]);                                                             \
    return final(r);                                                           \
  }

// Declare a reduction function which takes 2 arguments

#define DECLARE_REDUCTION_OPERATOR_2(fn, init, op, op2, final)                 \
                                                                               \
  template <typename T, int D>                                                 \
  inline T fn(const vect<T, D> &a, const vect<T, D> &b) {                      \
    T r(init);                                                                 \
    for (int d = 0; d < D; ++d)                                                \
      r op(a[d] op2 b[d]);                                                     \
    return final(r);                                                           \
  }

#endif // #ifndef VECT_HELPERS_HH
