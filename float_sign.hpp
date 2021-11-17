#pragma once

//----------------------------------------------------------------------------------------
// Copyright (c) 2021 Walter Mascarenhas
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// The MPL has a "Disclaimer of Warranty" and a "Limitation of Liability",
// and we provide no guarantee regarding the correctness of this source code.
//
//----------------------------------------------------------------------------------------

#include <fenv.h>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <queue>
#include <quadmath.h>

//----------------------------------------------------------------------------------------
namespace wm {

//----------------------------------------------------------------------------------------
namespace float_sign_detail {

template <class I>
int next_permutation(I* permutation, uint64_t n, int sign)
{
  // Find longest non-increasing suffix
  uint64_t i = n - 1;
  while( (i > 0) && (permutation[i - 1] >= permutation[i]) )
    --i;

  // Now i is the head index of the suffix

  // Are we at the last permutation already?
  if (i <= 0)
    return -2;

  // Let permutation[i - 1] be the pivot
  // Find rightmost element greater than the pivot

  uint64_t j = n - 1;
  while (permutation[j] <= permutation[i - 1])
    j--;

  // Now the value permutation[j] will become the new pivot
  // Assertion: j >= i

  // Swap the pivot with j
  std::swap(permutation[i - 1], permutation[j]); // 1 change in sign

  // Reverse the suffix
  j = n - 1;

  if( ((j - i + 3)/2) & 0x1 )
    sign = -sign;

  while (i < j)                                 // (j - i + 1) / 2 changes in sign
  {
    std::swap(permutation[i], permutation[j]);
    i++;
    j--;
  }

  // Successfully computed the next permutation
  return sign;
}

template <class T>
struct traits;

template <class T>
requires (std::is_floating_point_v<T>)
struct traits<T>
{
  using value_type = T;

  static constexpr T zero()
  {
    return value_type(0);
  }

  static constexpr value_type tau()
  {
    return 2 * std::numeric_limits<value_type>::min() / std::numeric_limits<value_type>::epsilon();
  }

  static constexpr value_type sigma()
  {
    return 2 / std::numeric_limits<value_type>::min();
  }

  static constexpr value_type tau_sigma()
  {
    return tau() * sigma();
  }

  static constexpr value_type inv_sigma()
  {
    return std::numeric_limits<value_type>::min() * value_type{0.5};
  }

  static constexpr bool is_positive(value_type x)
  {
    return x > value_type{0.0};
  }

  static constexpr bool is_zero(value_type x)
  {
    return x == value_type{0.0};
  }

  static constexpr bool is_negative(value_type x)
  {
    return x < value_type{0.0};
  }

  static constexpr bool isfinite(value_type x)
  {
    return std::isfinite(x);
  }

  static constexpr bool is_greater_than_tau(value_type x)
  {
    return x > tau();
  }

  static constexpr bool is_greater_than_tau_sigma(value_type x)
  {
    return x > tau_sigma();
  }

  static constexpr value_type mul_by_sigma(value_type x)
  {
    return x * sigma();
  }

  static constexpr value_type div_by_sigma(value_type x)
  {
    return x * inv_sigma();
  }

  static value_type fma(value_type a, value_type x, value_type y)
  {
    return std::fma(a, x, y);
  }
};

template <>
struct traits<__float128>
{
  using value_type = __float128;

  static constexpr value_type zero()
  {
    return value_type(0);
  }

  static constexpr value_type tau()
  {
    return 2 * FLT128_MIN / FLT128_EPSILON;
  }

  static constexpr value_type sigma()
  {
    return 2 / FLT128_MIN;
  }

  static constexpr value_type tau_sigma()
  {
    return tau() * sigma();
  }

  static constexpr value_type inv_sigma()
  {
    return FLT128_MIN * 0.5;
  }

  static bool isfinite(value_type x)
  {
    return !isinfq(x);
  }

  static constexpr bool is_zero(value_type x)
  {
    return x == value_type{0.0};
  }

  static constexpr bool is_positive(value_type x)
  {
    return x > value_type{0.0};
  }

  static constexpr bool is_negative(value_type x)
  {
    return x < value_type{0.0};
  }

  static constexpr bool is_greater_than_tau(value_type x)
  {
    return x > tau();
  }

  static constexpr bool is_greater_than_tau_sigma(value_type x)
  {
    return x > tau_sigma();
  }

  static constexpr value_type mul_by_sigma(value_type x)
  {
    return x * sigma();
  }

  static constexpr value_type div_by_sigma(value_type x)
  {
    return inv_sigma() * x;
  }

  static value_type fma(value_type a, value_type x, value_type y)
  {
    return fmaq(a, x, y);
  }
};

} // wm::float_sign_detail
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// quick_sign tries to compute the sign of a sum of products sum prod aij quickly.
// It may fail, in which case it returns 0, or suceed, in which case it returns +1 or -1.
//
// The quick_sign<T> object qs is used as follows:
//
// 1) we clear qs.
// 2) For as sign in {-1,1}, we add each prod sign * ai1 * ai2 ... aik by calling
//    qs.add<sign>(ai1, ai2, .... aik )
// 3) we compute the sign of the sum of products by calling qs.sign()
//
// In the process above, we assume that the rounding mode is upwards

template <class T = double>
class quick_sign
{
  public:

    using traits = float_sign_detail::traits<T>;

    quick_sign()
    : u_(traits::zero()), md_(traits::zero()){}

    void clear()
    {
      u_ = 0;
      md_ = 0;
    }

    int raw_sign()
    {
      return (traits::is_negative(md_) ? 1 : (traits::is_negative(u_) ? -1: 0));
    }

    template <int Sign = 1, class... X>
    void raw_add(T const& a, X const&... x) // adding the product a * x... to *this
    {
      constexpr int sign = Sign;
      if constexpr (sizeof...(X) == 0)
      {
        if constexpr (sign > 0)
        {
          u_  +=  a;
          md_ += -a;
        }
        else
        {
          u_  += -a;
          md_ +=  a;
        }
      }
      else
      {
        if( traits::is_positive(a) )
          add_impl<sign>(-a, a, x...);
        else
        {
          if( traits::is_negative(a) )
            add_impl<-sign>(a, -a, x...);
        }
      }
    }

private:

    template <int Sign, class... X>                     // prod = sign *|a0 ... ak| * ak+1
    void add_impl(T md, T u, T const& a, X const&... x) //  u = |a0 ... ak| + ek
    {                                                   // md = (ek - |a0 ....ak|) => ek * a - |a0... ak+1| + ek
      constexpr int sign = Sign;

      if( traits::is_positive(a) )
      {
        if constexpr (sizeof...(x) > 0)
          add_impl<sign>(md * a, u * a, x...);
        else
        {
          if constexpr (sign > 0)
          {
            u_  = traits::fma(u,  a,  u_);
            md_ = traits::fma(md, a, md_);
          }
          else
          {
            u_  = traits::fma(md, a,  u_);
            md_ = traits::fma( u, a, md_);
          }
        }
      }
      else
      {
        if( traits::is_negative(a) )
        {
          T ma = -a;
          if constexpr (sizeof...(x) > 0)
            add_impl<-sign>(md * ma, u * ma, x...);
          else
          {
            if constexpr (sign > 0)
            {
              md_ = traits::fma(u,  ma, md_);
              u_  = traits::fma(md, ma,  u_);
            }
            else
            {
              u_  = fma(u,  ma,   u_);
              md_ = fma(md, ma,  md_);
            }
          }
        }
      }
    }

    T u_;
    T md_;
};

//----------------------------------------------------------------------------------------
// quick_sign_for_rows computes the sign of
// S = sign * sum_{i=0}}^{rows - 1} prod_{j = 0}^{n-1} x[j + i * cols]
// As in the other functions in this file, it assumes that the rounding mode is upwards,
// It returns:
//
// -1 if S is certainly negative
//  0 if is cannot decide which sign S has.
//  1 if S is certainly positive

template <class T = double>
void raw_add_quick_product(T& u, T& md, T const* x, uint64_t n, int sign)
{
  assert( n >= 2 );

  T abs_p;
  T minus_abs_p;
  if( float_sign_detail::traits<T>::is_negative(*x) )
  {
    abs_p = -*x;
    sign = -sign;
    minus_abs_p = *x++;
  }
  else
  {
    abs_p = *x;
    minus_abs_p = -*x++;
  }

  for(uint64_t j = 2; j < n; ++j)
  {
    if( float_sign_detail::traits<T>::is_negative(*x) )
    {
      T mx = -*x++;
      sign = -sign;
      abs_p *= mx;
      minus_abs_p *= mx;
    }
    else
    {
      abs_p *= *x;
      minus_abs_p *= *x++;
    }
  }

  if( float_sign_detail::traits<T>::is_negative(*x) )
  {
    T mx = -*x++;
    //sign = -sign;
    if( sign < 0 )
    {
      u  = fma(abs_p,       mx, u);
      md = fma(minus_abs_p, mx, md);
//      return std::array<T,2>{abs_p, minus_abs_p};
    }
    else
    {
      u  = fma(minus_abs_p, mx, u);
      md = fma(abs_p,       mx, md);
//   return std::array<T,2>{minus_abs_p, abs_p};
    }
//    abs_p *= *x;
//    minus_abs_p *= *x++;
//    abs_p *= mx;
//    minus_abs_p *= mx;
  }
  else
  {
    if( sign > 0 )
    {
      u  = fma(abs_p,       *x, u);
      md = fma(minus_abs_p, *x, md);
//      return std::array<T,2>{abs_p, minus_abs_p};
    }
    else
    {
      u  = fma(minus_abs_p, *x, u);
      md = fma(abs_p,       *x, md);
//   return std::array<T,2>{minus_abs_p, abs_p};
    }
//    abs_p *= *x;
//    minus_abs_p *= *x++;
  }
}

template <class T = double>
int raw_quick_by_rows(T const* x, uint64_t rows, uint64_t cols, int sign = 1)
{
  T u( float_sign_detail::traits<T>::zero());             // sum of products,   rounded up
  T md(u);           // (-sum_of_products) rounded up = - (sum_of_products rounded down)

  if( cols < 2 )
  {
    if( cols == 0 )
      return 0;

    for(uint64_t i = 0; i < rows; ++i)
    {
      u  +=  x[i];
      md += -x[i];
    }
  }
  else
  {
    for(uint64_t i = 0; i < rows; ++i, x += cols)
      raw_add_quick_product(u, md, x, cols, sign);
  }

  if( float_sign_detail::traits<T>::is_negative(u) )
    return -1;

  if( float_sign_detail::traits<T>::is_negative(md) )
    return 1;

  return 0;
}

template <class T = double>
int quick_by_rows(T const* x, uint64_t rows, uint64_t cols, int sign = 1)
{
  int old = std::fegetround();
  std::fesetround(FE_UPWARD);
  int ret_sign = raw_quick_by_rows(x, rows, cols, sign);
  std::fesetround(old);
  return ret_sign;
}

template <class T = double>
int raw_quick_jagged(std::pair<T const*,uint64_t> const* x, uint64_t rows, int sign = 1)
{
  T u(float_sign_detail::traits<T>::zero());            // sum of products,   rounded up
  T md(u);           // (-sum_of_products) rounded up = - (sum_of_products rounded down)

  for(uint64_t i = 0; i < rows; ++i, ++x)
  {
    if( x->second < 2)
    {
      if( x->second == 1 )
      {
        u  +=  x->first[0];
        md += -x->first[0];
      }
    }
    else
    {
      raw_add_quick_product(u, md, x->first, x->second, sign);
    }
  }

  if( float_sign_detail::traits<T>::is_negative(u) )
    return -1;

  if( float_sign_detail::traits<T>::is_negative(md) )
    return 1;

  return 0;
}

template <class T = double>
int quick_jagged(std::pair<T const*,uint64_t> const* x, uint64_t rows, int sign = 1)
{
  int old = std::fegetround();
  std::fesetround(FE_UPWARD);
  int ret_sign = raw_quick_jagged(x, rows, sign);
  std::fesetround(old);
  return ret_sign;
}

template <class T>
int raw_quick_det_sign(T const* x, uint64_t n)
{
  if( n < 2)
  {
    if( n == 0)
      return 0;

    if( float_sign_detail::traits<T>::is_negative(x[0]) )
      return -1;

    if( float_sign_detail::traits<T>::is_positive(x[0]) )
      return 1;

    return 0;
  }

  T parcel[n];
  int permut_sign = 1;
  uint permutation[n];
  std::iota(permutation, permutation + n, 0);

  T u(float_sign_detail::traits<T>::zero());
  T md(u);

  do
  {
    for(uint k = 0; k < n; ++k)
      parcel[k] = x[ k + permutation[k] * n];

    raw_add_quick_product(u, md, parcel, n, permut_sign);
    permut_sign = float_sign_detail::next_permutation(permutation, n, permut_sign);
  }
  while( permut_sign > -2 );

  if( float_sign_detail::traits<T>::is_negative(u) )
    return -1;

  if( float_sign_detail::traits<T>::is_negative(md) )
    return 1;

  return 0;
}

template <class T>
int quick_det_sign(T const* x, uint64_t n)
{
  int old = std::fegetround();
  std::fesetround(FE_UPWARD);
  int sign = raw_quick_det_sign(x, n);
  std::fesetround(old);
  return sign;
}

//----------------------------------------------------------------------------------------
// float_sign computes the sign of a sum of products sum prod aij quickly.
// It returns -1, 0 or 1
//
// The float_sign<T> object fs is used as follows:
//
// 1) we clear fs.
// 2) For as sign in {-1,1},  we add each prod sign * ai1 * ai2 ... aik by calling
//    fs.add<sign>(ai1, ai2, .... aik )
// 3) we compute the sign of the sum of products by calling fs.sign()
//
// It can also be used to compute the sign of products in which the factors are rows
// of a matrix, that is, to compute sign of
//
// S =  sum_{i = 0}^(rows - 1) prod_{j=0}^{cols - 1} x[j + i * cols]
//
// This computation is performed by invoking:
//
// float_sign.by_rows(x, rows,cols)
//
// In the processes above, we assume that the rounding mode is upwards
//

template <class T = double>
class float_sign
{
  public:

    using traits = float_sign_detail::traits<T>;

    int raw_jagged(std::pair<T const*, uint64_t>* parcel, uint64_t n_parcels, int sign = 1)
    {
      clear();
      for(uint64_t i = 0; i < n_parcels; ++i, ++parcel)
      {
        if( parcel->second == 0 )
          continue;

        split_product(parcel->first, parcel->second, sign);
        for(auto& s : plus)
        {
          s.normalize();
          queue_[0].emplace(s);
        }

        for(auto& s : minus)
        {
          s.normalize();
          queue_[1].emplace(s);
        }
      }

      return this->raw_sign();
    }

    int jagged(std::pair<T const*, uint64_t>* parcel, uint64_t n_parcels, int sign = 1)
    {
      int old = std::fegetround();
      std::fesetround(FE_UPWARD);
      int ret_sign = raw_jagged(parcel, n_parcels, sign);
      std::fesetround(old);
      return ret_sign;
    }

    int raw_by_rows(T const* x, uint64_t rows, uint64_t cols, int sign = 1)
    {
     if( cols == 0 )
        return 0;

      clear();
      for(uint64_t i = 0; i < rows; ++i, x += cols)
      {
        split_product(x, cols, sign);
        for(auto& s : plus)
        {
          s.normalize();
          queue_[0].push(s);
        }

        for(auto& s : minus)
        {
          s.normalize();
          queue_[1].push(s);
        }
      }
      return this->raw_sign();
    }

    int by_rows(T const* x, uint64_t rows, uint64_t cols, int sign = 1)
    {
      int old = std::fegetround();
      std::fesetround(FE_UPWARD);
      int ret_sign = raw_by_rows(x, rows, cols, sign);
      std::fesetround(old);
      return ret_sign;
    }

    int raw_det_sign(T const* x, uint64_t n)
    {
      if( n < 2)
      {
        if( n == 0)
          return 0;

        if( traits::is_negative(x[0]) )
          return -1;

        if( traits::is_positive(x[0]) )
          return 1;

        return 0;
      }

      clear();
      T parcel[n];
      int permut_sign = 1;
      uint permutation[n];
      std::iota(permutation, permutation + n, 0);

      do
      {
        for(uint k = 0; k < n; ++k)
          parcel[k] = x[ k + permutation[k] * n];

        split_product(parcel, n, permut_sign);
        for(auto& s : plus)
        {
          s.normalize();
          queue_[0].push(s);
        }

        for(auto& s : minus)
        {
          s.normalize();
          queue_[1].push(s);
        }

        permut_sign = float_sign_detail::next_permutation(permutation, n, permut_sign);
      }
      while( permut_sign > -2 );
      return this->raw_sign();
    }

    int det_sign(T const* x, uint64_t n)
    {
      int old = std::fegetround();
      std::fesetround(FE_UPWARD);
      int ret_sign = raw_det_sign(x, n);
      std::fesetround(old);
      return ret_sign;
    }

    int raw_sign()
    {
      while( true )
      {
        if( queue_[0].empty() )
          return (queue_[1].empty() ? 0 : -1);

        if( queue_[1].empty() )
          return 1;

        scaled p = queue_[0].top();
        scaled m = queue_[1].top();

        if( p.exp >= m.exp )
        {
          if( p.exp == m.exp )
          {
            if( p.t > m.t )
            {
              if( (p.t > queue_[1].size() * m.t) || diff<0>(p.t, m.t, p.exp) )
                return 1;
            }
            else
            {
              if( (m.t > queue_[0].size() * p.t) || diff<1>(m.t, p.t, p.exp)  )
                return -1;
            }
          }
          else
          {
            if( (p.exp > m.exp + 1) || diff<0>(p.t, traits::div_by_sigma(m.t), p.exp) )
              return 1;
          }
        }
        else
        {
          if( (m.exp > p.exp + 1) || diff<1>(m.t, traits::div_by_sigma(p.t), m.exp) )
            return -1;
        }
      }
    }


    template <int Sign, class X, class... Y>
    requires ( ((Sign == -1) || (Sign ==1)) && requires(X x, Y... y)
    {
      {T(x)};
      {(T(y),...)};
    })
    void raw_add(X const& x, Y&... y)
    {
      constexpr int sign_param = Sign;
      T tx(x);
      int sign;
      plus.clear();
      minus.clear();

      if( traits::is_positive(tx) )
      {
        sign = sign_param;
        plus.emplace_back(std::move(tx), 0);
      }
      else
      {
        if( traits::is_zero(tx) )
          return;

        sign = -sign_param;
        plus.emplace_back(-tx, 0);
      }

      bool is_non_zero = (push_y(sign, y) && ...);
      if( !is_non_zero )
        return;

      if( sign < 0 )
        swap(plus, minus);

      for(auto& s : plus)
      {
        s.normalize();
        queue_[0].emplace( std::move(s) );
      }

      for(auto& s : minus)
      {
        s.normalize();
        queue_[1].emplace( std::move(s) );
      }
    }

    void clear()
    {
      queue_[0].clear();
      queue_[1].clear();
    }

  private:

    bool push_y(int& sign, T y)
    {
      if( traits::is_negative(y) )
      {
        y = -y;
        sign = -sign;
      }
      else
      {
        if( traits::is_zero(y) )
          return false;
      }

      uint64_t n_minus = minus.size();

      for(auto& s : plus)
        update(s, y, minus);

      for(uint64_t k = 0; k < n_minus; ++k)
         update(minus[k], y, plus);

      return true;
    }

    struct scaled // value = sigma^exp * t, with  tau < t <= sigma * tau
    {
      T t;
      int exp;

      bool operator<(scaled const& x) const
      {
        if( exp < x.exp )
          return true;

        if( exp > x.exp )
          return false;

        return t < x.t;
      }

      void normalize()
      {
        if( traits::is_greater_than_tau(t) )
        {
          if( traits::is_greater_than_tau_sigma(t) )
          {
            ++exp;
            t = traits::div_by_sigma(t);
          }
        }
        else  // x <= tau
        {
          --exp;
          t = traits::mul_by_sigma(t);
        }
      }
    };

    void update(scaled& s, T xk, std::vector<scaled>& opposite)
    {
      while( true )
      {
        T p = s.t * xk;
        if( traits::is_greater_than_tau(p) )
        {
          if( traits::isfinite(p) )
          {
            T d = traits::fma(-xk, s.t, p);
            if( !traits::is_zero(d) )
            {
              assert( d > 0 );
              opposite.emplace_back(std::move(d), s.exp);
            }
            s.t = p;
            return;
          }

          s.exp += 1;
          if( s.t >= xk )
            s.t = traits::div_by_sigma(s.t);
          else
            xk = traits::div_by_sigma(xk);
        }
        else
        {
          s.exp -= 1;
          if( s.t <= xk )
            s.t = traits::mul_by_sigma(s.t);
          else
            xk = traits::mul_by_sigma(xk);
        }
      }
    }

    void split_product(T const* x, uint64_t n, int sign)
    {
      T xk;
      plus.clear();
      minus.clear();

      if( traits::is_positive( x[0] ) )
        plus.emplace_back(x[0], 0);
      else
      {
        if( traits::is_zero( x[0] ) )
          return;

        sign = -sign;
        plus.emplace_back(-x[0], 0);
      }

      for(uint64_t j = 1; j < n; ++j)
      {
        if( traits::is_positive(x[j]) )
          xk =  x[j];
        else
        {
          if( traits::is_zero(x[j]) )
          {
            plus.clear();
            minus.clear();
            return;
          }

          xk = -x[j];
          sign = -sign;
        }

        uint64_t n_minus = minus.size();
        for(auto& s : plus )
          update(s, xk, minus);

        for(uint64_t k = 0; k < n_minus; ++k)
          update(minus[k], xk, plus);
      }

      if( sign < 0 )
        swap(plus, minus);
    }

    template <int Lead>
    bool diff(T const& b, T const& a, int exp)
    {
       constexpr int lead = Lead;

      T c = b - a; // c = b - a + eps
      T d = b - c; // d = a - eps

      if( traits::is_zero(d) )
        return true;

      T e = a - d; // e = eps => b - a = c - e
      queue_[0].pop();
      queue_[1].pop();

      if( !traits::is_zero(c) )
      {
        scaled sc{c, exp};
        sc.normalize();
        queue_[lead].emplace(std::move(sc));

        if( !traits::is_zero(e) )
        {
          scaled se{e, exp};
          se.normalize();
          queue_[(1 + lead) & 0x1].emplace(std::move(se));
        }
      }
      return false;
    }

    struct queue : public std::priority_queue<scaled>
    {
      void clear()
      {
        this->c.clear();
      }
    };

    queue queue_[2];
    std::vector<scaled> plus;
    std::vector<scaled> minus;
};

} // namespace wm
//----------------------------------------------------------------------------------------

