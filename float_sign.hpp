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

#include <cassert>
#include <cmath>
#include <cstdint>
#include <queue>
#include <quadmath.h>

//----------------------------------------------------------------------------------------
namespace wm {

//----------------------------------------------------------------------------------------
namespace float_sign_detail {

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

  static constexpr bool is_greater_than_tau(value_type x)
  {
    return x > tau();
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

  static constexpr value_type inv_sigma()
  {
    return FLT128_MIN * 0.5;
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

    int sign()
    {
      return (traits::is_negative(md_) ? 1 : (traits::is_negative(u_) ? -1: 0));
    }

    template <int Sign = 1, class... X>
    void add(T const& a, X const&... x) // adding the product a * x... to *this
    {
      constexpr int sign = Sign;
      if( traits::is_positive(a) )
        add_impl<sign>(-a, a, x...);
      else
      {
        if( traits::is_negative(a) )
          add_impl<-sign>(a, -a, x...);
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
            u_  += u  * a;
            md_ += md * a;
          }
          else
          {
            md_ += u  * a;
            u_  += md * a;
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
              md_ += u  * ma;
              u_  += md * ma;
            }
            else
            {
              u_  += u  * ma;
              md_ += md * ma;
            }
          }
        }
      }
    }

    T u_;
    T md_;
};

template <class T = double>
class float_sign
{
  public:

    using traits = float_sign_detail::traits<T>;

    template <int Sign>
    requires  ((Sign == 1) || (Sign == -1))
    void add(T const& a) // adding Sign * a to *this
    {
      if( traits::is_greater_than_tau( a ) )
        push< Sign, 0>( T{a} );
      else
      {
        if( traits::is_negative(a) )
        {
          T ma = -a;
          if( traits::is_greater_than_tau( ma ) )
            push< -Sign, 0>( std::move(ma) );
          else
            push< -Sign, 1>( traits::mul_by_sigma(ma) );
        }
        else
        {
          if( traits::is_positive( a ) )
            push< Sign, 1>( traits::mul_by_sigma(a) );
        }
      }
    }

    template <int Sign, class... U>
    requires ( ((Sign == 1) || (Sign == -1)) && (std::is_constructible_v<T,U> &&...) )
    void add(T const& a, T const& b, U const&... u) // adding the product a  * u... to *this
    {
      if( traits::is_positive(a) )
        add_impl< Sign, 0, positive, unknown>(a, b, u...);
      else
      {
        if( traits::is_negative(a) )
          add_impl<-Sign,0, positive, unknown>(-a, b, u...);
      }
    }

    int sign()
    {
      while( true )
      {
        if( q< -1 >().empty() )
        {
          if( q< 1 >().empty() )
            return check_qs();
          else
            return check_q<1>();
        }

        if( q< 1 >().empty() )
          return check_q<-1>();

        T const& p = q<  1 >().top();
        T const& m = q< -1 >().top();

        if( p > m )
        {
          if( diff< 1 >(p, m) )
            return 1;
        }
        else
        {
          if( diff< -1 >(m, p) )
            return -1;
        }
      }
    }

    void clear()
    {
      q<  -1 >().clear();
      q<   1 >().clear();
      qs< -1 >().clear();
      qs<  1 >().clear();
    }

  private:

    static constexpr int large = 2;
    static constexpr int unknown = 0;
    static constexpr int positive = 1;

    template <int Sign, int Exp>
    void push(T&& a)
    {
      if constexpr (Exp == 0)
        q< Sign >().emplace( std::move(a) );
      else
        qs< Sign >().emplace( std::move(a), Exp);
    }

    template <int Sign, int Exp, int A, int B = unknown>
    requires ((A == large) || (A == positive))
    void add_impl(T const& a)
    {
      constexpr int size_a = A;
      if( (size_a == large) || traits::is_greater_than_tau( a ) )
        push<Sign, Exp    >( T{a} );
      else
        push<Sign, Exp + 1>( traits::mul_by_sigma(a) );
    }

    template <int Sign, int Exp, int A, int B, class... U>
    requires ( ((A == large)   || ((A == positive) && (B != large))) &&
               ((B == unknown) || (B == positive) || (B == large) ) )
    void add_impl(T const& a, T const& b, U const&... u)
    {
      constexpr int size_a = A;
      constexpr int size_b = B;
      if constexpr (size_b == unknown)
      {
        if( traits::is_positive(b) )
          add_impl< Sign, Exp, A, positive>(a, b, u...);
        else
        {
          if( traits::is_negative(b) )
            add_impl<-Sign, Exp, A, positive>(a, -b, u...);
        }
      }
      else
      {
        T p = a * b;
        if( ((size_a == large) && (size_b == large)) || traits::is_greater_than_tau(p) )
        {
          T f = fma(-a, b, p); // -a * b + p

          if( traits::is_positive(f) )
            add_impl<-Sign, Exp, positive, unknown>(f, u...);

          add_impl< Sign, Exp, positive, unknown>(p, u...);
        }
        else
        {
          if constexpr (size_a == large)
            add_impl< Sign, Exp, large, large>(a, traits::mul_by_sigma(b), u...);
          else
          {
            if( a <= b )
              add_impl<Sign, Exp + 1, large, positive>(traits::mul_by_sigma(a), b, u...);
            else
              add_impl<Sign, Exp + 1, large, positive>(traits::mul_by_sigma(b), a, u...);
          }
        }
      }
    }

    template <int Sign>
    int diff(T const& b, T const& a)
    {
      T c = b - a; // c = b - a + eps
      T d = b - c; // d = a - eps

      if( traits::is_zero(d) )
        return true;

      T e = a - d; // eps
      q<  Sign >().pop();
      q< -Sign >().pop();

      if( traits::is_positive(c) )
        add_impl<Sign, 0, positive>(c);

      if( traits::is_positive(e) )
          add_impl< -Sign, 0, positive>(e);

      return false;
    }

    template <int Sign>
    int check_q()
    {
      while( !q< Sign >().empty() )
      {
        T const& p = q<Sign>().top();
        if( p > traits::tau() * q< -Sign >().size() )
          return Sign;

        scaled s = qs<-Sign>().top();
        if( s.exp > 1 )
          return Sign;

        q<  Sign>().pop();
        qs<-Sign>().pop();
        if( diff<Sign>(p, traits::div_by_sigma(s.t)) )
          return Sign;
      }

      return check_qs();
    }

    int check_qs()
    {
      while( true )
      {
        if( qs< -1 >().empty() )
          return (qs< 1 >().empty() ? 0 : 1);

        if( qs< 1 >().empty() )
          return -1;

        scaled p = qs<  1 >().top();
        scaled m = qs< -1 >().top();
        if( p.exp <= m.exp )
        {
          if( p.exp == m.exp )
          {
            if( p.t > m.t )
            {
              if( diff_s<1>(p.t, m.t, p.exp) )
                return 1;
            }
            else
            {
              if( diff_s<-1>(m.t, p.t, p.exp) )
                return 1;
            }
          }
          else
          {
            if( (p.exp > m.exp + 1) || diff_s< 1>(p.t, traits::div_by_sigma(m.t), p.exp) )
              return 1;
          }
        }
        else
        {
          if( (m.exp < p.exp + 1) || diff_s<-1>(m.t, traits::div_by_sigma(p.t), m.exp) )
            return -1;
        }
      }
    }

    template <int Sign>
    bool diff_s(T const& b, T const& a, int exp)
    {
      constexpr int sign = Sign;

      T c = b - a; // c = b - a + eps
      T d = b - c; // d = a - eps

      if( traits::is_zero(d) )
        return true;

      T e = a - d; // e = eps => b - a = c - e
      qs<  1>().pop();
      qs< -1>().pop();

      if( traits::is_greater_than_tau(c) )
        qs< sign >().emplace( std::move(c), exp);
      else
      {
        if( traits::is_positive(c) )
          qs< sign >().emplace( traits::mul_by_sigma(c), exp + 1);
      }

      if( traits::is_greater_than_tau(e) )
        qs< -sign >().emplace( std::move(e), exp);
      else
      {
        if( traits::is_positive(e) )
          qs< -sign >().emplace( traits::mul_by_sigma(e), exp + 1);
      }

      return false;
    }

    struct scaled // value = nu^exp * t,  with nu <= t <= 1
    {
      T t;
      int exp;

      bool operator<(scaled const& x) const
      {
        if( exp > x.exp )
          return true;

        if( exp < x.exp )
          return false;

        return t < x.t;
      }
    };

    struct queue : public std::priority_queue<T>
    {
      void clear()
      {
        this->c.clear();
      }
    };

    template <int Sign>
    queue& q()
    {
      return queue_[ (Sign + 1)/2 ];
    }

    struct queue_s : public std::priority_queue<scaled>
    {
      void clear()
      {
        this->c.clear();
      }
    };

    template <int Sign>
    queue_s& qs()
    {
      return queue_s_[ (Sign + 1)/2 ];
    }

    queue queue_[2];
    queue_s queue_s_[2];
};


} // namespace wm
//----------------------------------------------------------------------------------------

