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
//
// This file contains examples illustrating the use of the functions to compute
// the exact signs of sums of products.
//
//----------------------------------------------------------------------------------------

#include "gtest/gtest.h"
#include <cfenv>
#include <iostream>
#include "wm/wi/extra/float_sign.hpp"

//----------------------------------------------------------------------------------------
// In this sample we use "raw" functions, which assume that the rounding mode is upwards.
// They are more efficient than their "non raw" counterparts, but require attention to
// sertting/restoring the rounding mode, and lack attention to this detail may lead
// to bugs.

void raw_sample()
{
  int old = fegetround();
  std::fesetround(FE_UPWARD); // we need the rounding mode to be upwards

  float nu  = std::numeric_limits<float>::min();      // the smallest positive normal float number
  float eps = std::numeric_limits<float>::epsilon();  // the machine precission

  wm::quick_sign<float> quick;  // object for executing the quick algorithm
  wm::float_sign<float> robust; // object of  executing the robust algorithm

  float mat1[4 * 3] = {-1.0f, 4.0f * nu,  0.5f, 0.5f,
                        1.0f, 2.0f * nu, -1.0f, 0.5f * (1.0f + eps),
                        1.0f, 8.0f * nu,  0.25f,        1.0f + eps};

  std::pair<float const*,uint64_t> jagged1[3] =
  {
    std::pair<float const*,uint64_t>(mat1    , 4),
    std::pair<float const*,uint64_t>(mat1 + 4, 4),
    std::pair<float const*,uint64_t>(mat1 + 8, 4)
  };

  quick.clear();
  quick.raw_add<-1>(mat1[1], mat1[ 2], mat1[ 3]);  // -nu
  quick.raw_add< 1>(mat1[5], mat1[ 6], mat1[ 7]);  // -nu (1 + eps)
  quick.raw_add< 1>(mat1[9], mat1[10], mat1[11]);  //  2 * nu * (1 +eps)

  // computing the sign of  4 nu * 0.5 * 0.5 + 2 n *(-1) *0.5 * (1 + epsilon) +
  //                        8 * nu * 0.25 * (1 + epsilon)

  int sq = quick.raw_sign();
  if( sq != 0)       // 0 means that the quick algorithm failed.
    std::cout << "Bad quick sign!!" << std::endl;

  sq = wm::raw_quick_by_rows(mat1, 3, 4);
  if( sq != 0)      // 0 means that the quick algorithm failed.
    std::cout << "Bad quick sign by rows!!" << std::endl;

  sq = wm::raw_quick_jagged(jagged1, 3);
  if( sq != 0)      // 0 means that the quick algorithm failed.
    std::cout << "Bad quick sign by rows!!" << std::endl;

  //computing the same sign with the robust algorithm.

  robust.clear();
  robust.raw_add<-1>(mat1[1], mat1[ 2], mat1[ 3]);  // -nu
  robust.raw_add< 1>(mat1[5], mat1[ 6], mat1[ 7]);  // -nu (1 + eps)
  robust.raw_add< 1>(mat1[9], mat1[10], mat1[11]);  //  2 * nu * (1 +eps)

  int sr = robust.raw_sign();
  if( sr != 1 )    // the correct sign is +.
    std::cout << "Bad robust sign!!" << std::endl;

  sr = robust.raw_by_rows(mat1, 3, 4);
  if( sr != 1 )    // the correct sign is +.
    std::cout << "Bad robust by rows!!" << std::endl;

  sr = robust.raw_jagged(jagged1, 3);
  if( sr != 1 )    // the correct sign is +.
    std::cout << "Bad robust jagged!!" << std::endl;

  // one more time, changing signs and magnitudes

  float mat2[3* 4] = { 4.0f * nu,  0.5f,                0.5f,  1.0f,
                      -1.0f * nu, -1.0f,         (1.0f + eps), 1.0f,
                       8.0f * nu, -0.5f,  0.5f * (1.0f + eps), 1.0f};

  std::pair<float const*,uint64_t> jagged2[3] =
  {
    std::pair<float const*,uint64_t>(mat2    , 4),
    std::pair<float const*,uint64_t>(mat2 + 4, 4),
    std::pair<float const*,uint64_t>(mat2 + 8, 4)
  };

  quick.clear();
  quick.raw_add<1>(mat2[0], mat2[1], mat2[ 2]);  // -nu
  quick.raw_add<1>(mat2[4], mat2[5], mat2[ 6]);  // -nu (1 + eps)
  quick.raw_add<1>(mat2[8], mat2[9], mat2[10]);  //  2 * nu * (1 +eps)

  sq = quick.raw_sign();
  if( sq != 0)
     std::cout << "Bad quick sign!!" << std::endl;

  sq = wm::raw_quick_by_rows(mat2, 3, 4);
  if( sq != 0)
   std::cout << "Bad quick sign by rows!!" << std::endl;

  sq = wm::raw_quick_jagged(jagged2, 3);
  if( sq != 0)
   std::cout << "Bad quick jagged!!" << std::endl;

  robust.clear();
  robust.raw_add<1>(mat2[0], mat2[1], mat2[ 2]);  // -nu
  robust.raw_add<1>(mat2[4], mat2[5], mat2[ 6]);  // -nu (1 + eps)
  robust.raw_add<1>(mat2[8], mat2[9], mat2[10]);  //  2 * nu * (1 +eps)

  sr = robust.raw_sign();
  if( sr != -1 )
    std::cout << "Bad robust sign!!" << std::endl;

  sr = robust.raw_by_rows(mat2, 3, 4);
  if( sr != -1 )    // the correct sign is +.
    std::cout << "Bad robust by rows!!" << std::endl;

  sr = robust.raw_jagged(jagged2, 3);
  if( sr != -1 )    // the correct sign is +.
    std::cout << "Bad robust jagged!!" << std::endl;

  std::fesetround(old);
}

//----------------------------------------------------------------------------------------
// computing the sign of determinants. In this function, since we are not using
// raw functions, there is no need to set the rounding mode upwards

void det_sample()
{
  std::srand(10);
  constexpr uint n_max = 7;
  float mat[n_max * n_max];

  auto gen = []()
  {
    return (2.0f * std::rand())/RAND_MAX - 1.0;
  };

  wm::float_sign<float> robust;

  for(uint n = 1; n < n_max; ++n)
  {
    std::generate_n(mat, n * n, gen);
    int quick_det  = wm::quick_det_sign(mat, n);
    int robust_det = robust.det_sign(mat, n);
    if( quick_det && (robust_det != quick_det) )
      std::cout << "bad det sign" << std::endl;
  }
}

