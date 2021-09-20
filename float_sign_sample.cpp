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

#include <cfenv>
#include "wm/wi/extra/float_sign.hpp"

void sample()
{
  std::fesetround(FE_UPWARD); // we need the rounding mode to be upwards

  float nu  = std::numeric_limits<float>::min();
  float eps = std::numeric_limits<float>::epsilon();

  wm::quick_sign<float> quick;
  wm::float_sign<float> robust;

  quick.add<-1>(4.0f * nu,  0.5f,                0.5f);  // -nu
  quick.add< 1>(2.0f * nu, -1.0f, 0.5f * (1.0f + eps));  // -nu (1 + eps)
  quick.add< 1>(8.0f * nu,  0.25f,         1.0f + eps);  //  2 * nu * (1 +eps)

  int sq = quick.sign();
  assert( sq == 0);

  robust.add<-1>(4.0f * nu,  0.5f,                0.5f);  // -nu
  robust.add< 1>(2.0f * nu, -1.0f, 0.5f * (1.0f + eps));  // -nu (1 + eps)
  robust.add< 1>(8.0f * nu,  0.25f,         1.0f + eps);  //  2 * nu * (1 +eps)

  int sr = robust.sign();
  assert( sr == 1 );

  // one more time, changing signs and magnitudes

  quick.clear();
  robust.clear();

  quick.add< 1>( 4.0f * nu,  0.5f,                0.5f ); // nu
  quick.add< 1>(-1.0f * nu, -1.0f,         (1.0f + eps)); // nu (1 + eps)
  quick.add< 1>( 8.0f * nu, -0.5f,  0.5f * (1.0f + eps)); // -2 * nu * (1 +eps)

  sq = quick.sign();
  assert( sq == 0);

  robust.add< 1>( 4.0f * nu,  0.5f,               0.5f  ); // nu
  robust.add< 1>(-1.0f * nu, -1.0f,         (1.0f + eps)); // nu (1 + eps)
  robust.add< 1>(8.0f * nu,  -0.5f,  0.5f * (1.0f + eps)); // -2 * nu * (1 +eps)

  sr = robust.sign();
  assert( sr == -1 );

}
