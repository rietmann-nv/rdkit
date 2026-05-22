//
//  Copyright 2001-2021 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

#include "types.h"

namespace RDKit {
namespace detail {
const PropToken computedPropNameToken = PropToken(computedPropName);
}

namespace common_properties {
const PropToken _hasMassQueryToken = PropToken(_hasMassQuery);

const PropToken _ChiralityPossibleToken = PropToken(_ChiralityPossible);
const PropToken _chiralPermutationToken = PropToken(_chiralPermutation);
const PropToken _CIPCodeToken = PropToken(_CIPCode);
const PropToken _CIPRankToken = PropToken(_CIPRank);
const PropToken _isotopicHsToken = PropToken(_isotopicHs);
const PropToken _MolFileBondEndPtsToken = PropToken(_MolFileBondEndPts);
const PropToken _MolFileBondAttachToken = PropToken(_MolFileBondAttach);
const PropToken _MolFileRLabelToken = PropToken(_MolFileRLabel);
const PropToken _ringStereoAtomsAllToken = PropToken("_ringStereoAtomsAll");
const PropToken _ringStereoAtomsBeginsToken = PropToken("_ringStereoAtomsBegins");
const PropToken _ringStereoGroupToken = PropToken("_ringStereoGroup");
const PropToken _supplementalSmilesLabelToken = PropToken(_supplementalSmilesLabel);
const PropToken _UnknownStereoToken = PropToken(_UnknownStereo);
const PropToken dummyLabelToken = PropToken(dummyLabel);
const PropToken isImplicitToken = PropToken(isImplicit);
const PropToken molAtomMapNumberToken = PropToken(molAtomMapNumber);
const PropToken molFileAliasToken = PropToken(molFileAlias);
const PropToken molFileValueToken = PropToken(molFileValue);

}  // namespace common_properties

//  template <typename T>
//  T larger_of(T arg1,T arg2) { return arg1>arg2 ? arg1 : arg2; };

void Union(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res) {
  res.clear();
  res = r1;
  for (auto ri : r2) {
    if (std::find(res.begin(), res.end(), ri) == res.end()) {
      res.push_back(ri);
    }
  }
}

void Intersect(const INT_VECT &r1, const INT_VECT &r2, INT_VECT &res) {
  res.clear();
  for (auto ri : r1) {
    if (std::find(r2.begin(), r2.end(), ri) != r2.end()) {
      res.push_back(ri);
    }
  }
}

void Union(const VECT_INT_VECT &rings, INT_VECT &res, const INT_VECT *exclude) {
  res.clear();
  auto nrings = static_cast<int>(rings.size());
  for (int id = 0; id < nrings; ++id) {
    if (exclude) {
      if (std::find(exclude->begin(), exclude->end(), id) != exclude->end()) {
        continue;
      }
    }
    for (const auto &ri : rings[id]) {
      if (std::find(res.begin(), res.end(), ri) == res.end()) {
        res.push_back(ri);
      }
    }
  }
}

int nextCombination(INT_VECT &comb, int tot) {
  int nelem = static_cast<int>(comb.size());
  int celem = nelem - 1;

  while (comb[celem] == (tot - nelem + celem)) {
    celem--;
    if (celem < 0) {
      return -1;
    }
  }

  unsigned int i;
  comb[celem] += 1;
  for (i = celem + 1; i < comb.size(); i++) {
    comb[i] = comb[i - 1] + 1;
  }
  return celem;
}
}  // namespace RDKit
