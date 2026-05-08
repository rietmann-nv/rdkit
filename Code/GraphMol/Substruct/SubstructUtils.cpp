//
//  Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SubstructUtils.h"
#include <set>
#include <RDGeneral/utils.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/RDMol.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

namespace {

//! Native equivalent of Atom::Match(const Atom*) operating directly on
//! AtomData, avoiding the round-trip through the compatibility wrappers.
bool atomDataMatch(const AtomData &qAtom, const AtomData &mAtom) {
  if (qAtom.getAtomicNum() != mAtom.getAtomicNum()) {
    return false;
  }
  if (qAtom.getAtomicNum() == 0) {
    // dummy: only enforce isotope match when both sides specify one
    const int qIso = qAtom.getIsotope();
    const int mIso = mAtom.getIsotope();
    if (qIso && mIso && qIso != mIso) {
      return false;
    }
    return true;
  }
  if (qAtom.getFormalCharge() &&
      qAtom.getFormalCharge() != mAtom.getFormalCharge()) {
    return false;
  }
  if (qAtom.getIsotope() && qAtom.getIsotope() != mAtom.getIsotope()) {
    return false;
  }
  if (qAtom.getNumRadicalElectrons() &&
      qAtom.getNumRadicalElectrons() != mAtom.getNumRadicalElectrons()) {
    return false;
  }
  return true;
}

//! Compares the atom-property entries listed in `properties` between two atoms
//! by index. Mirrors the legacy propertyCompat<T> behavior.
bool atomPropertyCompat(const RDMol &qmol, atomindex_t qIdx, const RDMol &mmol,
                        atomindex_t mIdx,
                        const std::vector<std::string> &properties) {
  for (const auto &prop : properties) {
    PropToken token(prop);
    std::string qProp;
    std::string mProp;
    const bool qHas = qmol.getAtomPropIfPresent(token, qIdx, qProp);
    const bool mHas = mmol.getAtomPropIfPresent(token, mIdx, mProp);
    if (qHas && mHas) {
      if (qProp != mProp) {
        return false;
      }
    } else if (qHas || mHas) {
      return false;
    }
  }
  return true;
}

bool bondPropertyCompat(const RDMol &qmol, atomindex_t qIdx, const RDMol &mmol,
                        atomindex_t mIdx,
                        const std::vector<std::string> &properties) {
  for (const auto &prop : properties) {
    PropToken token(prop);
    std::string qProp;
    std::string mProp;
    const bool qHas = qmol.getBondPropIfPresent(token, qIdx, qProp);
    const bool mHas = mmol.getBondPropIfPresent(token, mIdx, mProp);
    if (qHas && mHas) {
      if (qProp != mProp) {
        return false;
      }
    } else if (qHas || mHas) {
      return false;
    }
  }
  return true;
}

}  // namespace

namespace detail {
// Helper class used by the sortMatchesByDegreeOfCoreSubstitution
// and getMostSubstitutedCoreMatch functions. A penalty of 1.0 is assigned
// to matches for each terminal dummy atom matching hydrogen.
// To make the sort stable in case of ties, a fraction of 1.0
// is added to each score based on match indices.
class ScoreMatchesByDegreeOfCoreSubstitution {
 public:
  typedef std::pair<unsigned int, double> IdxScorePair;
  ScoreMatchesByDegreeOfCoreSubstitution(
      const RDKit::ROMol &mol, const RDKit::ROMol &query,
      const std::vector<RDKit::MatchVectType> &matches)
      : d_mol(mol),
        d_query(query),
        d_matches(matches),
        d_sumIndices(0.0),
        d_minIdx(-1),
        d_isSorted(false) {
    PRECONDITION(!matches.empty(), "matches must not be empty");
    auto na = d_mol.getNumAtoms();
    d_sumIndices = static_cast<double>(na * (na + 1) / 2);
    unsigned int i = 0;
    d_matchIdxVsScore.reserve(d_matches.size());
    for (const auto &match : d_matches) {
      d_matchIdxVsScore.emplace_back(i++, computeScore(match));
    }
  }
  const RDKit::MatchVectType &getMostSubstitutedCoreMatch() {
    if (d_minIdx == -1) {
      d_minIdx = std::min_element(d_matchIdxVsScore.begin(),
                                  d_matchIdxVsScore.end(), compare)
                     ->first;
    }
    return d_matches.at(d_minIdx);
  }
  std::vector<MatchVectType> sortMatchesByDegreeOfCoreSubstitution() {
    if (!d_isSorted) {
      std::sort(d_matchIdxVsScore.begin(), d_matchIdxVsScore.end(), compare);
      d_isSorted = true;
      d_minIdx = d_matchIdxVsScore.front().first;
    }
    std::vector<MatchVectType> res(d_matches.size());
    std::transform(
        d_matchIdxVsScore.begin(), d_matchIdxVsScore.end(), res.begin(),
        [this](const IdxScorePair &pair) { return d_matches.at(pair.first); });
    return res;
  }

 private:
  static bool compare(const IdxScorePair &aPair, const IdxScorePair &bPair) {
    return (aPair.second < bPair.second);
  }
  bool doesRGroupMatchHydrogen(const std::pair<int, int> &pair) const {
    const auto queryAtom = d_query.getAtomWithIdx(pair.first);
    const auto molAtom = d_mol.getAtomWithIdx(pair.second);
    return (molAtom->getAtomicNum() == 1 &&
            isAtomTerminalRGroupOrQueryHydrogen(queryAtom));
  }
  double computeScore(const RDKit::MatchVectType &match) const {
    double penalty = 0.0;
    double i = 0.0;
    for (const auto &pair : match) {
      i += static_cast<double>(pair.second);
      if (doesRGroupMatchHydrogen(pair)) {
        penalty += 1.0;
      }
    }
    penalty += i / d_sumIndices;
    return penalty;
  }
  const RDKit::ROMol &d_mol;
  const RDKit::ROMol &d_query;
  const std::vector<RDKit::MatchVectType> &d_matches;
  std::vector<IdxScorePair> d_matchIdxVsScore;
  double d_sumIndices;
  int d_minIdx;
  bool d_isSorted;
};
}  // namespace detail

bool atomCompat(const RDMol &qmol, atomindex_t qIdx, const RDMol &mmol,
                atomindex_t mIdx, const SubstructMatchParameters &ps) {
  if (ps.extraAtomCheckOverridesDefaultCheck && ps.extraAtomCheck) {
    return ps.extraAtomCheck(qmol, qIdx, mmol, mIdx);
  }

  const bool qHasQuery = qmol.hasAtomQuery(qIdx);
  const bool mHasQuery = mmol.hasAtomQuery(mIdx);

  bool res = false;
  if (ps.useQueryQueryMatches && qHasQuery && mHasQuery) {
    // Query-query match path: structural compare on the legacy query tree.
    // queriesMatch is keyed on the Atom*-typed query form, so we bridge the
    // two QueryAtom wrappers; this path is rarely exercised.
    const Atom *qAtomPtr = qmol.asROMol().getAtomWithIdx(qIdx);
    const Atom *mAtomPtr = mmol.asROMol().getAtomWithIdx(mIdx);
    res = static_cast<const QueryAtom *>(qAtomPtr)->QueryMatch(
        static_cast<const QueryAtom *>(mAtomPtr));
  } else if (qHasQuery) {
    res = qmol.getAtomQuery(qIdx)->Match(ConstRDMolAtom{&mmol, mIdx});
  } else {
    res = atomDataMatch(qmol.getAtom(qIdx), mmol.getAtom(mIdx));
  }
  if (!res) {
    return false;
  }
  if (!ps.atomProperties.empty()) {
    if (!atomPropertyCompat(qmol, qIdx, mmol, mIdx, ps.atomProperties)) {
      return false;
    }
  }
  if (ps.extraAtomCheck && !ps.extraAtomCheck(qmol, qIdx, mmol, mIdx)) {
    return false;
  }

  return res;
}

bool bondCompat(const RDMol &qmol, atomindex_t qIdx, const RDMol &mmol,
                atomindex_t mIdx, const SubstructMatchParameters &ps) {
  if (ps.extraBondCheckOverridesDefaultCheck && ps.extraBondCheck) {
    return ps.extraBondCheck(qmol, qIdx, mmol, mIdx);
  }

  const BondData &qBond = qmol.getBond(qIdx);
  const BondData &mBond = mmol.getBond(mIdx);
  const bool qHasQuery = qmol.hasBondQuery(qIdx);
  const bool mHasQuery = mmol.hasBondQuery(mIdx);

  using BondEnums::BondType;
  auto isConjugatedSingleOrDouble = [&](const BondData &bond) {
    return bond.getIsConjugated() && (bond.getBondType() == BondType::SINGLE ||
                                      bond.getBondType() == BondType::DOUBLE);
  };
  auto isSingleOrDouble = [&](const BondData &bond) {
    return bond.getBondType() == BondType::SINGLE ||
           bond.getBondType() == BondType::DOUBLE;
  };

  bool res = false;
  if (ps.useQueryQueryMatches && qHasQuery && mHasQuery) {
    const Bond *qBondPtr = qmol.asROMol().getBondWithIdx(qIdx);
    const Bond *mBondPtr = mmol.asROMol().getBondWithIdx(mIdx);
    res = static_cast<const QueryBond *>(qBondPtr)->QueryMatch(
        static_cast<const QueryBond *>(mBondPtr));
  } else if (ps.aromaticMatchesConjugated && !qHasQuery && !mHasQuery &&
             ((qBond.getBondType() == BondType::AROMATIC &&
               mBond.getBondType() == BondType::AROMATIC) ||
              (qBond.getBondType() == BondType::AROMATIC &&
               isConjugatedSingleOrDouble(mBond)) ||
              (mBond.getBondType() == BondType::AROMATIC &&
               isConjugatedSingleOrDouble(qBond)))) {
    res = true;
  } else if (ps.aromaticMatchesSingleOrDouble && !qHasQuery && !mHasQuery &&
             ((qBond.getBondType() == BondType::AROMATIC &&
               mBond.getBondType() == BondType::AROMATIC) ||
              (qBond.getBondType() == BondType::AROMATIC &&
               isSingleOrDouble(mBond)) ||
              (mBond.getBondType() == BondType::AROMATIC &&
               isSingleOrDouble(qBond)))) {
    res = true;
  } else if (qHasQuery) {
    res = qmol.getBondQuery(qIdx)->Match(ConstRDMolBond{&mmol, mIdx});
  } else if (qBond.getBondType() == BondType::UNSPECIFIED ||
             mBond.getBondType() == BondType::UNSPECIFIED) {
    res = true;
  } else {
    res = qBond.getBondType() == mBond.getBondType();
  }
  if (!res) {
    return false;
  }
  if (qBond.getBondType() == BondType::DATIVE &&
      mBond.getBondType() == BondType::DATIVE) {
    // For dative bonds the direction must also match.
    if (!atomDataMatch(qmol.getAtom(qBond.getBeginAtomIdx()),
                       mmol.getAtom(mBond.getBeginAtomIdx())) ||
        !atomDataMatch(qmol.getAtom(qBond.getEndAtomIdx()),
                       mmol.getAtom(mBond.getEndAtomIdx()))) {
      return false;
    }
  }
  if (!ps.bondProperties.empty()) {
    if (!bondPropertyCompat(qmol, qIdx, mmol, mIdx, ps.bondProperties)) {
      return false;
    }
  }
  if (ps.extraBondCheck && !ps.extraBondCheck(qmol, qIdx, mmol, mIdx)) {
    return false;
  }

  return res;
}

void removeDuplicates(std::vector<MatchVectType> &matches,
                      unsigned int nAtoms) {
  //
  //  This works by tracking the indices of the atoms in each match vector.
  //  This can lead to unexpected behavior when looking at rings and queries
  //  that don't specify bond orders.  For example querying this molecule:
  //    C1CCC=1
  //  with the pattern constructed from SMARTS C~C~C~C will return a
  //  single match, despite the fact that there are 4 different paths
  //  when valence is considered.  The defense of this behavior is
  //  that the 4 paths are equivalent in the semantics of the query.
  //  Also, OELib returns the same results
  //
  std::unordered_set<std::string> seen;
  std::vector<MatchVectType> res;
  res.reserve(matches.size());
  seen.reserve(matches.size());
  for (const auto &match : matches) {
    std::string val(nAtoms, '0');
    for (const auto &ci : match) {
      val[ci.second] = '1';
    }
    const bool inserted = seen.insert(std::move(val)).second;
    if (inserted) {
      res.push_back(match);
    }
  }
  res.shrink_to_fit();
  matches = std::move(res);
}

const MatchVectType &getMostSubstitutedCoreMatch(
    const ROMol &mol, const ROMol &core,
    const std::vector<MatchVectType> &matches) {
  detail::ScoreMatchesByDegreeOfCoreSubstitution matchScorer(mol, core,
                                                             matches);
  return matchScorer.getMostSubstitutedCoreMatch();
}

std::vector<MatchVectType> sortMatchesByDegreeOfCoreSubstitution(
    const ROMol &mol, const ROMol &core,
    const std::vector<MatchVectType> &matches) {
  detail::ScoreMatchesByDegreeOfCoreSubstitution matchScorer(mol, core,
                                                             matches);
  return matchScorer.sortMatchesByDegreeOfCoreSubstitution();
}

bool isAtomTerminalRGroupOrQueryHydrogen(const Atom *atom) {
  return (atom->getDegree() == 1 && isAtomDummy(atom)) ||
         (atom->hasQuery() &&
          describeQuery(atom).find("AtomAtomicNum 1 = val") !=
              std::string::npos);
}

#define PT_OPT_GET(opt) params.opt = pt.get(#opt, params.opt)
#define PT_OPT_PUT(opt) pt.put(#opt, params.opt);

void updateSubstructMatchParamsFromJSON(SubstructMatchParameters &params,
                                        const std::string &json) {
  if (json.empty()) {
    return;
  }
  std::istringstream ss;
  ss.str(json);
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  PT_OPT_GET(useChirality);
  PT_OPT_GET(useEnhancedStereo);
  PT_OPT_GET(aromaticMatchesConjugated);
  PT_OPT_GET(useQueryQueryMatches);
  PT_OPT_GET(recursionPossible);
  PT_OPT_GET(uniquify);
  PT_OPT_GET(maxMatches);
  PT_OPT_GET(maxRecursiveMatches);
  PT_OPT_GET(numThreads);
  PT_OPT_GET(specifiedStereoQueryMatchesUnspecified);
  PT_OPT_GET(aromaticMatchesSingleOrDouble);
}

std::string substructMatchParamsToJSON(const SubstructMatchParameters &params) {
  boost::property_tree::ptree pt;

  PT_OPT_PUT(useChirality);
  PT_OPT_PUT(useEnhancedStereo);
  PT_OPT_PUT(aromaticMatchesConjugated);
  PT_OPT_PUT(useQueryQueryMatches);
  PT_OPT_PUT(recursionPossible);
  PT_OPT_PUT(uniquify);
  PT_OPT_PUT(maxMatches);
  PT_OPT_PUT(maxRecursiveMatches);
  PT_OPT_PUT(numThreads);
  PT_OPT_PUT(specifiedStereoQueryMatchesUnspecified);
  PT_OPT_PUT(aromaticMatchesSingleOrDouble);

  std::stringstream ss;
  boost::property_tree::json_parser::write_json(ss, pt);
  return ss.str();
}

}  // namespace RDKit
