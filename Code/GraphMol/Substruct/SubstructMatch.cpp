//
//  Copyright (C) 2001-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDThreads.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Resonance.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/Chirality.h>

#include "SubstructMatch.h"
#include "SubstructUtils.h"
#include <GraphMol/GenericGroups/GenericGroups.h>
#include <boost/smart_ptr.hpp>
#include <map>
#include <span>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#include <thread>
#include <future>
#endif

#include "vf2.hpp"

namespace RDKit {
namespace detail {

namespace {
bool hasChiralLabel(const AtomData &atom) {
  return atom.getChiralTag() == AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW ||
         atom.getChiralTag() == AtomEnums::ChiralType::CHI_TETRAHEDRAL_CCW;
}

bool enhancedStereoIsOK(
    const RDMol &mol, const RDMol &query,
    std::unordered_map<unsigned int, unsigned int> &q_to_mol,
    const std::unordered_map<unsigned int, StereoGroup const *>
        &molStereoGroups,
    const std::unordered_map<unsigned int, bool> &matches) {
  std::unordered_map<unsigned int, StereoGroup const *> molAtomsToQueryGroups;

  // If the query has stereo groups:
  // * OR only matches AND or OR (not absolute)
  // * AND only matches OR
  for (const auto &sg : query.asROMol().getStereoGroups()) {
    if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
      continue;
    }
    const bool is_and = sg.getGroupType() == StereoGroupType::STEREO_AND;
    for (const atomindex_t qAtomIdx : sg.getAtomIndices()) {
      const auto mol_group = molStereoGroups.find(q_to_mol[qAtomIdx]);
      if (mol_group == molStereoGroups.end()) {
        // group matching absolute. not ok.
        return false;
      } else if (is_and && mol_group->second->getGroupType() !=
                               StereoGroupType::STEREO_AND) {
        // AND matching OR. not ok.
        return false;
      }

      molAtomsToQueryGroups[q_to_mol[qAtomIdx]] = &sg;
    }
  }

  // If the mol has stereo groups:
  // * All atoms must either be the same or opposite, you can't mix
  // * Only one stereogroup must cover all matched atoms in the mol stereo group
  for (const auto &sg : mol.asROMol().getStereoGroups()) {
    if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
      continue;
    }
    bool doesMatch = false;
    bool seen = false;
    StereoGroup const *QGroup = nullptr;

    for (const atomindex_t mAtomIdx : sg.getAtomIndices()) {
      auto thisDoesMatch = matches.find(mAtomIdx);
      if (thisDoesMatch == matches.end()) {
        continue;
      }

      auto pos = molAtomsToQueryGroups.find(mAtomIdx);
      auto thisQGroup =
          pos == molAtomsToQueryGroups.end() ? nullptr : pos->second;
      if (!seen) {
        doesMatch = thisDoesMatch->second;
        QGroup = thisQGroup;
        seen = true;
      } else if (doesMatch != thisDoesMatch->second) {
        // diastereomer. not ok.
        return false;
      } else if (thisQGroup != QGroup) {
        // mix of groups in query. not ok.
        return false;
      }
    }
  }

  return true;
}

}  // namespace

typedef std::map<unsigned int, QueryAtom::QUERYATOM_QUERY *> SUBQUERY_MAP;

typedef struct {
  ResonanceMolSupplier &resMolSupplier;
  const RDMol &query;
  const SubstructMatchParameters &params;
} ResSubstructMatchHelperArgs_;

void MatchSubqueries(const RDMol &mol, QueryAtom::QUERYATOM_QUERY *q,
                     const SubstructMatchParameters &params,
                     SUBQUERY_MAP &subqueryMap,
                     std::vector<RecursiveStructureQuery *> &locked);

bool insertIfNeeded(std::set<MatchVectType> &matches, const MatchVectType &m) {
  bool shouldInsert = true;
  std::unordered_set<int> matchAsSet;
  std::transform(m.begin(), m.end(),
                 std::inserter(matchAsSet, matchAsSet.begin()),
                 [](const std::pair<int, int> &p) { return p.second; });
  for (auto it = matches.begin(); it != matches.end(); ++it) {
    std::unordered_set<int> existingMatchAsSet;
    std::transform(
        it->begin(), it->end(),
        std::inserter(existingMatchAsSet, existingMatchAsSet.begin()),
        [](const std::pair<int, int> &p) { return p.second; });
    if (matchAsSet == existingMatchAsSet) {
      if (m < *it) {
        matches.erase(it);
      } else {
        shouldInsert = false;
      }
      break;
    }
  }
  if (shouldInsert) {
    matches.insert(m);
  }
  return shouldInsert;
}

bool tryToInsert(std::set<MatchVectType> &matches, const MatchVectType &match,
                 const SubstructMatchParameters &params) {
  if (matches.size() == params.maxMatches) {
    return false;
  }
  if (!params.uniquify) {
    matches.insert(match);
  } else {
    insertIfNeeded(matches, match);
  }
  return true;
}

void ResSubstructMatchHelper_(const ResSubstructMatchHelperArgs_ &args,
                              std::set<MatchVectType> *matches, unsigned int bi,
                              unsigned int ei);

typedef std::vector<std::pair<std::uint32_t, std::uint32_t>> ssPairType;

}  // namespace detail

MolMatchFinalCheckFunctor::MolMatchFinalCheckFunctor(
    const RDMol &query, const RDMol &mol, const SubstructMatchParameters &ps)
    : d_query(query), d_mol(mol), d_params(ps) {
  if (d_params.useEnhancedStereo) {
    for (const auto &sg : d_mol.asROMol().getStereoGroups()) {
      if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
        continue;
      }
      for (const atomindex_t aIdx : sg.getAtomIndices()) {
        d_molStereoGroups[aIdx] = &sg;
      }
    }
  }
}

bool MolMatchFinalCheckFunctor::operator()(const std::uint32_t q_c[],
                                           const std::uint32_t m_c[]) {
  if (d_params.extraFinalCheck || d_params.useGenericMatchers) {
    const std::span<const std::uint32_t> aids(m_c, d_query.getNumAtoms());
    if (d_params.useGenericMatchers &&
        !GenericGroups::genericAtomMatcher(d_mol, d_query, aids)) {
      return false;
    }
    if (d_params.extraFinalCheck && !d_params.extraFinalCheck(d_mol, aids)) {
      return false;
    }
  }

  HashedStorageType match;
  if (d_params.uniquify) {
    match.resize(d_mol.getNumAtoms());
#ifdef RDK_INTERNAL_BITSET_HAS_HASH
    match.reset();
#else
    std::fill(match.begin(), match.end(), 0);
#endif
    for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
      match[m_c[i]] = 1;
    }
    if (matchesSeen.find(match) != matchesSeen.end()) {
      return false;
    }
  }

  if (!d_params.useChirality) {
    if (d_params.uniquify) {
      matchesSeen.insert(match);
    }
    return true;
  }

  std::unordered_map<unsigned int, bool> matches;

  // check chiral atoms:
  const std::uint32_t qNumAtoms = d_query.getNumAtoms();
  for (std::uint32_t i = 0; i < qNumAtoms; ++i) {
    const AtomData &qAt = d_query.getAtom(q_c[i]);
    const std::uint32_t qDegree = d_query.getAtomDegree(q_c[i]);

    // With less than 3 neighbors we can't establish CW/CCW parity, so query
    // will be a match if it has any kind of chirality.
    if (qDegree < 3 || !detail::hasChiralLabel(qAt)) {
      continue;
    }
    const AtomData &mAt = d_mol.getAtom(m_c[i]);
    const std::uint32_t mDegree = d_mol.getAtomDegree(m_c[i]);
    if (!detail::hasChiralLabel(mAt)) {
      if (d_params.specifiedStereoQueryMatchesUnspecified) {
        continue;
      }
      return false;
    }
    if (qDegree > mDegree) {
      return false;
    }

    INT_LIST qOrder;
    INT_LIST mOrder;
    for (std::uint32_t j = 0; j < qNumAtoms; ++j) {
      const std::uint32_t qBondIdx =
          d_query.getBondIndexBetweenAtoms(q_c[i], q_c[j]);
      const std::uint32_t mBondIdx =
          d_mol.getBondIndexBetweenAtoms(m_c[i], m_c[j]);
      if (qBondIdx != std::numeric_limits<std::uint32_t>::max() &&
          mBondIdx != std::numeric_limits<std::uint32_t>::max()) {
        mOrder.push_back(static_cast<int>(mBondIdx));
        qOrder.push_back(static_cast<int>(qBondIdx));
        if (mOrder.size() == qDegree) {
          break;
        }
      }
    }
    CHECK_INVARIANT(qOrder.size() == qDegree, "missing matches");
    CHECK_INVARIANT(qOrder.size() == mOrder.size(), "bad matches");
    const int qPermCount =
        d_query.asROMol().getAtomWithIdx(q_c[i])->getPerturbationOrder(qOrder);

    const unsigned unmatchedNeighbors = mDegree - mOrder.size();
    mOrder.insert(mOrder.end(), unmatchedNeighbors, -1);

    INT_LIST moOrder;
    auto [mBondsBegin, mBondsEnd] = d_mol.getAtomBonds(m_c[i]);
    for (auto bondIt = mBondsBegin; bondIt != mBondsEnd; ++bondIt) {
      const int dbidx = static_cast<int>(*bondIt);
      if (std::find(mOrder.begin(), mOrder.end(), dbidx) != mOrder.end()) {
        moOrder.push_back(dbidx);
      } else {
        moOrder.push_back(-1);
      }
    }

    const int mPermCount =
        static_cast<int>(countSwapsToInterconvert(moOrder, mOrder));

    const bool requireMatch = qPermCount % 2 == mPermCount % 2;
    const bool labelsMatch = qAt.getChiralTag() == mAt.getChiralTag();
    const bool matchOK = requireMatch == labelsMatch;

    // if this is not part of a stereogroup and doesn't match, return false
    const auto msg = d_molStereoGroups.find(m_c[i]);
    if (msg == d_molStereoGroups.end()) {
      if (!matchOK) {
        return false;
      }
    } else {
      matches[m_c[i]] = matchOK;
    }
  }

  std::unordered_map<unsigned int, unsigned int> q_to_mol;
  for (std::uint32_t j = 0; j < qNumAtoms; ++j) {
    q_to_mol[q_c[j]] = m_c[j];
  }

  if (d_params.useEnhancedStereo) {
    if (!detail::enhancedStereoIsOK(d_mol, d_query, q_to_mol, d_molStereoGroups,
                                    matches)) {
      return false;
    }
  }

  // now check double bonds
  using BondEnums::BondType;
  using BondEnums::BondStereo;
  const std::uint32_t qNumBonds = d_query.getNumBonds();
  for (std::uint32_t qBondIdx = 0; qBondIdx < qNumBonds; ++qBondIdx) {
    const BondData &qBnd = d_query.getBond(qBondIdx);
    if (qBnd.getBondType() != BondType::DOUBLE ||
        qBnd.getStereo() <= BondStereo::STEREOANY) {
      continue;
    }

    if (!d_query.hasBondStereoAtoms(qBondIdx)) {
      continue;
    }
    const atomindex_t *qStereoAtoms = d_query.getBondStereoAtoms(qBondIdx);
    if (qStereoAtoms[0] == atomindex_t(-1) ||
        qStereoAtoms[1] == atomindex_t(-1)) {
      continue;
    }

    const std::uint32_t mBondIdx = d_mol.getBondIndexBetweenAtoms(
        q_to_mol[qBnd.getBeginAtomIdx()], q_to_mol[qBnd.getEndAtomIdx()]);
    CHECK_INVARIANT(mBondIdx != std::numeric_limits<std::uint32_t>::max(),
                    "Matching bond not found");
    const BondData &mBnd = d_mol.getBond(mBondIdx);
    if (mBnd.getBondType() != BondType::DOUBLE) {
      continue;
    }

    if (!d_params.specifiedStereoQueryMatchesUnspecified &&
        mBnd.getStereo() <= BondStereo::STEREOANY) {
      return false;
    }

    if (!d_mol.hasBondStereoAtoms(mBondIdx)) {
      continue;
    }
    const atomindex_t *mStereoAtoms = d_mol.getBondStereoAtoms(mBondIdx);
    if (mStereoAtoms[0] == atomindex_t(-1) ||
        mStereoAtoms[1] == atomindex_t(-1)) {
      continue;
    }

    unsigned int end1Matches = 0;
    unsigned int end2Matches = 0;
    if (q_to_mol[qBnd.getBeginAtomIdx()] == mBnd.getBeginAtomIdx()) {
      // query Begin == mol Begin
      if (q_to_mol[qStereoAtoms[0]] == mStereoAtoms[0]) {
        end1Matches = 1;
      }
      if (q_to_mol[qStereoAtoms[1]] == mStereoAtoms[1]) {
        end2Matches = 1;
      }
    } else {
      // query End == mol Begin
      if (q_to_mol[qStereoAtoms[0]] == mStereoAtoms[1]) {
        end1Matches = 1;
      }
      if (q_to_mol[qStereoAtoms[1]] == mStereoAtoms[0]) {
        end2Matches = 1;
      }
    }

    const unsigned totalMatches = end1Matches + end2Matches;
    const auto mStereo = Chirality::translateEZLabelToCisTrans(
        static_cast<Bond::BondStereo>(mBnd.getStereo()));
    const auto qStereo = Chirality::translateEZLabelToCisTrans(
        static_cast<Bond::BondStereo>(qBnd.getStereo()));

    if (mStereo == qStereo && totalMatches == 1) {
      return false;
    }
    if (mStereo != qStereo && totalMatches != 1) {
      return false;
    }
  }
  if (d_params.uniquify) {
    matchesSeen.insert(match);
  }
  return true;
}

namespace detail {

class AtomLabelFunctor {
 public:
  AtomLabelFunctor(const RDMol &query, const RDMol &mol,
                   const SubstructMatchParameters &ps)
      : d_query(query), d_mol(mol), d_params(ps) {};

  bool operator()(unsigned int i, unsigned int j) const {
    if (d_params.useChirality) {
      const AtomData &qAt = d_query.getAtom(i);
      if (qAt.getChiralTag() == AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW ||
          qAt.getChiralTag() == AtomEnums::ChiralType::CHI_TETRAHEDRAL_CCW) {
        const AtomData &mAt = d_mol.getAtom(j);
        if (!d_params.specifiedStereoQueryMatchesUnspecified &&
            mAt.getChiralTag() != AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW &&
            mAt.getChiralTag() != AtomEnums::ChiralType::CHI_TETRAHEDRAL_CCW) {
          return false;
        }
      }
    }
    return atomCompat(d_query, i, d_mol, j, d_params);
  }

 private:
  const RDMol &d_query;
  const RDMol &d_mol;
  const SubstructMatchParameters &d_params;
};
class BondLabelFunctor {
 public:
  BondLabelFunctor(const RDMol &query, const RDMol &mol,
                   const SubstructMatchParameters &ps)
      : d_query(query), d_mol(mol), d_params(ps) {};
  bool operator()(::RDKit::RDMol::edge_descriptor i,
                  ::RDKit::RDMol::edge_descriptor j) const {
    using BondEnums::BondType;
    using BondEnums::BondStereo;
    if (d_params.useChirality) {
      const BondData &qBnd = d_query.getBond(*i);
      if (qBnd.getBondType() == BondType::DOUBLE &&
          qBnd.getStereo() > BondStereo::STEREOANY) {
        const BondData &mBnd = d_mol.getBond(*j);
        if (mBnd.getBondType() == BondType::DOUBLE &&
            !d_params.specifiedStereoQueryMatchesUnspecified &&
            mBnd.getStereo() <= BondStereo::STEREOANY) {
          return false;
        }
      }
    }
    return bondCompat(d_query, *i, d_mol, *j, d_params);
  }

 private:
  const RDMol &d_query;
  const RDMol &d_mol;
  const SubstructMatchParameters &d_params;
};
void ResSubstructMatchHelper_(const ResSubstructMatchHelperArgs_ &args,
                              std::set<MatchVectType> *matches, unsigned int bi,
                              unsigned int ei) {
  for (unsigned int i = bi;
       (matches->size() < args.params.maxMatches) && (i < ei); ++i) {
    std::unique_ptr<ROMol> mol{args.resMolSupplier[i]};
    std::vector<MatchVectType> matchesTmp =
        SubstructMatch(mol->asRDMol(), args.query, args.params);
    for (const auto &match : matchesTmp) {
      if (!tryToInsert(*matches, match, args.params)) {
        break;
      }
    }
  }
};

struct RecursiveLocker {
  std::vector<RecursiveStructureQuery *> locked;
  RecursiveLocker(const RDMol &query, const bool recursionPossible) {
    if (recursionPossible) {
      locked.reserve(query.getNumAtoms());
    }
  }

  ~RecursiveLocker() {
    for (auto v : locked) {
      v->clear();
#ifdef RDK_BUILD_THREADSAFE_SSS
      v->d_mutex.unlock();
#endif
    }
  }
};

// A minimal container which satisfies the vf2_all() output-sequence interface
// but only counts matches instead of storing them.
struct MatchCounter {
  using value_type = ssPairType;

  void clear() { d_count = 0; }
  void resize(size_t) { d_count = 0; }
  void reserve(size_t) {}

  bool empty() const { return d_count == 0; }
  size_t size() const { return d_count; }

  void push_back(const value_type &) { ++d_count; }

 private:
  size_t d_count = 0;
};
}  // namespace detail

// ----------------------------------------------
//
// find all matches
std::vector<MatchVectType> SubstructMatch(
    const RDMol &mol, const RDMol &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> matches;
  const std::uint32_t mNumAtoms = mol.getNumAtoms();
  const std::uint32_t qNumAtoms = query.getNumAtoms();
  if (!mNumAtoms || !qNumAtoms || qNumAtoms > mNumAtoms) {
    return matches;
  }

  detail::RecursiveLocker locker(query, params.recursionPossible);

  if (params.recursionPossible) {
    detail::SUBQUERY_MAP subqueryMap;
    const ROMol &qROMol = query.asROMol();
    for (const auto atom : qROMol.atoms()) {
      if (atom->hasQuery()) {
        detail::MatchSubqueries(mol, atom->getQuery(), params, subqueryMap,
                                locker.locked);
      }
    }
  }

  detail::AtomLabelFunctor atomLabeler(query, mol, params);
  detail::BondLabelFunctor bondLabeler(query, mol, params);
  MolMatchFinalCheckFunctor matchChecker(query, mol, params);

  std::vector<detail::ssPairType> pms;
  bool found =
      boost::vf2_all(query.getTopology(), mol.getTopology(), atomLabeler,
                     bondLabeler, matchChecker, pms, params.maxMatches);
  if (found) {
    matches.reserve(pms.size());
    MatchVectType matchVect(qNumAtoms);
    for (const auto &pairs : pms) {
      for (const auto &pair : pairs) {
        matchVect[pair.first] =
            std::make_pair(static_cast<int>(pair.first),
                           static_cast<int>(pair.second));
      }
      matches.push_back(matchVect);
    }
  }
  return matches;
}

std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const ROMol &query,
    const SubstructMatchParameters &params) {
  return SubstructMatch(mol.asRDMol(), query.asRDMol(), params);
}

unsigned int SubstructMatchCount(const RDMol &mol, const RDMol &query,
                                 const SubstructMatchParameters &params) {
  if (!mol.getNumAtoms() || !query.getNumAtoms()) {
    return 0;
  }

  detail::RecursiveLocker locker(query, params.recursionPossible);

  if (params.recursionPossible) {
    detail::SUBQUERY_MAP subqueryMap;
    const ROMol &qROMol = query.asROMol();
    for (const auto atom : qROMol.atoms()) {
      if (atom->hasQuery()) {
        detail::MatchSubqueries(mol, atom->getQuery(), params, subqueryMap,
                                locker.locked);
      }
    }
  }

  detail::AtomLabelFunctor atomLabeler(query, mol, params);
  detail::BondLabelFunctor bondLabeler(query, mol, params);
  MolMatchFinalCheckFunctor matchChecker(query, mol, params);

  detail::MatchCounter counter;
  boost::vf2_all(query.getTopology(), mol.getTopology(), atomLabeler,
                 bondLabeler, matchChecker, counter, params.maxMatches);
  return static_cast<unsigned int>(counter.size());
}

unsigned int SubstructMatchCount(const ROMol &mol, const ROMol &query,
                                 const SubstructMatchParameters &params) {
  return SubstructMatchCount(mol.asRDMol(), query.asRDMol(), params);
}

std::vector<MatchVectType> SubstructMatch(
    const MolBundle &bundle, const ROMol &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> res;
  for (unsigned int i = 0; i < bundle.size() && res.empty(); ++i) {
    res = SubstructMatch(*bundle[i], query, params);
  }
  return res;
}

std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const MolBundle &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> res;
  for (unsigned int i = 0; i < query.size() && res.empty(); ++i) {
    res = SubstructMatch(mol, *query[i], params);
  }
  return res;
}

std::vector<MatchVectType> SubstructMatch(
    const MolBundle &mol, const MolBundle &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> res;
  for (unsigned int i = 0; i < mol.size() && res.empty(); ++i) {
    for (unsigned int j = 0; j < query.size() && res.empty(); ++j) {
      res = SubstructMatch(*mol[i], *query[j], params);
    }
  }
  return res;
}

// ----------------------------------------------
//
// find all matches in a ResonanceMolSupplier object
//
//
std::vector<MatchVectType> SubstructMatch(
    ResonanceMolSupplier &resMolSupplier, const ROMol &query,
    const SubstructMatchParameters &params) {
  std::set<MatchVectType> matches;
  detail::ResSubstructMatchHelperArgs_ args = {resMolSupplier, query.asRDMol(),
                                               params};
  unsigned int nt =
      std::min(resMolSupplier.length(), getNumThreadsToUse(params.numThreads));
  if (nt == 1) {
    detail::ResSubstructMatchHelper_(args, &matches, 0,
                                     resMolSupplier.length());
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    std::vector<std::future<void>> tg;
    std::vector<std::unique_ptr<std::set<MatchVectType>>> matchesThread(nt);
    unsigned int ei = 0;
    double dpt =
        static_cast<double>(resMolSupplier.length()) / static_cast<double>(nt);
    double dc = 0.0;
    for (unsigned int ti = 0; ti < nt; ++ti) {
      matchesThread[ti] = std::make_unique<std::set<MatchVectType>>();
      unsigned int bi = ei;
      dc += dpt;
      ei = static_cast<unsigned int>(floor(dc));
      tg.emplace_back(std::async(std::launch::async,
                                 detail::ResSubstructMatchHelper_, args,
                                 matchesThread[ti].get(), bi, ei));
    }
    for (auto &fut : tg) {
      fut.get();
    }

    for (unsigned int ti = 0; ti < nt; ++ti) {
      for (const auto &match : *matchesThread[ti]) {
        if (!detail::tryToInsert(matches, match, args.params)) {
          break;
        }
      }
    }
  }
#endif
  return std::vector<MatchVectType>(matches.begin(), matches.end());
}

namespace detail {
unsigned int RecursiveMatcher(const RDMol &mol, const RDMol &query,
                              std::vector<int> &matches,
                              SUBQUERY_MAP &subqueryMap,
                              const SubstructMatchParameters &params,
                              std::vector<RecursiveStructureQuery *> &locked) {
  SubstructMatchParameters lparams = params;
  lparams.maxMatches = std::max(params.maxRecursiveMatches, params.maxMatches);
  lparams.uniquify = false;
  const ROMol &qROMol = query.asROMol();
  for (auto qAtom : qROMol.atoms()) {
    if (qAtom->hasQuery()) {
      MatchSubqueries(mol, qAtom->getQuery(), lparams, subqueryMap, locked);
    }
  }

  detail::AtomLabelFunctor atomLabeler(query, mol, lparams);
  detail::BondLabelFunctor bondLabeler(query, mol, lparams);
  MolMatchFinalCheckFunctor matchChecker(query, mol, lparams);

  matches.clear();
  matches.resize(0);
  std::vector<detail::ssPairType> pms;
  bool found =
      boost::vf2_all(query.getTopology(), mol.getTopology(), atomLabeler,
                     bondLabeler, matchChecker, pms, lparams.maxMatches);
  unsigned int res = 0;
  if (found) {
    matches.reserve(pms.size());
    int rootIdx = -1;
    bool hasRootAtom = query.getMolPropIfPresent(
        PropToken(common_properties::_queryRootAtom), rootIdx);
    for (const auto &pairs : pms) {
      if (!hasRootAtom) {
        matches.push_back(pairs.begin()->second);
      } else {
        bool foundRoot = false;
        for (const auto &pairIter : pairs) {
          if (pairIter.first == static_cast<unsigned int>(rootIdx)) {
            matches.push_back(pairIter.second);
            foundRoot = true;
            break;
          }
        }
        if (!foundRoot) {
          BOOST_LOG(rdErrorLog)
              << "no match found for queryRootAtom" << std::endl;
        }
      }
      if (matches.size() == lparams.maxMatches) {
        break;
      }
    }
    res = matches.size();
  }
  return res;
}

void MatchSubqueries(const RDMol &mol, QueryAtom::QUERYATOM_QUERY *query,
                     const SubstructMatchParameters &params,
                     SUBQUERY_MAP &subqueryMap,
                     std::vector<RecursiveStructureQuery *> &locked) {
  PRECONDITION(query, "bad query");
  if (query->getDescription() == "RecursiveStructure") {
    auto *rsq = (RecursiveStructureQuery *)query;
#ifdef RDK_BUILD_THREADSAFE_SSS
    rsq->d_mutex.lock();
#endif
    locked.push_back(rsq);
    rsq->clear();
    bool matchDone = false;
    if (rsq->getSerialNumber() &&
        subqueryMap.find(rsq->getSerialNumber()) != subqueryMap.end()) {
      // we've matched an equivalent serial number before, just
      // copy in the matches:
      matchDone = true;
      auto orsq =
          (const RecursiveStructureQuery *)subqueryMap[rsq->getSerialNumber()];
      for (auto setIter = orsq->beginSet(); setIter != orsq->endSet();
           ++setIter) {
        rsq->insert(*setIter);
      }
    }

    if (!matchDone) {
      ROMol const *queryMol = rsq->getQueryMol();
      if (queryMol) {
        std::vector<int> matchStarts;
        unsigned int res =
            RecursiveMatcher(mol, queryMol->asRDMol(), matchStarts,
                             subqueryMap, params, locked);
        if (res) {
          for (int &matchStart : matchStarts) {
            rsq->insert(matchStart);
          }
        }
      }
      if (rsq->getSerialNumber()) {
        subqueryMap[rsq->getSerialNumber()] = query;
      }
    }
  }

  // now recurse over our children (these things can be nested)
  for (auto childIt = query->beginChildren(); childIt != query->endChildren();
       ++childIt) {
    MatchSubqueries(mol, childIt->get(), params, subqueryMap, locked);
  }
  // std::cout << "<<- back " << (int)query << std::endl;
}

}  // end of namespace detail

bool AtomCoordsMatchFunctor::operator()(const RDMol &queryMol,
                                        atomindex_t queryAtomIdx,
                                        const RDMol &targetMol,
                                        atomindex_t targetAtomIdx) const {
  // Conformers still live on the compat side of the molecule object.
  const ROMol &queryROMol = queryMol.asROMol();
  const ROMol &targetROMol = targetMol.asROMol();
  if (!queryROMol.getNumConformers() || !targetROMol.getNumConformers()) {
    return false;
  }
  const auto &queryPos =
      queryROMol.getConformer(d_queryConfId).getAtomPos(queryAtomIdx);
  const auto &targetPos =
      targetROMol.getConformer(d_refConfId).getAtomPos(targetAtomIdx);
  return (queryPos - targetPos).lengthSq() <= d_tol2;
};

bool AtomCoordsMatchFunctor::operator()(const Atom &queryAtom,
                                        const Atom &targetAtom) const {
  return (*this)(queryAtom.getOwningMol().asRDMol(), queryAtom.getIdx(),
                 targetAtom.getOwningMol().asRDMol(), targetAtom.getIdx());
};

}  // namespace RDKit
