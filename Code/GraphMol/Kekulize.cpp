//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/Rings.h>
#include <GraphMol/SanitException.h>
#include <RDGeneral/RDLog.h>
#include <boost/dynamic_bitset.hpp>
#include <numeric>
#include <utility>

namespace RDKit {
// Local utility namespace
namespace {

void backTrack(RDMol &mol, INT_INT_DEQ_MAP &, int lastOpt, INT_VECT &done,
               INT_DEQUE &aqueue, boost::dynamic_bitset<> &dBndCands,
               boost::dynamic_bitset<> &dBndAdds) {
  // so we made a wrong turn at the lastOpt
  // remove on done list that comes after the lastOpt including itself

  auto ei = std::find(done.begin(), done.end(), lastOpt);
  INT_VECT tdone(done.begin(), ei);

  auto eri = std::find(done.rbegin(), done.rend(), lastOpt);
  ++eri;
  // and push them back onto the stack
  for (auto ri = done.rbegin(); ri != eri; ++ri) {
    aqueue.push_front(*ri);
  }

  // remove any double bonds that were added since we passed through lastOpt
  auto &bondVec = mol.getBondDataVector();
  for (uint32_t bi = 0, nbnds = uint32_t(bondVec.size()); bi < nbnds; ++bi) {
    if (dBndAdds[bi]) {
      BondData &bnd = bondVec[bi];
      const int aid1 = int(bnd.getBeginAtomIdx());
      const int aid2 = int(bnd.getEndAtomIdx());
      // if one of these atoms has been dealt with before lastOpt
      // we don't have to change the double bond addition
      if ((std::find(tdone.begin(), tdone.end(), aid1) == tdone.end()) &&
          (std::find(tdone.begin(), tdone.end(), aid2) == tdone.end())) {
        // otherwise strip the double bond and set it back to single
        // and add the atoms to candidate for double bonds
        dBndAdds[bi] = 0;
        bnd.setBondType(BondEnums::BondType::SINGLE);
        dBndCands[aid1] = 1;
        dBndCands[aid2] = 1;
      }
    }
  }
  done = tdone;
}

void markDbondCands(RDMol &mol, const INT_VECT &allAtms,
                    boost::dynamic_bitset<> &dBndCands, INT_VECT &questions,
                    INT_VECT &done) {
  // ok this function does more than mark atoms that are candidates for
  // double bonds during kekulization
  // - check that a non-aromatic atom does not have any aromatic bonds
  // - marks all aromatic bonds to single bonds
  // - marks atoms that can take a double bond
  const auto &atomVec = mol.getAtomDataVector();
  auto &bondVec = mol.getBondDataVector();
  const RingInfoCache &ringInfo = mol.getRingInfo();

  const bool hasAromaticOrDummyAtom = std::any_of(
      allAtms.begin(), allAtms.end(), [&atomVec, &mol](int allAtm) {
        const AtomData &a = atomVec[allAtm];
        return a.getAtomicNum() == 0 || mol.isAromaticAtom(atomindex_t(allAtm));
      });
  // if there's not at least one atom in the ring that's
  // marked as being aromatic or a dummy,
  // there's no point in continuing:
  if (!hasAromaticOrDummyAtom) {
    return;
  }
  // mark rings which are not candidates for double bonds
  // i.e. that have at least one atom which is in a single ring
  // and is not aromatic
  const uint32_t nRings = ringInfo.numRings();
  boost::dynamic_bitset<> isRingNotCand(nRings);
  for (uint32_t ri = 0; ri < nRings; ++ri) {
    isRingNotCand.set(ri);
    const uint32_t ringBegin = ringInfo.ringBegins[ri];
    const uint32_t ringEnd = ringInfo.ringBegins[ri + 1];
    for (uint32_t pos = ringBegin; pos < ringEnd; ++pos) {
      const uint32_t ai = ringInfo.atomsInRings[pos];
      if (mol.isAromaticAtom(atomindex_t(ai)) &&
          ringInfo.numAtomRings(ai) == 1) {
        isRingNotCand.reset(ri);
        break;
      }
    }
  }
  std::vector<uint32_t> makeSingle;

  boost::dynamic_bitset<> inAllAtms(mol.getNumAtoms());
  for (int allAtm : allAtms) {
    inAllAtms.set(allAtm);
    const atomindex_t atomIdx = atomindex_t(allAtm);
    const AtomData &at = atomVec[atomIdx];

    if (at.getAtomicNum() && !mol.isAromaticAtom(atomIdx)) {
      done.push_back(allAtm);
      continue;
    }

    // count the number of neighbors connected with single,
    // double, or aromatic bonds. Along the way, mark
    // bonds that we will later mark as being single:
    int sbo = 0;
    unsigned nToIgnore = 0;
    unsigned int nonArNonDummyNbr = 0;
    auto [bondBegin, bondEnd] = mol.getAtomBonds(atomIdx);
    for (auto bIt = bondBegin; bIt != bondEnd; ++bIt) {
      const uint32_t bondIdx = *bIt;
      const BondData &bond = bondVec[bondIdx];
      const atomindex_t otherIdx = bond.getOtherAtomIdx(atomIdx);
      const AtomData &otherAt = atomVec[otherIdx];
      if (otherAt.getAtomicNum() && !otherAt.getIsAromatic() &&
          inAllAtms.test(otherIdx)) {
        ++nonArNonDummyNbr;
      }
      if (bond.getIsAromatic() &&
          (bond.getBondType() == BondEnums::BondType::SINGLE ||
           bond.getBondType() == BondEnums::BondType::DOUBLE ||
           bond.getBondType() == BondEnums::BondType::AROMATIC)) {
        ++sbo;
        // Defer setting this bond to single: doing it now could mess with
        // the valence calculation we use below to determine implicit Hs.
        makeSingle.push_back(bondIdx);
      } else {
        const int bondContrib =
            std::lround(bond.getValenceContrib(atomIdx));
        sbo += bondContrib;
        if (!bondContrib) {
          ++nToIgnore;
        }
      }
    }

    const uint32_t numAtomRings = ringInfo.numAtomRings(atomIdx);
    const uint32_t memBegin = ringInfo.atomMembershipBegins[atomIdx];
    const uint32_t memEnd = ringInfo.atomMembershipBegins[atomIdx + 1];
    size_t numNonCandRings = 0;
    for (uint32_t pos = memBegin; pos < memEnd; ++pos) {
      if (isRingNotCand.test(ringInfo.atomMemberships[pos])) {
        ++numNonCandRings;
      }
    }
    if (!at.getAtomicNum() && nonArNonDummyNbr < numAtomRings &&
        numNonCandRings < numAtomRings) {
      // dummies always start as candidates to have a double bond:
      dBndCands[allAtm] = 1;
      // but they don't have to have one, so mark them as questionable:
      questions.push_back(allAtm);
    } else {
      // for non dummies, it's a bit more work to figure out if they
      // can take a double bond:

      sbo += mol.getAtomTotalNumHs(atomIdx);
      auto dv =
          PeriodicTable::getTable()->getDefaultValence(at.getAtomicNum());
      auto chrg = at.getFormalCharge();
      if (isEarlyAtom(at.getAtomicNum())) {
        chrg = -chrg;  // fix for GitHub #65
      }
      // special case for carbon - see GitHub #539
      if (at.getAtomicNum() == 6 && chrg > 0) {
        chrg = -chrg;
      }
      dv += chrg;
      const int tbo =
          int(at.getValence(AtomData::ValenceType::EXPLICIT) +
              at.getValence(AtomData::ValenceType::IMPLICIT));
      const int nRadicals = at.getNumRadicalElectrons();
      const int totalDegree =
          mol.getAtomDegree(atomIdx) +
          at.getValence(AtomData::ValenceType::IMPLICIT) - nToIgnore;

      const auto &valList =
          PeriodicTable::getTable()->getValenceList(at.getAtomicNum());
      unsigned int vi = 1;

      while (tbo > dv && vi < valList.size() && valList[vi] > 0) {
        dv = valList[vi] + chrg;
        ++vi;
      }

      // Kekulize aromatic N-oxides, such as O=n1ccccc1
      // These only reach here if SANITIZE_CLEANUP is disabled.
      if (tbo == 5 && sbo == 4 && dv == 3 && totalDegree == 3 &&
          nRadicals == 0 && chrg == 0 && mol.getAtomTotalNumHs(atomIdx) == 0) {
        switch (at.getAtomicNum()) {
          case 7:   // N
          case 15:  // P
          case 33:  // As
            dv = 5;
            break;
        }
      }
      if (totalDegree + nRadicals >= dv) {
        // if our degree + nRadicals exceeds the default valence,
        // there's no way we can take a double bond, just continue.
        continue;
      }

      // we're a candidate if our total current bond order + nRadicals + 1
      // matches the valence state
      // (including nRadicals here was SF.net issue 3349243)
      if (dv == (sbo + 1 + nRadicals)) {
        dBndCands[allAtm] = 1;
      } else if (!nRadicals && at.getNoImplicit() && dv == (sbo + 2)) {
        // special case: there is currently no radical on the atom, but if
        // if we allow one then this is a candidate:
        dBndCands[allAtm] = 1;
      }
    }
  }  // loop over all atoms in the fused system

  // now turn all the aromatic bonds in this fused system to single
  for (uint32_t bondIdx : makeSingle) {
    bondVec[bondIdx].setBondType(BondEnums::BondType::SINGLE);
  }
}

bool kekulizeWorker(RDMol &mol, const INT_VECT &allAtms,
                    boost::dynamic_bitset<> dBndCands,
                    boost::dynamic_bitset<> dBndAdds, INT_VECT done,
                    const UINT_VECT &atomRanks, unsigned int maxBackTracks) {
  const auto &atomVec = mol.getAtomDataVector();
  auto &bondVec = mol.getBondDataVector();
  INT_DEQUE astack;
  INT_INT_DEQ_MAP options;
  int lastOpt = -1;
  boost::dynamic_bitset<> localBondsAdded(mol.getNumBonds());
  boost::dynamic_bitset<> inAllAtms(mol.getNumAtoms());
  for (int allAtm : allAtms) {
    inAllAtms.set(allAtm);
  }

  auto lessByRank = [&atomRanks](int a, int b) {
    const auto ra = atomRanks.at(static_cast<unsigned int>(a));
    const auto rb = atomRanks.at(static_cast<unsigned int>(b));
    return (ra < rb) || (ra == rb && a < b);
  };

  // Prefer starting traversal at atoms which are the *end* of wedged/dashed
  // bonds. Wedged bonds encode stereo and must remain single bonds; by starting
  // the kekulization walk at wedge-end atoms we assign their double bond to a
  // *different* neighbor first, giving the algorithm more freedom to keep the
  // wedged bond single.
  boost::dynamic_bitset<> wedgeEndAtoms(mol.getNumAtoms());
  for (const BondData &bond : mol.getBondDataVector()) {
    if (bond.getBondDir() == BondEnums::BondDir::BEGINWEDGE ||
        bond.getBondDir() == BondEnums::BondDir::BEGINDASH) {
      const atomindex_t endIdx = bond.getEndAtomIdx();
      if (inAllAtms.test(endIdx)) {
        wedgeEndAtoms.set(endIdx);
      }
    }
  }

  // Pre-sort allAtms: wedge-end atoms first, then by canonical rank.
  // This way the first not-yet-done atom is always the best starting point.
  INT_VECT sortedAtms(allAtms);
  std::sort(sortedAtms.begin(), sortedAtms.end(),
            [&wedgeEndAtoms, &lessByRank](int a, int b) {
              const bool wa = wedgeEndAtoms.test(a);
              const bool wb = wedgeEndAtoms.test(b);
              if (wa != wb) {
                return wa;  // wedge-end atoms come first
              }
              return lessByRank(a, b);
            });

  // ok the algorithm goes something like this
  // - start with an atom that has been marked aromatic before
  // - check if it can have a double bond
  // - add its neighbors to the stack
  // - check if one of its neighbors can also have a double bond
  // - if yes add a double bond.
  // - if multiple neighbors can have double bonds - add them to a
  //   options stack we may have to retrace out path if we chose the
  //   wrong neighbor to add the double bond
  // - if double bond added update the candidates for double bond
  // - move to the next atom on the stack and repeat the process
  // - if an atom that can have multiple a double bond has no
  //   neighbors that can take double bond - we made a mistake
  //   earlier by picking a wrong candidate for double bond
  // - in this case back track to where we made the mistake

  int curr = -1;
  INT_DEQUE btmoves;
  unsigned int numBT = 0;  // number of back tracks so far
  while ((done.size() < sortedAtms.size()) || !astack.empty()) {
    // pick a curr atom to work with
    if (astack.size() > 0) {
      curr = astack.front();
      astack.pop_front();
    } else {
      curr = -1;
      for (int allAtm : sortedAtms) {
        if (std::find(done.begin(), done.end(), allAtm) == done.end()) {
          curr = allAtm;
          break;
        }
      }
    }
    CHECK_INVARIANT(curr >= 0, "starting point not found");
    done.push_back(curr);

    // loop over the neighbors if we can add double bonds or
    // simply push them onto the stack
    INT_DEQUE opts;
    bool cCand = false;
    if (dBndCands[curr]) {
      cCand = true;
    }
    int ncnd;
    // if we are here because of backtracking
    if (options.find(curr) != options.end()) {
      opts = options[curr];
      CHECK_INVARIANT(opts.size() > 0, "");
    } else {
      INT_DEQUE lstack;
      std::vector<int> optsV;
      std::vector<int> wedgedOptsV;
      std::vector<int> nbrs;
      auto [nbrBegin, nbrEnd] = mol.getAtomNeighbors(atomindex_t(curr));
      for (auto nIt = nbrBegin; nIt != nbrEnd; ++nIt) {
        const int nbrIdx = static_cast<int>(*nIt);
        // ignore if the neighbor is not part of the fused system
        if (!inAllAtms.test(nbrIdx)) {
          continue;
        }
        // ignore if the neighbor has already been dealt with before
        if (std::find(done.begin(), done.end(), nbrIdx) != done.end()) {
          continue;
        }
        nbrs.push_back(nbrIdx);
      }

      std::sort(nbrs.begin(), nbrs.end(), lessByRank);

      for (int nbrIdx : nbrs) {
        const uint32_t nbrBondIdx = mol.getBondIndexBetweenAtoms(
            atomindex_t(curr), atomindex_t(nbrIdx));
        const BondData &nbrBond = bondVec[nbrBondIdx];

        // if the neighbor is not on the stack add it
        if (std::find(astack.begin(), astack.end(), nbrIdx) == astack.end()) {
          lstack.push_back(nbrIdx);
        }

        // check if the neighbor is also a candidate for a double bond
        // the refinement that we'll make to the candidate check we've already
        // done is to make sure that the bond is either flagged as aromatic
        // or involves a dummy atom. This was Issue 3525076.
        // This fix is not really 100% of the way there: a situation like
        // that for Issue 3525076 but involving a dummy atom in the cage
        // could lead to the same failure. The full fix would require
        // a fairly detailed analysis of all bonds in the molecule to determine
        // which of them is eligible to be converted.
        if (cCand && dBndCands[nbrIdx] &&
            (nbrBond.getIsAromatic() ||
             atomVec[curr].getAtomicNum() == 0 ||
             atomVec[nbrIdx].getAtomicNum() == 0)) {
          // in order to try and avoid making wedged bonds double, we will add
          // this neighbor at the back of the options after this loop if the
          // bond is wedged. otherwise we append it to the options directly
          if (nbrBond.getBondDir() == BondEnums::BondDir::BEGINWEDGE ||
              nbrBond.getBondDir() == BondEnums::BondDir::BEGINDASH) {
            wedgedOptsV.push_back(nbrIdx);
          } else {
            optsV.push_back(nbrIdx);
          }
        }  // end of curr atoms can have a double bond
      }  // end of looping over neighbors

      // Non-wedged options first, then wedged — both already in rank order
      // because nbrs was pre-sorted by lessByRank above.
      for (int v : optsV) {
        opts.push_back(v);
      }
      for (int v : wedgedOptsV) {
        opts.push_back(v);
      }
      astack.insert(astack.end(), lstack.begin(), lstack.end());
    }
    // now add a double bond from current to one of the neighbors if we can
    if (cCand) {
      if (!opts.empty()) {
        ncnd = opts.front();
        opts.pop_front();
        const uint32_t bondIdx = mol.getBondIndexBetweenAtoms(
            atomindex_t(curr), atomindex_t(ncnd));
        BondData &bnd = bondVec[bondIdx];
        bnd.setBondType(BondEnums::BondType::DOUBLE);
        if (bnd.getBondDir() != BondEnums::BondDir::NONE) {
          bnd.setBondDir(BondEnums::BondDir::NONE);
        }

        // remove current and the neighbor from the dBndCands list
        dBndCands[curr] = 0;
        dBndCands[ncnd] = 0;

        // add them to the list of bonds to which have been made double
        dBndAdds[bondIdx] = 1;
        localBondsAdded[bondIdx] = 1;

        // if this is an atom we previously visted and picked we
        // simply tried a different option now, overwrite the options
        // stored for this atoms
        if (options.find(curr) != options.end()) {
          if (opts.size() == 0) {
            options.erase(curr);
            btmoves.pop_back();
            if (btmoves.size() > 0) {
              lastOpt = btmoves.back();
            } else {
              lastOpt = -1;
            }
          } else {
            options[curr] = opts;
          }
        } else {
          // this is new atoms we are trying and have other
          // neighbors as options to add double bond store this to
          // the options stack, we may have made a mistake in
          // which one we chose and have to return here
          if (opts.size() > 0) {
            lastOpt = curr;
            btmoves.push_back(lastOpt);
            options[curr] = opts;
          }
        }

      }  // end of adding a double bond
      else if (atomVec[curr].getAtomicNum()) {
        // we have a non-dummy atom that should be getting a double
        // bond but none of the neighbors can take one. Most likely
        // because of a wrong choice earlier so back track
        if ((lastOpt >= 0) && (numBT < maxBackTracks)) {
          backTrack(mol, options, lastOpt, done, astack, dBndCands, dBndAdds);
          ++numBT;
        } else {
          // undo any remaining changes we made while here
          // this was github #962
          for (uint32_t bidx = 0, nbnds = uint32_t(bondVec.size()); bidx < nbnds;
               ++bidx) {
            if (localBondsAdded[bidx]) {
              bondVec[bidx].setBondType(BondEnums::BondType::SINGLE);
            }
          }
          return false;
        }
      }  // end of else try to backtrack
    }    // end of curr atom atom being a cand for double bond
  }      // end of while we are not done with all atoms
  return true;
}

class QuestionEnumerator {
 public:
  QuestionEnumerator(INT_VECT questions)
      : d_questions(std::move(questions)), d_pos(1){};
  INT_VECT next() {
    INT_VECT res;
    if (d_pos >= (0x1u << d_questions.size())) {
      return res;
    }
    for (unsigned int i = 0; i < d_questions.size(); ++i) {
      if (d_pos & (0x1u << i)) {
        res.push_back(d_questions[i]);
      }
    }
    ++d_pos;
    return res;
  };

 private:
  INT_VECT d_questions;
  unsigned int d_pos;
};

bool permuteDummiesAndKekulize(RDMol &mol, const INT_VECT &allAtms,
                               boost::dynamic_bitset<> dBndCands,
                               INT_VECT &questions,
                               const UINT_VECT &atomRanks,
                               unsigned int maxBackTracks) {
  boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
  for (int allAtm : allAtms) {
    atomsInPlay[allAtm] = 1;
  }
  bool kekulized = false;
  QuestionEnumerator qEnum(questions);
  while (!kekulized && questions.size()) {
    boost::dynamic_bitset<> dBndAdds(mol.getNumBonds());
    INT_VECT done;
    // reset the state: all aromatic bonds are remarked to single:
    for (BondData &bond : mol.getBondDataVector()) {
      if (bond.getIsAromatic() &&
          bond.getBondType() != BondEnums::BondType::SINGLE &&
          atomsInPlay[bond.getBeginAtomIdx()] &&
          atomsInPlay[bond.getEndAtomIdx()]) {
        bond.setBondType(BondEnums::BondType::SINGLE);
      }
    }
    // pick a new permutation of the questionable atoms:
    const auto &switchOff = qEnum.next();
    if (!switchOff.size()) {
      break;
    }
    auto tCands = dBndCands;
    for (int it : switchOff) {
      tCands[it] = 0;
    }
    // try kekulizing again:
    kekulized =
        kekulizeWorker(mol, allAtms, tCands, dBndAdds, done, atomRanks,
                       maxBackTracks);
  }
  return kekulized;
}

void kekulizeFused(RDMol &mol, const VECT_INT_VECT &arings,
                   const UINT_VECT &atomRanks, unsigned int maxBackTracks) {
  // get all the atoms in the ring system
  INT_VECT allAtms;
  Union(arings, allAtms);
  // get all the atoms that are candidates to receive a double bond
  // also mark atoms in the fused system that are not aromatic to begin with
  // as done. Mark all the bonds that are part of the aromatic system
  // to be single bonds
  INT_VECT done;
  INT_VECT questions;
  const uint32_t nats = mol.getNumAtoms();
  const uint32_t nbnds = mol.getNumBonds();
  boost::dynamic_bitset<> dBndCands(nats);
  boost::dynamic_bitset<> dBndAdds(nbnds);
  markDbondCands(mol, allAtms, dBndCands, questions, done);

  auto kekulized =
      kekulizeWorker(mol, allAtms, dBndCands, dBndAdds, done, atomRanks,
                     maxBackTracks);
  if (!kekulized && questions.size()) {
    // we failed, but there are some dummy atoms we can try permuting.
    kekulized = permuteDummiesAndKekulize(mol, allAtms, dBndCands, questions,
                                          atomRanks, maxBackTracks);
  }
  if (!kekulized) {
    // we exhausted all option (or crossed the allowed
    // number of backTracks) and we still need to backtrack
    // can't kekulize this thing
    std::vector<unsigned int> problemAtoms;
    std::ostringstream errout;
    errout << "Can't kekulize mol.";
    errout << "  Unkekulized atoms:";
    for (unsigned int i = 0; i < nats; ++i) {
      if (dBndCands[i]) {
        errout << " " << i;
        problemAtoms.push_back(i);
      }
    }
    std::string msg = errout.str();
    BOOST_LOG(rdErrorLog) << msg << std::endl;
    throw KekulizeException(msg, problemAtoms);
  }
}
}  // namespace

namespace MolOps {
namespace details {
void KekulizeFragment(RDMol &mol, const boost::dynamic_bitset<> &atomsToUse,
                      boost::dynamic_bitset<> bondsToUse, bool markAtomsBonds,
                      bool canonical, unsigned int maxBackTracks) {
  PRECONDITION(atomsToUse.size() == mol.getNumAtoms(),
               "atomsToUse is wrong size");
  PRECONDITION(bondsToUse.size() == mol.getNumBonds(),
               "bondsToUse is wrong size");
  // if there are no atoms to use we can directly return
  if (atomsToUse.none()) {
    return;
  }

  // there's no point doing kekulization if there are no aromatic bonds
  // without queries:
  bool foundAromatic = false;
  auto &bondVec = mol.getBondDataVector();
  for (uint32_t bondIdx = 0, nbnds = uint32_t(bondVec.size()); bondIdx < nbnds;
       ++bondIdx) {
    if (bondsToUse[bondIdx]) {
      if (QueryOps::hasBondTypeQuery(ConstRDMolBond(&mol, bondIdx))) {
        // we don't kekulize bonds with bond type queries
        bondsToUse[bondIdx] = 0;
      } else if (bondVec[bondIdx].getIsAromatic()) {
        foundAromatic = true;
      }
    }
  }

  // before everything do implicit valence calculation and store them
  // we will repeat after kekulization and compare for the sake of error
  // checking
  const uint32_t numAtoms = mol.getNumAtoms();
  INT_VECT valences(numAtoms);
  boost::dynamic_bitset<> dummyAts(numAtoms);
  const auto &atomVec = mol.getAtomDataVector();

  for (uint32_t atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
    if (!atomsToUse[atomIdx]) {
      continue;
    }
    mol.calcAtomImplicitValence(atomIdx, false);
    const AtomData &atom = atomVec[atomIdx];
    valences[atomIdx] =
        int(atom.getValence(AtomData::ValenceType::EXPLICIT) +
            atom.getValence(AtomData::ValenceType::IMPLICIT));
    if (mol.isAromaticAtom(atomIdx)) {
      foundAromatic = true;
    }
    if (!atom.getAtomicNum()) {
      dummyAts[atomIdx] = 1;
    }
  }
  if (!foundAromatic) {
    return;
  }
  UINT_VECT atomRanks(numAtoms);
  if (canonical) {
    // Canon::rankFragmentAtoms is part of the canonical ranking module
    // (new_canon.h) which is still ROMol-shaped. Bridge through asROMol()
    // for this single call. Documented as the deferred Canon module port
    // in EDGE_CASES.md; the canonical=true path is not exercised by the
    // sanitize hot path (sanitize calls Kekulize(mol, true, false)).
    Canon::rankFragmentAtoms(mol.asROMol(), atomRanks, atomsToUse, bondsToUse);
  } else {
    // When canonical=false (e.g. during sanitization), we skip the
    // expensive ranking step and use atom indices directly.  This is
    // appropriate because sanitization runs *before* stereo perception:
    // canonical ranking would be based on incomplete chemistry and the
    // "deterministic" result would be meaningless.  Callers who need a
    // canonical Kekulé form should call Kekulize() with canonical=true
    // after the molecule is fully sanitized and stereo has been assigned.
    std::iota(atomRanks.begin(), atomRanks.end(), 0u);
  }
  // if any bonds to kekulize then give it a try:
  if (bondsToUse.any()) {
    // mark atoms at the beginning of wedged bonds
    boost::dynamic_bitset<> wedgedAtoms(numAtoms);
    for (uint32_t bondIdx = 0, nbnds = uint32_t(bondVec.size()); bondIdx < nbnds;
         ++bondIdx) {
      const BondData &bond = bondVec[bondIdx];
      if (bondsToUse[bondIdx] &&
          (bond.getBondDir() == BondEnums::BondDir::BEGINWEDGE ||
           bond.getBondDir() == BondEnums::BondDir::BEGINDASH)) {
        wedgedAtoms.set(bond.getBeginAtomIdx());
      }
    }

    // A bit on the state of the molecule at this point
    // - aromatic and non aromatic atoms and bonds may be mixed up

    // - for all aromatic bonds it is assumed that that both the following
    //   are true:
    //       - getIsAromatic returns true
    //       - getBondType return aromatic
    // - all aromatic atoms return true for "getIsAromatic"

    // first find all the simple rings in the molecule that are not
    // completely composed of dummy atoms
    if (!mol.getRingInfo().isInitialized()) {
      MolOps::findSSSR(mol, mol.getRingInfo());
    }
    // Build the rings vector-of-vectors view from the native CSR ring info.
    const RingInfoCache &ringInfoCache = mol.getRingInfo();
    VECT_INT_VECT allrings;
    allrings.reserve(ringInfoCache.numRings());
    for (uint32_t ri = 0, nRings = ringInfoCache.numRings(); ri < nRings;
         ++ri) {
      const uint32_t rb = ringInfoCache.ringBegins[ri];
      const uint32_t re = ringInfoCache.ringBegins[ri + 1];
      INT_VECT ring(re - rb);
      for (uint32_t pos = rb; pos < re; ++pos) {
        ring[pos - rb] = int(ringInfoCache.atomsInRings[pos]);
      }
      allrings.push_back(std::move(ring));
    }

    std::deque<INT_VECT> tmpRings;
    auto containsNonDummy = [&atomsToUse, &dummyAts](const INT_VECT &ring) {
      bool ringOk = false;
      for (auto ai : ring) {
        if (!atomsToUse[ai]) {
          return false;
        }
        if (!dummyAts[ai]) {
          ringOk = true;
        }
      }
      return ringOk;
    };
    // we can't just copy the rings over: we're going to rearrange them so that
    // we try to favor starting the traversal of any ring from an atom that is
    // at the end of a wedged ring bond. This is part of our attempt to avoid
    // assigning double bonds to bonds with wedging
    for (const auto &ring : allrings) {
      if (containsNonDummy(ring)) {
        unsigned int startPos = 0;
        bool hasWedge = false;
        for (auto ri = 0u; ri < ring.size(); ++ri) {
          if (wedgedAtoms[ring[ri]]) {
            startPos = ri;
            hasWedge = true;
            break;
          }
        }
        INT_VECT nring(ring.size());
        for (auto ri = 0u; ri < ring.size(); ++ri) {
          nring[ri] = ring.at((ri + startPos) % ring.size());
        }
        if (!hasWedge) {
          tmpRings.push_back(nring);
        } else {
          tmpRings.push_front(nring);
        }
      }
    }
    VECT_INT_VECT arings;
    arings.reserve(allrings.size());
    arings.insert(arings.end(), tmpRings.begin(), tmpRings.end());
    VECT_INT_VECT allbrings;
    RingUtils::convertToBonds(arings, allbrings, mol);
    VECT_INT_VECT brings;
    brings.reserve(allbrings.size());
    auto copyBondRingsWithinFragment = [&bondsToUse](const INT_VECT &ring) {
      return std::all_of(ring.begin(), ring.end(), [&bondsToUse](const int bi) {
        return bondsToUse[bi];
      });
    };
    VECT_INT_VECT aringsRemaining;
    aringsRemaining.reserve(arings.size());
    for (unsigned i = 0; i < allbrings.size(); ++i) {
      if (copyBondRingsWithinFragment(allbrings[i])) {
        brings.push_back(allbrings[i]);
        aringsRemaining.push_back(arings[i]);
      }
    }
    arings = std::move(aringsRemaining);

    // make a neighbor map for the rings i.e. a ring is a
    // neighbor to another candidate ring if it shares at least
    // one bond
    // useful to figure out fused systems
    INT_INT_VECT_MAP neighMap;
    RingUtils::makeRingNeighborMap(brings, neighMap);

    int curr = 0;
    int cnrs = rdcast<int>(arings.size());
    boost::dynamic_bitset<> fusDone(cnrs);
    while (curr < cnrs) {
      INT_VECT fused;
      RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
      VECT_INT_VECT frings(fused.size());
      std::transform(fused.begin(), fused.end(), frings.begin(),
                     [&arings](const int ri) { return arings[ri]; });
      kekulizeFused(mol, frings, atomRanks, maxBackTracks);
      int rix;
      for (rix = 0; rix < cnrs; ++rix) {
        if (!fusDone[rix]) {
          curr = rix;
          break;
        }
      }
      if (rix == cnrs) {
        break;
      }
    }
  }
  if (markAtomsBonds) {
    // if we want the atoms and bonds to be marked non-aromatic do
    // that here.
    if (!mol.getRingInfo().isInitialized()) {
      MolOps::findSSSR(mol, mol.getRingInfo());
    }
    auto &atomVec = mol.getAtomDataVector();
    for (uint32_t bondIdx = 0, nbnds = uint32_t(bondVec.size());
         bondIdx < nbnds; ++bondIdx) {
      if (bondsToUse[bondIdx]) {
        bondVec[bondIdx].setIsAromatic(false);
      }
    }
    const RingInfoCache &ringInfoCache = mol.getRingInfo();
    for (uint32_t atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
      AtomData &atom = atomVec[atomIdx];
      if (atomsToUse[atomIdx] && atom.getIsAromatic()) {
        // if we're doing the full molecule and there are aromatic atoms not in
        // a ring, throw an exception
        if (atomsToUse.all() && bondsToUse.all() &&
            !ringInfoCache.numAtomRings(atomIdx)) {
          std::ostringstream errout;
          errout << "non-ring atom " << atomIdx << " marked aromatic";
          auto msg = errout.str();
          BOOST_LOG(rdErrorLog) << msg << std::endl;
          throw AtomKekulizeException(msg, atomIdx);
        }
        atom.setIsAromatic(false);
        // make sure "explicit" Hs on things like pyrroles don't hang around
        // this was Github Issue 141
        if ((atom.getAtomicNum() == 7 || atom.getAtomicNum() == 15) &&
            atom.getFormalCharge() == 0 && atom.getNumExplicitHs() == 1) {
          atom.setNoImplicit(false);
          atom.setNumExplicitHs(0);
          mol.updateAtomPropertyCache(atomIdx, false);
        }
      }
    }
  }

  // ok some error checking here force a implicit valence
  // calculation that should do some error checking by itself. In
  // addition compare them to what they were before kekulizing
  for (uint32_t atomIdx = 0; atomIdx < numAtoms; ++atomIdx) {
    if (!atomsToUse[atomIdx]) {
      continue;
    }
    const AtomData &atom = atomVec[atomIdx];
    const int val = int(atom.getValence(AtomData::ValenceType::EXPLICIT) +
                        atom.getValence(AtomData::ValenceType::IMPLICIT));
    if (val != valences[atomIdx]) {
      std::ostringstream errout;
      errout << "Kekulization somehow screwed up valence on " << atomIdx
             << ": " << val << "!=" << valences[atomIdx] << std::endl;
      auto msg = errout.str();
      BOOST_LOG(rdErrorLog) << msg << std::endl;
      throw AtomKekulizeException(msg, atomIdx);
    }
  }
}

void KekulizeFragment(RWMol &mol, const boost::dynamic_bitset<> &atomsToUse,
                      boost::dynamic_bitset<> bondsToUse, bool markAtomsBonds,
                      bool canonical, unsigned int maxBackTracks) {
  KekulizeFragment(mol.asRDMol(), atomsToUse, std::move(bondsToUse),
                   markAtomsBonds, canonical, maxBackTracks);
}
}  // namespace details
void Kekulize(RDMol &mol, bool markAtomsBonds, bool canonical,
              unsigned int maxBackTracks) {
  boost::dynamic_bitset<> atomsToUse(mol.getNumAtoms());
  atomsToUse.set();
  boost::dynamic_bitset<> bondsToUse(mol.getNumBonds());
  bondsToUse.set();
  details::KekulizeFragment(mol, atomsToUse, bondsToUse, markAtomsBonds,
                            canonical, maxBackTracks);
}

void Kekulize(RWMol &mol, bool markAtomsBonds, bool canonical,
              unsigned int maxBackTracks) {
  Kekulize(mol.asRDMol(), markAtomsBonds, canonical, maxBackTracks);
}

bool KekulizeIfPossible(RDMol &mol, bool markAtomsBonds, bool canonical,
                        unsigned int maxBackTracks) {
  auto &bondVec = mol.getBondDataVector();
  boost::dynamic_bitset<> aromaticBonds(bondVec.size());
  for (uint32_t bondIdx = 0, nbnds = uint32_t(bondVec.size()); bondIdx < nbnds;
       ++bondIdx) {
    if (bondVec[bondIdx].getIsAromatic()) {
      aromaticBonds.set(bondIdx);
    }
  }
  boost::dynamic_bitset<> aromaticAtoms(mol.getNumAtoms());
  for (uint32_t atomIdx = 0, numAtoms = mol.getNumAtoms(); atomIdx < numAtoms;
       ++atomIdx) {
    if (mol.isAromaticAtom(atomIdx)) {
      aromaticAtoms.set(atomIdx);
    }
  }
  bool res = true;
  try {
    Kekulize(mol, markAtomsBonds, canonical, maxBackTracks);
  } catch (const MolSanitizeException &) {
    res = false;
    for (uint32_t i = 0, nbnds = uint32_t(bondVec.size()); i < nbnds; ++i) {
      if (aromaticBonds[i]) {
        BondData &bond = bondVec[i];
        bond.setIsAromatic(true);
        bond.setBondType(BondEnums::BondType::AROMATIC);
      }
    }
    auto &atomVec = mol.getAtomDataVector();
    for (uint32_t i = 0, nats = uint32_t(atomVec.size()); i < nats; ++i) {
      if (aromaticAtoms[i]) {
        atomVec[i].setIsAromatic(true);
      }
    }
  }
  return res;
}

bool KekulizeIfPossible(RWMol &mol, bool markAtomsBonds, bool canonical,
                        unsigned int maxBackTracks) {
  return KekulizeIfPossible(mol.asRDMol(), markAtomsBonds, canonical,
                            maxBackTracks);
}
}  // namespace MolOps
}  // namespace RDKit
