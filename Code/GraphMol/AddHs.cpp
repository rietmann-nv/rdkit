//
//  Copyright (C) 2003-2025 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDKitBase.h"
#include <list>
#include "QueryAtom.h"
#include "QueryOps.h"
#include "MonomerInfo.h"
#include "Chirality.h"
#include <Geometry/Transform3D.h>
#include <Geometry/point.h>
#include <boost/algorithm/string/classification.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/range/iterator_range.hpp>

constexpr double sq_dist_zero_tol = 1.e-4;

namespace RDKit {

// Local utility functionality:
namespace {
Atom *getAtomNeighborNot(ROMol *mol, const Atom *atom, const Atom *other) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(atom, "bad atom");
  PRECONDITION(atom->getDegree() > 1, "bad degree");
  PRECONDITION(other, "bad atom");
  Atom *res = nullptr;

  ROMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = mol->getAtomNeighbors(atom);
  while (nbrIdx != endNbrs) {
    if (*nbrIdx != other->getIdx()) {
      res = mol->getAtomWithIdx(*nbrIdx);
      break;
    }
    ++nbrIdx;
  }

  POSTCONDITION(res, "no neighbor found");
  return res;
}

void AssignHsResidueInfo(RWMol &mol) {
  int max_serial = 0;
  unsigned int stopIdx = mol.getNumAtoms();
  for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    auto *info =
        (AtomPDBResidueInfo *)(mol.getAtomWithIdx(aidx)->getMonomerInfo());
    if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE &&
        info->getSerialNumber() > max_serial) {
      max_serial = info->getSerialNumber();
    }
  }

  AtomPDBResidueInfo *current_info = nullptr;
  int current_h_id = 0;
  for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    Atom *newAt = mol.getAtomWithIdx(aidx);
    auto *info = (AtomPDBResidueInfo *)(newAt->getMonomerInfo());
    if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
      ROMol::ADJ_ITER begin, end;
      boost::tie(begin, end) = mol.getAtomNeighbors(newAt);
      while (begin != end) {
        if (mol.getAtomWithIdx(*begin)->getAtomicNum() == 1) {
          // Make all Hs unique - increment id even for existing
          ++current_h_id;
          // skip if hydrogen already has PDB info
          auto *h_info = (AtomPDBResidueInfo *)mol.getAtomWithIdx(*begin)
                             ->getMonomerInfo();
          if (h_info &&
              h_info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
            continue;
          }
          // the hydrogens have unique names on residue basis (H1, H2, ...)
          if (!current_info ||
              current_info->getResidueNumber() != info->getResidueNumber() ||
              current_info->getChainId() != info->getChainId()) {
            current_h_id = 1;
            current_info = info;
          }
          std::string h_label = std::to_string(current_h_id);
          if (h_label.length() > 3) {
            h_label = h_label.substr(h_label.length() - 3, 3);
          }
          while (h_label.length() < 3) {
            h_label = h_label + " ";
          }
          h_label = "H" + h_label;
          // wrap around id to '3H12'
          h_label = h_label.substr(3, 1) + h_label.substr(0, 3);
          AtomPDBResidueInfo *newInfo = new AtomPDBResidueInfo(
              h_label, max_serial, "", info->getResidueName(),
              info->getResidueNumber(), info->getChainId(), "", 1.0, 0.0,
              info->getIsHeteroAtom());
          mol.getAtomWithIdx(*begin)->setMonomerInfo(newInfo);

          ++max_serial;
        }
        ++begin;
      }
    }
  }
}

std::map<unsigned int, std::vector<unsigned int>> getIsoMap(RDMol &mol) {
  std::map<unsigned int, std::vector<unsigned int>> isoMap;
  // Wipe any existing _isotopicHs entries; they will be recomputed from
  // the molecule's current isotope-bearing H neighbors below.
  mol.clearAtomPropIfPresent(common_properties::_isotopicHsToken);
  auto &bondVec = mol.getBondDataVector();
  const auto &atomVec = mol.getAtomDataVector();
  for (uint32_t bondIdx = 0, numBonds = uint32_t(bondVec.size());
       bondIdx < numBonds; ++bondIdx) {
    const BondData &bond = bondVec[bondIdx];
    const atomindex_t baIdx = bond.getBeginAtomIdx();
    const atomindex_t eaIdx = bond.getEndAtomIdx();
    const AtomData &ba = atomVec[baIdx];
    const AtomData &ea = atomVec[eaIdx];
    int ha = -1;
    unsigned int iso = 0;
    if (ba.getAtomicNum() == 1 && ba.getIsotope() && ea.getAtomicNum() != 1) {
      ha = eaIdx;
      iso = ba.getIsotope();
    } else if (ea.getAtomicNum() == 1 && ea.getIsotope() &&
               ba.getAtomicNum() != 1) {
      ha = baIdx;
      iso = ea.getIsotope();
    }
    if (ha == -1) {
      continue;
    }
    isoMap[ha].push_back(iso);
  }
  return isoMap;
}

int atomPerturbationOrder(const RDMol &mol, atomindex_t atomIdx,
                          const INT_LIST &probe) {
  // Replicates Atom::getPerturbationOrder against the flat adjacency.
  const size_t numBonds = mol.getAtomDegree(atomIdx);
  PRECONDITION(numBonds == probe.size(), "size mismatch");
  auto [beginBonds, endBonds] = mol.getAtomBonds(atomIdx);
  std::vector<int> copy(probe.begin(), probe.end());

  // The native bond list is sorted by construction; leave the safety check
  // in place to match the legacy fallback.
  bool isSorted = true;
  for (size_t i = 0; numBonds > 0 && i < numBonds - 1; ++i) {
    isSorted = isSorted && (beginBonds[i] < beginBonds[i + 1]);
  }
  if (isSorted) {
    int nSwaps = 0;
    for (size_t desti = 0, n = copy.size(); desti < n; ++desti) {
      size_t besti = desti;
      for (size_t srci = desti + 1; srci < n; ++srci) {
        if (copy[srci] < copy[besti]) {
          besti = srci;
        }
      }
      if (besti != desti) {
        ++nSwaps;
        std::swap(copy[besti], copy[desti]);
      }
    }
    return nSwaps;
  }
  std::vector<int> bondsCopy(beginBonds, endBonds);
  return static_cast<int>(countSwapsToInterconvert(bondsCopy, copy));
}

bool may_need_extra_H(const RDMol &mol, atomindex_t atomIdx) {
  unsigned single_bonds = 0;
  unsigned aromatic_bonds = 0;
  const auto &bondVec = mol.getBondDataVector();
  auto [bondBegin, bondEnd] = mol.getAtomBonds(atomIdx);
  for (auto it = bondBegin; it != bondEnd; ++it) {
    const BondData &bond = bondVec[*it];
    if (bond.getBondType() == BondEnums::BondType::SINGLE) {
      ++single_bonds;
    } else if (bond.getBondType() == BondEnums::BondType::AROMATIC) {
      ++aromatic_bonds;
    } else {
      return false;
    }
  }
  const AtomData &atom = mol.getAtomDataVector()[atomIdx];
  const unsigned int totalValence =
      atom.getValence(AtomData::ValenceType::EXPLICIT) +
      atom.getValence(AtomData::ValenceType::IMPLICIT);
  return single_bonds == 1 && aromatic_bonds == 2 && totalValence == 3;
}

}  // end of unnamed namespace

namespace MolOps {

namespace {
RDGeom::Point3D pickBisector(const RDGeom::Point3D &nbr1Vect,
                             const RDGeom::Point3D &nbr2Vect,
                             const RDGeom::Point3D &nbr3Vect) {
  auto dirVect = nbr2Vect + nbr3Vect;
  if (dirVect.lengthSq() < sq_dist_zero_tol) {
    // nbr2Vect and nbr3Vect are anti-parallel (was #3854)
    dirVect = nbr2Vect;
    std::swap(dirVect.x, dirVect.y);
    dirVect.x *= -1;
  }
  if (dirVect.dotProduct(nbr1Vect) < 0) {
    dirVect *= -1;
  }

  return dirVect;
}
}  // namespace

void setTerminalAtomCoords(ROMol &mol, unsigned int idx,
                           unsigned int otherIdx) {
  // we will loop over all the coordinates
  PRECONDITION(otherIdx != idx, "degenerate atoms");
  Atom *atom = mol.getAtomWithIdx(idx);
  PRECONDITION(mol.getAtomDegree(atom) == 1, "bad atom degree");
  const Bond *bond = mol.getBondBetweenAtoms(otherIdx, idx);
  PRECONDITION(bond, "no bond between atoms");

  const Atom *otherAtom = mol.getAtomWithIdx(otherIdx);
  double bondLength =
      PeriodicTable::getTable()->getRb0(1) +
      PeriodicTable::getTable()->getRb0(otherAtom->getAtomicNum());

  RDGeom::Point3D dirVect(0, 0, 0);

  RDGeom::Point3D perpVect, rotnAxis, nbrPerp;
  RDGeom::Point3D nbr1Vect, nbr2Vect, nbr3Vect;
  RDGeom::Transform3D tform;
  RDGeom::Point3D otherPos, atomPos;

  const Atom *nbr1 = nullptr, *nbr2 = nullptr, *nbr3 = nullptr;
  const Bond *nbrBond;
  ROMol::ADJ_ITER nbrIdx, endNbrs;

  switch (otherAtom->getDegree()) {
    case 1:
      // --------------------------------------------------------------------------
      //   No other atoms present:
      // --------------------------------------------------------------------------
      // loop over the conformations and set the coordinates
      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           cfi++) {
        if ((*cfi)->is3D()) {
          dirVect.z = 1;
        } else {
          dirVect.x = 1;
        }
        otherPos = (*cfi)->getAtomPos(otherIdx);
        atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
        (*cfi)->setAtomPos(idx, atomPos);
      }
      break;

    case 2:
      // --------------------------------------------------------------------------
      //  One other neighbor:
      // --------------------------------------------------------------------------
      nbr1 = getAtomNeighborNot(&mol, otherAtom, atom);
      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           ++cfi) {
        otherPos = (*cfi)->getAtomPos(otherIdx);
        RDGeom::Point3D nbr1Pos = (*cfi)->getAtomPos(nbr1->getIdx());
        // get a normalized vector pointing away from the neighbor:
        nbr1Vect = nbr1Pos - otherPos;
        if (nbr1Vect.lengthSq() < sq_dist_zero_tol) {
          // no difference, which likely indicates that we have redundant atoms.
          // just put it on top of the heavy atom. This was #678
          (*cfi)->setAtomPos(idx, otherPos);
          continue;
        }
        nbr1Vect.normalize();
        nbr1Vect *= -1;

        // ok, nbr1Vect points away from the other atom, figure out where
        // this H goes:
        switch (otherAtom->getHybridization()) {
          case Atom::SP3:
            // get a perpendicular to nbr1Vect:
            if ((*cfi)->is3D()) {
              perpVect = nbr1Vect.getPerpendicular();
            } else {
              perpVect.z = 1.0;
            }
            // and move off it:
            tform.SetRotation((180 - 109.471) * M_PI / 180., perpVect);
            dirVect = tform * nbr1Vect;
            atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(idx, atomPos);
            break;
          case Atom::SP2:
            // default 3D position is to just take an arbitrary perpendicular
            // for 2D we take the normal to the xy plane
            if ((*cfi)->is3D()) {
              perpVect = nbr1Vect.getPerpendicular();
            } else {
              perpVect.z = 1.0;
            }
            if (nbr1->getDegree() > 1) {
              // can we use the neighboring atom to establish a perpendicular?
              nbrBond = mol.getBondBetweenAtoms(otherIdx, nbr1->getIdx());
              if (nbrBond->getIsAromatic() ||
                  nbrBond->getBondType() == Bond::DOUBLE ||
                  nbrBond->getIsConjugated()) {
                nbr2 = getAtomNeighborNot(&mol, nbr1, otherAtom);
                nbr2Vect =
                    nbr1Pos.directionVector((*cfi)->getAtomPos(nbr2->getIdx()));
                auto crossProd = nbr2Vect.crossProduct(nbr1Vect);

                // if nbr1 and nbr2 are aligned, the perpendicular will be null,
                // and we'll just keep the default calculated above. Otherwise
                // we use the cross product
                if (crossProd.lengthSq() >= sq_dist_zero_tol) {
                  perpVect = crossProd;
                }
              }
            }
            perpVect.normalize();
            // rotate the nbr1Vect 60 degrees about perpVect and we're done:
            tform.SetRotation(60. * M_PI / 180., perpVect);
            dirVect = tform * nbr1Vect;
            atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(idx, atomPos);
            break;
          case Atom::SP:
            // just lay the H along the vector:
            dirVect = nbr1Vect;
            atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(idx, atomPos);
            break;
          default:
            // FIX: handle other hybridizations
            // for now, just lay the H along the vector:
            dirVect = nbr1Vect;
            atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
            (*cfi)->setAtomPos(idx, atomPos);
        }
      }
      break;
    case 3:
      // --------------------------------------------------------------------------
      // Two other neighbors:
      // --------------------------------------------------------------------------
      boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(otherAtom);
      while (nbrIdx != endNbrs) {
        if (*nbrIdx != idx) {
          if (!nbr1) {
            nbr1 = mol.getAtomWithIdx(*nbrIdx);
          } else {
            nbr2 = mol.getAtomWithIdx(*nbrIdx);
          }
        }
        ++nbrIdx;
      }
      TEST_ASSERT(nbr1);
      TEST_ASSERT(nbr2);
      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           ++cfi) {
        // start along the average of the two vectors:
        otherPos = (*cfi)->getAtomPos(otherIdx);
        nbr1Vect = otherPos - (*cfi)->getAtomPos(nbr1->getIdx());
        nbr2Vect = otherPos - (*cfi)->getAtomPos(nbr2->getIdx());
        if (nbr1Vect.lengthSq() < sq_dist_zero_tol ||
            nbr2Vect.lengthSq() < sq_dist_zero_tol) {
          // no difference, which likely indicates that we have redundant atoms.
          // just put it on top of the heavy atom. This was #678
          (*cfi)->setAtomPos(idx, otherPos);
          continue;
        }
        nbr1Vect.normalize();
        nbr2Vect.normalize();
        dirVect = nbr1Vect + nbr2Vect;

        if (dirVect.lengthSq() < sq_dist_zero_tol) {
          // nbr1Vect and nbr2Vect are non-null, but they may
          // still cancel each other out
          continue;
        }
        dirVect.normalize();
        if ((*cfi)->is3D()) {
          switch (otherAtom->getHybridization()) {
            case Atom::SP3:
              // get the perpendicular to the neighbors:
              nbrPerp = nbr1Vect.crossProduct(nbr2Vect);
              // and the perpendicular to that:
              rotnAxis = nbrPerp.crossProduct(dirVect);
              // and then rotate about that:
              rotnAxis.normalize();
              tform.SetRotation((109.471 / 2) * M_PI / 180., rotnAxis);
              dirVect = tform * dirVect;
              atomPos =
                  otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
              (*cfi)->setAtomPos(idx, atomPos);
              break;
            case Atom::SP2:
              // don't need to do anything here, the H atom goes right on the
              // direction vector
              atomPos =
                  otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
              (*cfi)->setAtomPos(idx, atomPos);
              break;
            default:
              // FIX: handle other hybridizations
              // for now, just lay the H along the neighbor vector;
              atomPos =
                  otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
              (*cfi)->setAtomPos(idx, atomPos);
              break;
          }
        } else {
          // don't need to do anything here, the H atom goes right on the
          // direction vector
          atomPos = otherPos + dirVect;
          (*cfi)->setAtomPos(idx, atomPos);
        }
      }
      break;
    case 4:
      // --------------------------------------------------------------------------
      // Three other neighbors:
      // --------------------------------------------------------------------------
      boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(otherAtom);

      // We're using chiral tag for checking chirality, so we just take the
      // initial order
      while (nbrIdx != endNbrs) {
        if (*nbrIdx != idx) {
          if (!nbr1) {
            nbr1 = mol.getAtomWithIdx(*nbrIdx);
          } else if (!nbr2) {
            nbr2 = mol.getAtomWithIdx(*nbrIdx);
          } else {
            nbr3 = mol.getAtomWithIdx(*nbrIdx);
          }
        }
        ++nbrIdx;
      }

      TEST_ASSERT(nbr1);
      TEST_ASSERT(nbr2);
      TEST_ASSERT(nbr3);

      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           ++cfi) {
        otherPos = (*cfi)->getAtomPos(otherIdx);
        nbr1Vect = otherPos - (*cfi)->getAtomPos(nbr1->getIdx());
        nbr2Vect = otherPos - (*cfi)->getAtomPos(nbr2->getIdx());
        nbr3Vect = otherPos - (*cfi)->getAtomPos(nbr3->getIdx());
        if (nbr1Vect.lengthSq() < sq_dist_zero_tol ||
            nbr2Vect.lengthSq() < sq_dist_zero_tol ||
            nbr3Vect.lengthSq() < sq_dist_zero_tol) {
          // no difference, which likely indicates that we have redundant atoms.
          // just put it on top of the heavy atom. This was #678
          (*cfi)->setAtomPos(idx, otherPos);
          continue;
        }
        nbr1Vect.normalize();
        nbr2Vect.normalize();
        nbr3Vect.normalize();

        // if three neighboring atoms are more or less planar, this
        // is going to be in a quasi-random (but almost definitely bad)
        // direction...
        // correct for this (issue 2951221):
        if ((*cfi)->is3D()) {
          if (fabs(nbr3Vect.dotProduct(nbr1Vect.crossProduct(nbr2Vect))) <
              0.1) {
            // compute the normal:
            dirVect = nbr1Vect.crossProduct(nbr2Vect);

            // Each of the nbr vectors is non-null, but there might be pairs
            // that cancel each other out. Try to find a direction from atoms
            // that do not overlap.
            if (dirVect.lengthSq() < sq_dist_zero_tol) {
              // This definition of dirVect reverses the parity around otherIdx
              // the change of sign restores it
              dirVect = nbr1Vect.crossProduct(nbr3Vect) * -1;
            }
            if (dirVect.lengthSq() < sq_dist_zero_tol) {
              dirVect = nbr2Vect.crossProduct(nbr3Vect);
            }
            // We couldn't find a good direction
            if (dirVect.lengthSq() < sq_dist_zero_tol) {
              continue;
            }

            std::string cipCode;
            if (otherAtom->getPropIfPresent(common_properties::_CIPCode,
                                            cipCode)) {
              // the heavy atom is a chiral center, make sure
              // that we went go the right direction to preserve
              // its chirality. We use the chiral volume for this:
              RDGeom::Point3D v1 = dirVect - nbr3Vect;
              RDGeom::Point3D v2 = nbr1Vect - nbr3Vect;
              RDGeom::Point3D v3 = nbr2Vect - nbr3Vect;
              double vol = v1.dotProduct(v2.crossProduct(v3));

              if ((otherAtom->getChiralTag() ==
                       Atom::ChiralType::CHI_TETRAHEDRAL_CCW &&
                   vol < 0) ||
                  (otherAtom->getChiralTag() ==
                       Atom::ChiralType::CHI_TETRAHEDRAL_CW &&
                   vol > 0)) {
                dirVect *= -1;
              }
            }
          } else {
            dirVect = nbr1Vect + nbr2Vect + nbr3Vect;
          }
        } else {
          // we're in flatland

          // github #3879 and #908: find the two neighbors with the largest
          // outer angle between them and then place the H to bisect that angle
          // This is recommendation ST-1.1.4 from the 2006 IUPAC "Graphical
          // representation of stereochemical configuration" guideline
          auto angle12 = nbr1Vect.angleTo(nbr2Vect);
          auto angle13 = nbr1Vect.angleTo(nbr3Vect);
          auto angle23 = nbr2Vect.angleTo(nbr3Vect);
          auto accum1 = angle12 + angle13;
          auto accum2 = angle12 + angle23;
          auto accum3 = angle13 + angle23;
          if (accum1 <= accum2 && accum1 <= accum3) {
            dirVect = pickBisector(nbr1Vect, nbr2Vect, nbr3Vect);
          } else if (accum2 <= accum1 && accum2 <= accum3) {
            dirVect = pickBisector(nbr2Vect, nbr1Vect, nbr3Vect);
          } else {
            dirVect = pickBisector(nbr3Vect, nbr1Vect, nbr2Vect);
          }
        }

        dirVect.normalize();
        atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
        (*cfi)->setAtomPos(idx, atomPos);
      }
      break;
    default:
      // --------------------------------------------------------------------------
      // FIX: figure out what to do here
      // --------------------------------------------------------------------------
      atomPos = otherPos + dirVect * bondLength;
      for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
           ++cfi) {
        (*cfi)->setAtomPos(idx, atomPos);
      }
      break;
  }
}

namespace {
bool isQueryAtom(const RDMol &mol, atomindex_t atomIdx) {
  if (mol.hasAtomQuery(atomIdx)) {
    return true;
  }
  auto [bondBegin, bondEnd] = mol.getAtomBonds(atomIdx);
  for (auto it = bondBegin; it != bondEnd; ++it) {
    if (mol.hasBondQuery(*it)) {
      return true;
    }
  }
  return false;
}
}  // namespace
void addHs(RDMol &mol, const AddHsParameters &params,
           const UINT_VECT *onlyOnAtoms) {
  // when we hit each atom, clear its computed properties
  // NOTE: it is essential that we not clear the ring info in the
  // molecule's computed properties.  We don't want to have to
  // regenerate that.  This caused Issue210 and Issue212:
  mol.clearComputedProps(false);

  const uint32_t numHeavyAtoms = mol.getNumAtoms();
  boost::dynamic_bitset<> onAtoms(numHeavyAtoms);
  if (onlyOnAtoms) {
    for (auto atIdx : *onlyOnAtoms) {
      onAtoms.set(atIdx);
    }
  } else {
    onAtoms.set();
  }
  std::vector<unsigned int> numExplicitHs(numHeavyAtoms, 0);
  std::vector<unsigned int> numImplicitHs(numHeavyAtoms, 0);

  // Ensure all atoms have implicit valence calculated (e.g., after unpickling)
  // Skip query atoms (or atoms connected to query bonds) as they are handled
  // separately
  {
    auto &atomVec = mol.getAtomDataVector();
    for (uint32_t aidx = 0; aidx < numHeavyAtoms; ++aidx) {
      if (!atomVec[aidx].getNoImplicit() && !isQueryAtom(mol, aidx)) {
        mol.updateAtomPropertyCache(aidx, false);
      }
    }
  }

  // Cache H counts before adding bonds (which will clear property cache per
  // #8934). Total Hs precomputed so we can pre-size the conformer storage.
  unsigned int numAddHyds = 0;
  {
    const auto &atomVec = mol.getAtomDataVector();
    for (uint32_t aidx = 0; aidx < numHeavyAtoms; ++aidx) {
      const AtomData &atom = atomVec[aidx];
      numExplicitHs[aidx] = atom.getNumExplicitHs();
      numImplicitHs[aidx] = atom.getNumImplicitHs();
      if (onAtoms[aidx]) {
        if (params.skipQueries && isQueryAtom(mol, aidx)) {
          onAtoms.set(aidx, 0);
          continue;
        }
        numAddHyds += atom.getNumExplicitHs();
        if (!params.explicitOnly) {
          numAddHyds += atom.getNumImplicitHs();
        }
      }
    }
  }
  const unsigned int nSize = numHeavyAtoms + numAddHyds;
  if (mol.getNumConformers() > 0) {
    mol.allocateConformers(mol.getNumConformers(), nSize);
  }

  // Add Hs one heavy atom at a time. Note: pulling atomVec/bondVec references
  // before the loop is unsafe — addAtom/addBond can reallocate, invalidating
  // references — so we re-fetch as needed inside the loop.
  for (uint32_t aidx = 0; aidx < numHeavyAtoms; ++aidx) {
    if (!onAtoms[aidx]) {
      continue;
    }

    std::vector<unsigned int> isoHs;
    if (mol.getAtomPropIfPresent(common_properties::_isotopicHsToken, aidx,
                                 isoHs)) {
      mol.clearSingleAtomProp(common_properties::_isotopicHsToken, aidx);
    }
    auto isoH = isoHs.begin();

    // Snapshot the explicit-H count up front, then zero it on the heavy atom
    // before adding bonds (the addBond loop below will refresh derived
    // properties; we need the explicit count cleared so updateAtomPropertyCache
    // recomputes the new implicit valence correctly).
    const unsigned int onumexpl = numExplicitHs[aidx];
    for (unsigned int i = 0; i < onumexpl; i++) {
      AtomData &hAtomData = mol.addAtom();
      hAtomData.setAtomicNum(1);
      const atomindex_t newIdx = atomindex_t(mol.getNumAtoms() - 1);
      mol.addBond(aidx, newIdx, BondEnums::BondType::SINGLE);
      mol.updateAtomPropertyCache(newIdx, true);
      if (params.addCoords) {
        ROMol &romol = mol.asROMol();
        setTerminalAtomCoords(romol, newIdx, aidx);
      }
      if (isoH != isoHs.end()) {
        mol.getAtomDataVector()[newIdx].setIsotope(*isoH);
        ++isoH;
      }
    }
    mol.getAtomDataVector()[aidx].setNumExplicitHs(0);

    if (!params.explicitOnly) {
      for (unsigned int i = 0; i < numImplicitHs[aidx]; i++) {
        AtomData &hAtomData = mol.addAtom();
        hAtomData.setAtomicNum(1);
        const atomindex_t newIdx = atomindex_t(mol.getNumAtoms() - 1);
        mol.addBond(aidx, newIdx, BondEnums::BondType::SINGLE);
        // set the isImplicit label so that we can strip these back
        // off later if need be.
        mol.setSingleAtomProp(common_properties::isImplicitToken, newIdx, 1);
        mol.updateAtomPropertyCache(newIdx, true);
        if (params.addCoords) {
          ROMol &romol = mol.asROMol();
          setTerminalAtomCoords(romol, newIdx, aidx);
        }
        if (isoH != isoHs.end()) {
          mol.getAtomDataVector()[newIdx].setIsotope(*isoH);
          ++isoH;
        }
      }
    }
    // refresh the heavy atom's derived properties (valence count, etc.) —
    // non-strict per github #2782
    mol.updateAtomPropertyCache(aidx, false);
    if (isoH != isoHs.end()) {
      BOOST_LOG(rdWarningLog)
          << "extra H isotope information found on atom " << aidx << std::endl;
    }
  }
  // take care of AtomPDBResidueInfo for Hs if root atom has it
  if (params.addResidueInfo) {
    RWMol &rwmol = static_cast<RWMol &>(mol.asROMol());
    AssignHsResidueInfo(rwmol);
  }
}

void addHs(RWMol &mol, const AddHsParameters &params,
           const UINT_VECT *onlyOnAtoms) {
  addHs(mol.asRDMol(), params, onlyOnAtoms);
}

namespace {
// returns whether or not an adjustment was made, in case we want that info
bool adjustStereoAtomsIfRequired(RDMol &mol, atomindex_t atomIdx,
                                 atomindex_t heavyAtomIdx) {
  // nothing we can do if the degree is only 2 (and we should have covered
  // that earlier anyway)
  if (mol.getAtomDegree(heavyAtomIdx) == 2) {
    return false;
  }
  if (mol.getBondIndexBetweenAtoms(atomIdx, heavyAtomIdx) ==
      std::numeric_limits<uint32_t>::max()) {
    return false;
  }
  auto [hbBegin, hbEnd] = mol.getAtomBonds(heavyAtomIdx);
  for (auto bIt = hbBegin; bIt != hbEnd; ++bIt) {
    BondData &bnd = mol.getBond(*bIt);
    if (bnd.getBondType() != BondEnums::BondType::DOUBLE ||
        bnd.getStereo() <= BondEnums::BondStereo::STEREOANY) {
      continue;
    }
    if (!mol.hasBondStereoAtoms(*bIt)) {
      continue;
    }
    const atomindex_t *stereoAtoms = mol.getBondStereoAtoms(*bIt);
    const atomindex_t storedA = stereoAtoms[0];
    const atomindex_t storedB = stereoAtoms[1];
    int slot = -1;
    if (storedA == atomIdx) {
      slot = 0;
    } else if (storedB == atomIdx) {
      slot = 1;
    } else {
      continue;
    }
    // find another neighbor of the heavy atom (not the double-bond partner,
    // not the H being removed) and substitute it into the stereo-atom slot.
    const atomindex_t dblNbrIdx = bnd.getOtherAtomIdx(heavyAtomIdx);
    auto [nbrBegin, nbrEnd] = mol.getAtomNeighbors(heavyAtomIdx);
    for (auto nIt = nbrBegin; nIt != nbrEnd; ++nIt) {
      const uint32_t nbrIdx = *nIt;
      if (nbrIdx == dblNbrIdx || nbrIdx == atomIdx) {
        continue;
      }
      const atomindex_t newA =
          (slot == 0) ? atomindex_t(nbrIdx) : storedA;
      const atomindex_t newB =
          (slot == 1) ? atomindex_t(nbrIdx) : storedB;
      mol.setBondStereoAtoms(*bIt, newA, newB);
      bool madeAdjustment = true;
      switch (bnd.getStereo()) {
        case BondEnums::BondStereo::STEREOCIS:
          bnd.setStereo(BondEnums::BondStereo::STEREOTRANS);
          break;
        case BondEnums::BondStereo::STEREOTRANS:
          bnd.setStereo(BondEnums::BondStereo::STEREOCIS);
          break;
        default:
          // I think we shouldn't need to do anything with E and Z...
          madeAdjustment = false;
          break;
      }
      return madeAdjustment;
    }
  }
  return false;
}

void molRemoveH(RDMol &mol, atomindex_t idx, bool updateExplicitCount) {
  AtomData &atom = mol.getAtom(idx);
  PRECONDITION(atom.getAtomicNum() == 1, "idx corresponds to a non-Hydrogen");

  // Snapshot the bond indices for this H atom; we'll iterate over them
  // while mutating bond / atom state but the adjacency itself isn't
  // changing until removeAtom at the end.
  auto [hBondBegin, hBondEnd] = mol.getAtomBonds(idx);
  std::vector<uint32_t> hBondIndices(hBondBegin, hBondEnd);

  for (const uint32_t bondIdx : hBondIndices) {
    BondData &bond = mol.getBond(bondIdx);
    const atomindex_t heavyAtomIdx = bond.getOtherAtomIdx(idx);
    AtomData &heavyAtom = mol.getAtom(heavyAtomIdx);
    const int heavyAtomNum = heavyAtom.getAtomicNum();

    // we'll update the neighbor's explicit H count if we were told to
    // *or* if the neighbor is chiral, in which case the H is needed
    // in order to complete the coordination
    // *or* if the neighbor has the noImplicit flag set:
    if (updateExplicitCount || heavyAtom.getNoImplicit() ||
        heavyAtom.getChiralTag() != AtomEnums::ChiralType::CHI_UNSPECIFIED) {
      heavyAtom.setNumExplicitHs(heavyAtom.getNumExplicitHs() + 1);
    } else {
      // this is a special case related to Issue 228 and the
      // "disappearing Hydrogen" problem discussed in MolOps::adjustHs
      //
      // If we remove a hydrogen from an aromatic N or P, or if
      // the heavy atom it is connected to is not in its default
      // valence state, we need to be *sure* to increment the
      // explicit count, even if the H itself isn't marked as explicit
      const INT_VECT &defaultVs =
          PeriodicTable::getTable()->getValenceList(heavyAtomNum);
      const bool aromatic = mol.isAromaticAtom(heavyAtomIdx);
      const unsigned int heavyTotalValence =
          heavyAtom.getValence(AtomData::ValenceType::EXPLICIT) +
          heavyAtom.getValence(AtomData::ValenceType::IMPLICIT);
      if (((heavyAtomNum == 7 || heavyAtomNum == 15 ||
            may_need_extra_H(mol, heavyAtomIdx)) &&
           aromatic) ||
          (std::find(defaultVs.begin() + 1, defaultVs.end(),
                     int(heavyTotalValence)) != defaultVs.end())) {
        heavyAtom.setNumExplicitHs(heavyAtom.getNumExplicitHs() + 1);
      }
    }

    // One other consequence of removing the H from the graph is that we
    // may change the ordering of the bonds about a chiral center, which
    // may change the chiral label at that atom. Handle it explicitly.
    if (heavyAtom.getChiralTag() != AtomEnums::ChiralType::CHI_UNSPECIFIED) {
      INT_LIST neighborBondIndices;
      auto [habBegin, habEnd] = mol.getAtomBonds(heavyAtomIdx);
      for (auto it = habBegin; it != habEnd; ++it) {
        if (*it != bondIdx) {
          neighborBondIndices.push_back(int(*it));
        }
      }
      neighborBondIndices.push_back(int(bondIdx));
      const int nSwaps =
          atomPerturbationOrder(mol, heavyAtomIdx, neighborBondIndices);
      if (nSwaps % 2) {
        // invert chirality
        switch (heavyAtom.getChiralTag()) {
          case AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW:
            heavyAtom.setChiralTag(
                AtomEnums::ChiralType::CHI_TETRAHEDRAL_CCW);
            break;
          case AtomEnums::ChiralType::CHI_TETRAHEDRAL_CCW:
            heavyAtom.setChiralTag(
                AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW);
            break;
          default:
            break;
        }
      }
    }

    // If we are removing a H atom that defines bond stereo (e.g. imines),
    // also remove the bond stereo information, as it is no longer valid.
    if (mol.getAtomDegree(heavyAtomIdx) == 2) {
      auto [habBegin, habEnd] = mol.getAtomBonds(heavyAtomIdx);
      for (auto it = habBegin; it != habEnd; ++it) {
        if (*it != bondIdx) {
          BondData &nbnd = mol.getBond(*it);
          if (nbnd.getStereo() > BondEnums::BondStereo::STEREOANY) {
            nbnd.setStereo(BondEnums::BondStereo::STEREONONE);
            mol.clearBondStereoAtoms(*it);
          }
          break;
        }
      }
    }

    // wavy bond → mark the beginning atom with _UnknownStereo
    if (bond.getBondDir() == BondEnums::BondDir::UNKNOWN &&
        bond.getBeginAtomIdx() == heavyAtomIdx) {
      mol.setSingleAtomProp(common_properties::_UnknownStereoToken,
                            heavyAtomIdx, 1);
    } else if (bond.getBondDir() == BondEnums::BondDir::ENDDOWNRIGHT ||
               bond.getBondDir() == BondEnums::BondDir::ENDUPRIGHT) {
      // if the direction is set on this bond and the atom it's connected to
      // has no other single bonds with directions set, then we need to set
      // direction on one of the other neighbors in order to avoid double
      // bond stereochemistry possibly being lost. This was github #754
      bool foundADir = false;
      uint32_t oBondIdx = std::numeric_limits<uint32_t>::max();
      auto [habBegin, habEnd] = mol.getAtomBonds(heavyAtomIdx);
      for (auto it = habBegin; it != habEnd; ++it) {
        if (*it == bondIdx) {
          continue;
        }
        BondData &nbnd = mol.getBond(*it);
        if (nbnd.getBondType() == BondEnums::BondType::SINGLE) {
          if (nbnd.getBondDir() == BondEnums::BondDir::NONE) {
            oBondIdx = *it;
          } else {
            foundADir = true;
          }
        }
      }
      if (!foundADir && oBondIdx != std::numeric_limits<uint32_t>::max()) {
        BondData &oBond = mol.getBond(oBondIdx);
        const bool flipIt = (oBond.getBeginAtomIdx() == heavyAtomIdx) &&
                            (bond.getBeginAtomIdx() == heavyAtomIdx);
        if (flipIt) {
          oBond.setBondDir(
              bond.getBondDir() == BondEnums::BondDir::ENDDOWNRIGHT
                  ? BondEnums::BondDir::ENDUPRIGHT
                  : BondEnums::BondDir::ENDDOWNRIGHT);
        } else {
          oBond.setBondDir(bond.getBondDir());
        }
      }
      adjustStereoAtomsIfRequired(mol, idx, heavyAtomIdx);
    } else {
      adjustStereoAtomsIfRequired(mol, idx, heavyAtomIdx);
    }

    // remove the bond from any SGroups that might include it.
    for (auto &sg : getSubstanceGroups(mol)) {
      sg.removeBondWithIdx(bondIdx);
    }
  }

  // Finally, remove the atom from any SGroups that might include it, so the
  // SGroups don't get removed in removeAtom(). Since we allow removing
  // SGroup SAP lvidx H atoms, check for those and update them.
  for (auto &sg : getSubstanceGroups(mol)) {
    sg.removeAtomWithIdx(idx);
    sg.removeParentAtomWithIdx(idx);
    for (auto &sap : sg.getAttachPoints()) {
      if (sap.lvIdx == static_cast<int>(idx)) {
        sap.lvIdx = -1;
      }
    }
  }
  // computed properties will be cleared after all hydrogens are removed
  const bool clearProps = false;
  mol.removeAtom(idx, clearProps);
}

bool shouldRemoveH(const RDMol &mol, atomindex_t atomIdx,
                   const RemoveHsParameters &ps) {
  const AtomData &atom = mol.getAtom(atomIdx);
  if (atom.getAtomicNum() != 1) {
    return false;
  }
  if (!ps.removeWithQuery && mol.getAtomQuery(atomIdx) != nullptr) {
    return false;
  }
  const uint32_t degree = mol.getAtomDegree(atomIdx);
  if (!ps.removeDegreeZero && degree == 0) {
    if (ps.showWarnings) {
      BOOST_LOG(rdWarningLog)
          << "WARNING: not removing hydrogen atom without neighbors"
          << std::endl;
    }
    return false;
  }
  if (!ps.removeHigherDegrees && degree > 1) {
    return false;
  }
  if (!ps.removeIsotopes && !ps.removeAndTrackIsotopes && atom.getIsotope()) {
    return false;
  }
  if (!ps.removeNonimplicit &&
      !mol.hasAtomProp(common_properties::isImplicitToken, atomIdx)) {
    return false;
  }
  if (!ps.removeMapped && mol.getAtomMapNum(atomIdx)) {
    return false;
  }

  if (ps.removeInSGroups) {
    // If removing H in SGroups, do not remove H atoms in special
    // roles in the SGroup
    for (const auto &sg : getSubstanceGroups(mol)) {
      // The H atom is one of the "caps" of the SGroup. Technically,
      // it's not part of the group, but it defines its boundaries.
      for (const auto &bond_idx : sg.getBonds()) {
        if (sg.getBondType(bond_idx) == SubstanceGroup::BondType::XBOND) {
          const BondData &bond = mol.getBond(bond_idx);
          if (bond.getBeginAtomIdx() == atomIdx ||
              bond.getEndAtomIdx() == atomIdx) {
            return false;
          }
        }
      }

      for (const auto &sap : sg.getAttachPoints()) {
        // The H atom is an attach point. This would be weird, but is possible.
        // (if it is a 'leaving atom' we don't care, though)
        if (sap.aIdx == atomIdx) {
          return false;
        }
      }

      for (const auto &cs : sg.getCStates()) {
        // The bond to the H atom defines a CState
        const BondData &bond = mol.getBond(cs.bondIdx);
        if (bond.getBeginAtomIdx() == atomIdx ||
            bond.getEndAtomIdx() == atomIdx) {
          return false;
        }
      }
    }
  } else {
    for (const auto &sg : getSubstanceGroups(mol)) {
      if (sg.includesAtom(atomIdx)) {
        return false;
      }
    }
  }
  if (!ps.removeHydrides && atom.getFormalCharge() == -1) {
    return false;
  }
  bool removeIt = true;
  if (degree &&
      (!ps.removeDummyNeighbors || !ps.removeDefiningBondStereo ||
       !ps.removeOnlyHNeighbors || !ps.removeNontetrahedralNeighbors ||
       !ps.removeWithWedgedBond)) {
    bool onlyHNeighbors = true;
    auto [nbrBegin, nbrEnd] = mol.getAtomNeighbors(atomIdx);
    for (auto it = nbrBegin; it != nbrEnd; ++it) {
      const uint32_t nbrIdx = *it;
      const AtomData &nbr = mol.getAtom(nbrIdx);
      // is it a dummy?
      if (!ps.removeDummyNeighbors && nbr.getAtomicNum() < 1) {
        if (ps.showWarnings) {
          BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
                                     "with dummy atom neighbors"
                                  << std::endl;
        }
        return false;
      }
      // does it have non-tetrahedral stereo:
      if (!ps.removeNontetrahedralNeighbors &&
          Chirality::hasNonTetrahedralStereo(nbr)) {
        if (ps.showWarnings) {
          BOOST_LOG(rdWarningLog)
              << "WARNING: not removing hydrogen atom "
                 "with neighbor that has non-tetrahedral stereochemistry"
              << std::endl;
        }
        return false;
      }
      if (!ps.removeOnlyHNeighbors && nbr.getAtomicNum() != 1) {
        onlyHNeighbors = false;
      }
      const uint32_t bondToNbrIdx =
          mol.getBondIndexBetweenAtoms(atomIdx, nbrIdx);
      const BondData &bondToNbr = mol.getBond(bondToNbrIdx);
      if (!ps.removeWithWedgedBond) {
        if (bondToNbr.getBondDir() == BondEnums::BondDir::BEGINDASH ||
            bondToNbr.getBondDir() == BondEnums::BondDir::BEGINWEDGE) {
          if (ps.showWarnings) {
            BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
                                       "with wedged bond"
                                    << std::endl;
          }
          return false;
        }
      }
      // Check to see if the neighbor has a double bond and we're the only
      // neighbor at this end.  This was part of github #1810
      if (!ps.removeDefiningBondStereo && mol.getAtomDegree(nbrIdx) == 2) {
        auto [nbrBondBegin, nbrBondEnd] = mol.getAtomBonds(nbrIdx);
        for (auto bIt = nbrBondBegin; bIt != nbrBondEnd; ++bIt) {
          const BondData &nbrBnd = mol.getBond(*bIt);
          if (nbrBnd.getBondType() == BondEnums::BondType::DOUBLE &&
              (nbrBnd.getStereo() > BondEnums::BondStereo::STEREOANY ||
               bondToNbr.getBondDir() > BondEnums::BondDir::NONE)) {
            return false;
          }
        }
      }
    }
    if (removeIt && (!ps.removeOnlyHNeighbors && onlyHNeighbors)) {
      return false;
    }
  }
  return removeIt;
}

// Do not remove H atoms that are part of SGroups that only contain H atoms.
void filter_sgroup_emptying_hydrogens(const RDMol &mol,
                                      boost::dynamic_bitset<> &atomsToRemove) {
  for (const auto &sg : getSubstanceGroups(mol)) {
    const auto &atoms = sg.getAtoms();
    const auto &patoms = sg.getParentAtoms();

    // If the SGroup already didn't have atoms, we don't care about it
    if (atoms.empty() && patoms.empty()) {
      continue;
    }

    auto would_remove_atom = [&atomsToRemove](const auto idx) {
      return atomsToRemove[idx];
    };

    auto no_atoms = atoms.empty() ||
                    std::all_of(atoms.begin(), atoms.end(), would_remove_atom);
    if (no_atoms) {
      auto no_patoms =
          patoms.empty() ||
          std::all_of(patoms.begin(), patoms.end(), would_remove_atom);
      if (no_patoms) {
        for (auto atom : atoms) {
          atomsToRemove.set(atom, false);
        }
        for (auto patom : patoms) {
          atomsToRemove.set(patom, false);
        }
      }
    }
  }
}

}  // end of anonymous namespace

void removeHs(RDMol &mol, const RemoveHsParameters &ps,
              MolOps::SanitizeTemp &temp, bool sanitize) {
  // The SanitizeTemp workspace is currently unused by the body itself; it's
  // threaded through the API so that an eventual sanitizeMol(RDMol&,
  // SanitizeTemp&) port can reuse the buffers without changing this
  // function's signature.
  (void)temp;

  if (ps.removeAndTrackIsotopes) {
    // if there are any non-isotopic Hs remove them first to make sure
    // chirality is preserved.
    bool needRemoveHs = false;
    for (const auto &atom : mol.getAtomDataVector()) {
      if (atom.getAtomicNum() == 1 && atom.getIsotope() == 0) {
        needRemoveHs = true;
        break;
      }
    }
    if (needRemoveHs) {
      RemoveHsParameters psCopy(ps);
      psCopy.removeAndTrackIsotopes = false;
      psCopy.removeIsotopes = false;
      removeHs(mol, psCopy, temp, false);
    }
  }
  for (uint32_t atomIdx = 0, numAtoms = mol.getNumAtoms(); atomIdx < numAtoms;
       ++atomIdx) {
    mol.updateAtomPropertyCache(atomIdx, false);
  }
  if (ps.removeAndTrackIsotopes) {
    for (const auto &pair : getIsoMap(mol)) {
      mol.setSingleAtomProp(common_properties::_isotopicHsToken,
                            atomindex_t(pair.first), pair.second);
    }
  }
  boost::dynamic_bitset<> atomsToRemove{mol.getNumAtoms(), 0};

  for (uint32_t atomIdx = 0, numAtoms = mol.getNumAtoms(); atomIdx < numAtoms;
       ++atomIdx) {
    if (shouldRemoveH(mol, atomIdx, ps)) {
      atomsToRemove.set(atomIdx);
    }
  }

  // Once we know which H atoms would be removed, filter out those that
  // would cause any SGroups to become empty.
  if (ps.removeInSGroups) {
    filter_sgroup_emptying_hydrogens(mol, atomsToRemove);
  }

  // Remove the marked atoms one at a time. NOTE: stereochemistry handling
  // makes batch editing unsafe here.
  for (int idx = int(mol.getNumAtoms()) - 1; idx >= 0; --idx) {
    if (atomsToRemove[idx]) {
      molRemoveH(mol, atomindex_t(idx), ps.updateExplicitCount);
    }
  }
  mol.clearComputedProps(true);

  // If we removed non-implicit Hs, atom indices may have shifted in ways
  // that invalidate derived properties (ring membership, etc); re-sanitize.
  if (!atomsToRemove.empty() && ps.removeNonimplicit && sanitize) {
    sanitizeMol(mol);
  }

  // If we removed Hs and any chiral atom ended up with more than one
  // explicit H, drop those redundant counts.
  if (!atomsToRemove.empty()) {
    auto &atomVec = mol.getAtomDataVector();
    for (uint32_t atomIdx = 0, numAtoms = uint32_t(atomVec.size());
         atomIdx < numAtoms; ++atomIdx) {
      AtomData &atom = atomVec[atomIdx];
      if (!atom.getNoImplicit() &&
          atom.getChiralTag() != AtomEnums::ChiralType::CHI_UNSPECIFIED) {
        if (atom.getNumExplicitHs() > 1) {
          atom.setNumExplicitHs(0);
          mol.updateAtomPropertyCache(atomIdx, false);
        }
      }
    }
  }
}

void removeHs(RWMol &mol, const RemoveHsParameters &ps, bool sanitize) {
  MolOps::SanitizeTemp temp;
  removeHs(mol.asRDMol(), ps, temp, sanitize);
}
ROMol *removeHs(const ROMol &mol, const RemoveHsParameters &ps, bool sanitize) {
  auto *res = new RWMol(mol);
  try {
    removeHs(*res, ps, sanitize);
  } catch (const MolSanitizeException &) {
    delete res;
    throw;
  }
  return static_cast<ROMol *>(res);
}
void removeHs(RDMol &mol, MolOps::SanitizeTemp &temp, bool implicitOnly,
              bool updateExplicitCount, bool sanitize) {
  RemoveHsParameters ps;
  ps.removeNonimplicit = !implicitOnly;
  ps.updateExplicitCount = updateExplicitCount;
  removeHs(mol, ps, temp, sanitize);
}
void removeHs(RWMol &mol, bool implicitOnly, bool updateExplicitCount,
              bool sanitize) {
  RemoveHsParameters ps;
  ps.removeNonimplicit = !implicitOnly;
  ps.updateExplicitCount = updateExplicitCount;
  removeHs(mol, ps, sanitize);
};
ROMol *removeHs(const ROMol &mol, bool implicitOnly, bool updateExplicitCount,
                bool sanitize) {
  auto *res = new RWMol(mol);
  RemoveHsParameters ps;
  ps.removeNonimplicit = !implicitOnly;
  ps.updateExplicitCount = updateExplicitCount;
  try {
    removeHs(*res, ps, sanitize);
  } catch (const MolSanitizeException &) {
    delete res;
    throw;
  }
  return static_cast<ROMol *>(res);
}

void removeAllHs(RWMol &mol, bool sanitize) {
  RemoveHsParameters ps;
  ps.removeDegreeZero = true;
  ps.removeHigherDegrees = true;
  ps.removeOnlyHNeighbors = true;
  ps.removeIsotopes = true;
  ps.removeDummyNeighbors = true;
  ps.removeDefiningBondStereo = true;
  ps.removeWithWedgedBond = true;
  ps.removeWithQuery = true;
  ps.removeNonimplicit = true;
  ps.removeInSGroups = true;
  ps.showWarnings = false;
  ps.removeHydrides = true;
  ps.removeNontetrahedralNeighbors = true;
  removeHs(mol, ps, sanitize);
};
ROMol *removeAllHs(const ROMol &mol, bool sanitize) {
  auto *res = new RWMol(mol);
  try {
    removeAllHs(*res, sanitize);
  } catch (const MolSanitizeException &) {
    delete res;
    throw;
  }
  return static_cast<ROMol *>(res);
}

namespace {
enum class HydrogenType {
  NotAHydrogen,
  UnMergableQueryHydrogen,
  QueryHydrogen
};

template <class Q>
std::pair<bool, bool> queryHasHs(Q queryAtom, bool inor = false) {
  for (auto childit = queryAtom->beginChildren();
       childit != queryAtom->endChildren(); ++childit) {
    QueryAtom::QUERYATOM_QUERY::CHILD_TYPE query = *childit;
    if (query->getDescription() == "AtomOr") {
      return queryHasHs(query, true);
    } else if (query->getDescription() == "AtomAtomicNum") {
      if (static_cast<ATOM_EQUALS_QUERY *>(query.get())->getVal() == 1 &&
          !query->getNegation()) {
        return std::make_pair(true, inor);
      }
    } else if (query->getDescription() == "AtomType") {
      auto val = static_cast<ATOM_EQUALS_QUERY *>(query.get())->getVal();
      // 1001 == aromtic hydrogen (not a thing, really)
      // 1 == aliphatic hydrogen
      if ((val == 1001 || val == 1) && !query->getNegation()) {
        return std::make_pair(true, inor);
      }
    }
  }
  return std::make_pair(false, inor);
  ;
}

HydrogenType isQueryH(const Atom *atom) {
  PRECONDITION(atom, "bogus atom");
  if (atom->getAtomicNum() == 1) {
    // the simple case: the atom is flagged as being an H and
    // has no query
    if (!atom->hasQuery() ||
        (!atom->getQuery()->getNegation() &&
         atom->getQuery()->getDescription() == "AtomAtomicNum")) {
      return HydrogenType::QueryHydrogen;
    }
  }

  if (!(atom->getDegree() <= 1)) {
    // bonded and unbonded H atoms will continue rest will be returned
    return HydrogenType::NotAHydrogen;
  }

  if (atom->hasQuery() && atom->getQuery()->getNegation()) {
    // we will not merge negated queries
    return HydrogenType::NotAHydrogen;
  }

  if (atom->hasQuery()) {
    std::pair<bool, bool> res = std::make_pair(false, false);
    if (atom->getQuery()->getDescription() == "AtomOr") {
      res = queryHasHs(atom->getQuery(), true);
    } else if (atom->getQuery()->getDescription() == "AtomAnd") {
      res = queryHasHs(atom->getQuery(), false);
    }
    if (res.first) {     // hasH
      if (res.second) {  // inOr
        BOOST_LOG(rdWarningLog)
            << "WARNING: merging explicit H queries involved "
               "in ORs is not supported. This query will not "
               "be merged"
            << std::endl;
        return HydrogenType::UnMergableQueryHydrogen;
      } else {
        return HydrogenType::QueryHydrogen;
      }
    }
  }
  return HydrogenType::NotAHydrogen;
}
}  // namespace

//
//  This routine removes explicit hydrogens (and bonds to them) from
//  the molecular graph and adds them as queries to the heavy atoms
//  to which they are bound.  If the heavy atoms (or atom queries)
//  already have hydrogen-count queries, they will be updated.
//
//  NOTE:
//   - Hydrogens which aren't connected to a heavy atom will not be
//     removed.  This prevents molecules like "[H][H]" from having
//     all atoms removed.
//
//   - By default all hydrogens are removed, however if
//     merge_unmapped_only is true, any hydrogen participating
//     in an atom map will be retained
void mergeQueryHs(RWMol &mol, bool mergeUnmappedOnly, bool mergeIsotopes) {
  std::vector<unsigned int> atomsToRemove;

  boost::dynamic_bitset<> hatoms(mol.getNumAtoms());
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    hatoms[i] = isQueryH(mol.getAtomWithIdx(i)) == HydrogenType::QueryHydrogen;
  }
  unsigned int currIdx = 0, stopIdx = mol.getNumAtoms();
  while (currIdx < stopIdx) {
    Atom *atom = mol.getAtomWithIdx(currIdx);
    if (!hatoms[currIdx]) {
      unsigned int numHsToRemove = 0;
      ROMol::ADJ_ITER begin, end;
      boost::tie(begin, end) = mol.getAtomNeighbors(atom);

      while (begin != end) {
        if (hatoms[*begin]) {
          Atom &bgn = *mol.getAtomWithIdx(*begin);
          bool checkUnmapped =
              !mergeUnmappedOnly ||
              !bgn.hasProp(common_properties::molAtomMapNumber);
          bool checkIsotope = mergeIsotopes || bgn.getIsotope() == 0;
          if (checkUnmapped && checkIsotope) {
            atomsToRemove.push_back(rdcast<unsigned int>(*begin));
            ++numHsToRemove;
          }
        }
        ++begin;
      }
      if (numHsToRemove) {
        //
        //  We have H neighbors:
        //   Add the appropriate queries to compensate for their removal.
        //
        //  Examples:
        //    C[H] -> [C;!H0]
        //    C([H])[H] -> [C;!H0;!H1]
        //
        //  It would be more efficient to do this using range queries like:
        //    C([H])[H] -> [C;H{2-}]
        //  but that would produce non-standard SMARTS without the user
        //  having started with a non-standard SMARTS.
        //
        if (!atom->hasQuery()) {
          // it wasn't a query atom, we need to replace it so that we can add
          // a query:
          ATOM_EQUALS_QUERY *tmp = makeAtomNumQuery(atom->getAtomicNum());
          auto *newAt = new QueryAtom;
          newAt->setQuery(tmp);
          newAt->updateProps(*atom);
          mol.replaceAtom(atom->getIdx(), newAt);
          delete newAt;
          atom = mol.getAtomWithIdx(currIdx);
        }
        for (unsigned int i = 0; i < numHsToRemove; ++i) {
          ATOM_EQUALS_QUERY *tmp = makeAtomHCountQuery(i);
          tmp->setNegation(true);
          atom->expandQuery(tmp);
        }
      }  // end of numHsToRemove test

      // recurse if needed (was github isusue 544)
      if (atom->hasQuery()) {
        if (atom->getQuery()->getDescription() == "RecursiveStructure") {
          auto *rsq = dynamic_cast<RecursiveStructureQuery *>(atom->getQuery());
          CHECK_INVARIANT(rsq, "could not convert recursive structure query");
          RWMol *rqm = new RWMol(*rsq->getQueryMol());
          mergeQueryHs(*rqm, mergeUnmappedOnly, mergeIsotopes);
          rsq->setQueryMol(rqm);
        }

        // FIX: shouldn't be repeating this code here
        std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack(
            atom->getQuery()->beginChildren(), atom->getQuery()->endChildren());
        while (childStack.size()) {
          QueryAtom::QUERYATOM_QUERY::CHILD_TYPE qry = childStack.front();
          childStack.pop_front();
          if (qry->getDescription() == "RecursiveStructure") {
            auto *rsq = dynamic_cast<RecursiveStructureQuery *>(qry.get());
            CHECK_INVARIANT(rsq, "could not convert recursive structure query");
            RWMol *rqm = new RWMol(*rsq->getQueryMol());
            mergeQueryHs(*rqm, mergeUnmappedOnly, mergeIsotopes);
            rsq->setQueryMol(rqm);
          } else if (qry->beginChildren() != qry->endChildren()) {
            childStack.insert(childStack.end(), qry->beginChildren(),
                              qry->endChildren());
          }
        }
      }  // end of recursion loop
    }
    ++currIdx;
  }
  mol.beginBatchEdit();
  for (auto aidx : atomsToRemove) {
    mol.removeAtom(aidx);
  }
  mol.commitBatchEdit();
};
ROMol *mergeQueryHs(const ROMol &mol, bool mergeUnmappedOnly,
                    bool mergeIsotopes) {
  auto *res = new RWMol(mol);
  mergeQueryHs(*res, mergeUnmappedOnly, mergeIsotopes);
  return static_cast<ROMol *>(res);
};

bool needsHs(const ROMol &mol) {
  for (const auto atom : mol.atoms()) {
    bool includeNeighbors = false;
    if (atom->getTotalNumHs(includeNeighbors)) {
      return true;
    }
  }
  return false;
}

std::pair<bool, bool> hasQueryHs(const ROMol &mol) {
  bool queryHs = false;
  // We don't care about announcing ORs or other items during isQueryH
  RDLog::LogStateSetter blocker;

  for (const auto atom : mol.atoms()) {
    switch (isQueryH(atom)) {
      case HydrogenType::UnMergableQueryHydrogen:
        return std::make_pair(true, true);
      case HydrogenType::QueryHydrogen:
        queryHs = true;
        break;
      default:  // HydrogenType::NotAHydrogen:
        break;
    }
    if (atom->hasQuery()) {
      if (atom->getQuery()->getDescription() == "RecursiveStructure") {
        auto *rsq = dynamic_cast<RecursiveStructureQuery *>(atom->getQuery());
        CHECK_INVARIANT(rsq, "could not convert recursive structure query");
        auto res = hasQueryHs(*rsq->getQueryMol());
        if (res.second) {  // unmergableH implies queryH
          return res;
        }
        queryHs |= res.first;
      }

      // FIX: shouldn't be repeating this code here -- yet again!
      std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack(
          atom->getQuery()->beginChildren(), atom->getQuery()->endChildren());
      while (!childStack.empty()) {
        QueryAtom::QUERYATOM_QUERY::CHILD_TYPE qry = childStack.front();
        childStack.pop_front();
        if (qry->getDescription() == "RecursiveStructure") {
          auto *rsq = dynamic_cast<RecursiveStructureQuery *>(qry.get());
          CHECK_INVARIANT(rsq, "could not convert recursive structure query");
          auto res = hasQueryHs(*rsq->getQueryMol());
          if (res.second) {
            return res;
          }
          queryHs |= res.first;
        } else {
          childStack.insert(childStack.end(), qry->beginChildren(),
                            qry->endChildren());
        }
      }
    }
  }  // end of recursion loop

  return std::make_pair(queryHs, false);
}

}  // namespace MolOps
}  // namespace RDKit
