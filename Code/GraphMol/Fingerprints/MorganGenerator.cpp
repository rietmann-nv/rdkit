//
//  Copyright (C) 2018-2025 Boran Adas and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <RDGeneral/hash/hash.hpp>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <tuple>

#include <GraphMol/Fingerprints/FingerprintUtil.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace MorganFingerprint {

using namespace MorganFingerprints;

MorganAtomInvGenerator::MorganAtomInvGenerator(const bool includeRingMembership)
    : df_includeRingMembership(includeRingMembership) {}

std::vector<std::uint32_t> *MorganAtomInvGenerator::getAtomInvariants(
    const RDMol &mol) const {
  const std::uint32_t nAtoms = mol.getNumAtoms();
  std::unique_ptr<std::vector<std::uint32_t>> atomInvariants(
      new std::vector<std::uint32_t>(nAtoms));
  getConnectivityInvariants(mol, *atomInvariants, df_includeRingMembership);
  return atomInvariants.release();
}

std::string MorganAtomInvGenerator::infoString() const {
  return "MorganInvariantGenerator includeRingMembership=" +
         std::to_string(df_includeRingMembership);
}
void MorganAtomInvGenerator::toJSON(boost::property_tree::ptree &pt) const {
  pt.put("type", "MorganAtomInvGenerator");
  pt.put("includeRingMembership", df_includeRingMembership);
  AtomInvariantsGenerator::toJSON(pt);
}
void MorganAtomInvGenerator::fromJSON(const boost::property_tree::ptree &pt) {
  df_includeRingMembership =
      pt.get<bool>("includeRingMembership", df_includeRingMembership);
  AtomInvariantsGenerator::fromJSON(pt);
}

MorganAtomInvGenerator *MorganAtomInvGenerator::clone() const {
  return new MorganAtomInvGenerator(df_includeRingMembership);
}

MorganFeatureAtomInvGenerator::MorganFeatureAtomInvGenerator(
    const std::vector<const ROMol *> *patterns) {
  if (patterns) {
    dp_patterns = new std::vector<const ROMol *>;
    dp_patterns->reserve(patterns->size());
    for (auto pattern : *patterns) {
      dp_patterns->push_back(new ROMol(*pattern));
    }
  }
}

void MorganFeatureAtomInvGenerator::cleanUpPatterns() {
  if (dp_patterns) {
    for (auto mol : *dp_patterns) {
      delete mol;
    }
    delete dp_patterns;
    dp_patterns = nullptr;
  }
}

MorganFeatureAtomInvGenerator::~MorganFeatureAtomInvGenerator() {
  cleanUpPatterns();
}

std::string MorganFeatureAtomInvGenerator::infoString() const {
  return "MorganFeatureInvariantGenerator";
}
void MorganFeatureAtomInvGenerator::toJSON(
    boost::property_tree::ptree &pt) const {
  pt.put("type", "MorganFeatureAtomInvGenerator");
  if (dp_patterns) {
    boost::property_tree::ptree patternsNode;
    for (const auto &pattern : *dp_patterns) {
      boost::property_tree::ptree patternNode;
      std::string smarts = MolToSmarts(*pattern);
      patternNode.put("", smarts);
      patternsNode.push_back(std::make_pair("", patternNode));
    }
    pt.add_child("patternSMARTS", patternsNode);
  }
  AtomInvariantsGenerator::toJSON(pt);
}
void MorganFeatureAtomInvGenerator::fromJSON(
    const boost::property_tree::ptree &pt) {
  if (pt.get_child_optional("patternSMARTS")) {
    const auto &patternsNode = pt.get_child("patternSMARTS");
    cleanUpPatterns();
    dp_patterns = new std::vector<const ROMol *>();
    for (const auto &patternNode : patternsNode) {
      std::string smarts = patternNode.second.get_value<std::string>();
      ROMol *patternMol = SmartsToMol(smarts);
      if (patternMol) {
        dp_patterns->push_back(patternMol);
      }
    }
  }
  AtomInvariantsGenerator::fromJSON(pt);
}

MorganFeatureAtomInvGenerator *MorganFeatureAtomInvGenerator::clone() const {
  return new MorganFeatureAtomInvGenerator(dp_patterns);
}

std::vector<std::uint32_t> *MorganFeatureAtomInvGenerator::getAtomInvariants(
    const RDMol &mol) const {
  const std::uint32_t nAtoms = mol.getNumAtoms();
  std::vector<std::uint32_t> *result = new std::vector<std::uint32_t>(nAtoms);
  getFeatureInvariants(mol, *result, dp_patterns);
  return result;
}

MorganBondInvGenerator::MorganBondInvGenerator(const bool useBondTypes,
                                               const bool useChirality)
    : df_useBondTypes(useBondTypes), df_useChirality(useChirality) {}

std::vector<std::uint32_t> *MorganBondInvGenerator::getBondInvariants(
    const RDMol &rdmol) const {
  using BondEnums::BondType;
  using BondEnums::BondStereo;
  const auto &bondData = rdmol.getBondDataVector();
  const std::uint32_t numBonds = static_cast<std::uint32_t>(bondData.size());
  auto *result = new std::vector<std::uint32_t>(numBonds);
  for (std::uint32_t i = 0; i < numBonds; ++i) {
    const BondData &bond = bondData[i];
    std::int32_t bondInvariant = 1;
    if (df_useBondTypes) {
      if (!df_useChirality || bond.getBondType() != BondType::DOUBLE ||
          bond.getStereo() == BondStereo::STEREONONE) {
        bondInvariant = static_cast<std::int32_t>(bond.getBondType());
      } else {
        auto bondStereo = static_cast<std::int32_t>(bond.getStereo());
        if (!Chirality::getUseLegacyStereoPerception()) {
          // CIPLabeler is currently ROMol-only; bridge once if the labels are
          // not yet computed. See EDGE_CASES.md (CIPLabeler) for the deferred
          // native port that would remove this `mol.asROMol()` call.
          if (!rdmol.asROMol().hasProp(common_properties::_CIPComputed)) {
            CIPLabeler::assignCIPLabels(const_cast<ROMol &>(rdmol.asROMol()));
          }
          std::string cipCode;
          if (rdmol.getBondPropIfPresent(common_properties::_CIPCodeToken, i,
                                         cipCode)) {
            if (cipCode == "E") {
              bondStereo = static_cast<std::int32_t>(BondStereo::STEREOE);
            } else if (cipCode == "Z") {
              bondStereo = static_cast<std::int32_t>(BondStereo::STEREOZ);
            }
          }
        }
        const std::int32_t stereoOffset = 100;
        const std::int32_t bondTypeOffset = 10;
        bondInvariant =
            stereoOffset +
            bondTypeOffset * static_cast<std::int32_t>(bond.getBondType()) +
            bondStereo;
      }
    }
    (*result)[i] = static_cast<std::uint32_t>(bondInvariant);
  }
  return result;
}

std::string MorganBondInvGenerator::infoString() const {
  return "MorganInvariantGenerator useBondTypes=" +
         std::to_string(df_useBondTypes) +
         " useChirality=" + std::to_string(df_useChirality);
}
void MorganBondInvGenerator::toJSON(boost::property_tree::ptree &pt) const {
  pt.put("type", "MorganBondInvGenerator");
  pt.put("useBondTypes", df_useBondTypes);
  pt.put("useChirality", df_useChirality);
  BondInvariantsGenerator::toJSON(pt);
}
void MorganBondInvGenerator::fromJSON(const boost::property_tree::ptree &pt) {
  df_useBondTypes = pt.get<bool>("useBondTypes", df_useBondTypes);
  df_useChirality = pt.get<bool>("useChirality", df_useChirality);
  BondInvariantsGenerator::fromJSON(pt);
}

MorganBondInvGenerator *MorganBondInvGenerator::clone() const {
  return new MorganBondInvGenerator(df_useBondTypes, df_useChirality);
}

template <typename OutputType>
OutputType MorganEnvGenerator<OutputType>::getResultSize() const {
  return std::numeric_limits<OutputType>::max();
}

std::string MorganArguments::infoString() const {
  return "MorganArguments onlyNonzeroInvariants=" +
         std::to_string(df_onlyNonzeroInvariants) +
         " radius=" + std::to_string(d_radius);
}
void MorganArguments::toJSON(boost::property_tree::ptree &pt) const {
  pt.put("type", "MorganArguments");
  pt.put("onlyNonzeroInvariants", df_onlyNonzeroInvariants);
  pt.put("radius", d_radius);
  FingerprintArguments::toJSON(pt);
}
void MorganArguments::fromJSON(const boost::property_tree::ptree &pt) {
  d_radius = pt.get<std::uint32_t>("radius", d_radius);
  df_onlyNonzeroInvariants =
      pt.get<bool>("onlyNonzeroInvariants", df_onlyNonzeroInvariants);
  FingerprintArguments::fromJSON(pt);
}

template <typename OutputType>
void MorganAtomEnv<OutputType>::updateAdditionalOutput(
    AdditionalOutput *additionalOutput, size_t bitId) const {
  PRECONDITION(additionalOutput, "bad output pointer");
  PRECONDITION(d_mol, "bad mol pointer");
  if (additionalOutput->bitInfoMap) {
    (*additionalOutput->bitInfoMap)[bitId].emplace_back(d_atomId, d_layer);
  }
  if (additionalOutput->atomCounts) {
    (*additionalOutput->atomCounts)[d_atomId]++;
  }
  if (additionalOutput->atomToBits) {
    (*additionalOutput->atomToBits)[d_atomId].push_back(bitId);
  }
  if (additionalOutput->atomsPerBit) {
    std::vector<int> atomsInvolved;
    atomsInvolved.push_back(d_atomId);
    if (d_layer > 0) {
      const auto dm = MolOps::getDistanceMat(*d_mol);
      const unsigned int nAtoms = d_mol->getNumAtoms();
      for (unsigned int i = 0; i < nAtoms; ++i) {
        if (static_cast<unsigned int>(dm[d_atomId * nAtoms + i] + .1) <=
                d_layer &&
            i != d_atomId) {
          atomsInvolved.push_back(i);
        }
      }
    }
    (*additionalOutput->atomsPerBit)[bitId].push_back(std::move(atomsInvolved));
  }
}

template <typename OutputType>
OutputType MorganAtomEnv<OutputType>::getBitId(
    FingerprintArguments *,              // arguments
    const std::vector<std::uint32_t> *,  // atomInvariants
    const std::vector<std::uint32_t> *,  // bondInvariants
    AdditionalOutput *,                  // additional Output
    const bool,                          // hashResults
    const std::uint64_t                  // fpSize
) const {
  return d_code;
}  // namespace MorganFingerprint

template <typename OutputType>
std::vector<AtomEnvironment<OutputType> *>
MorganEnvGenerator<OutputType>::getEnvironments(
    const RDMol &mol, FingerprintArguments *arguments,
    const std::vector<std::uint32_t> *fromAtoms,
    const std::vector<std::uint32_t> *,  // ignoreAtoms
    const int,                           // confId
    const AdditionalOutput *,            // additionalOutput
    const std::vector<std::uint32_t> *atomInvariants,
    const std::vector<std::uint32_t> *bondInvariants,
    const bool  // hashResults
) const {
  PRECONDITION(atomInvariants && (atomInvariants->size() >= mol.getNumAtoms()),
               "bad atom invariants size");
  PRECONDITION(bondInvariants && (bondInvariants->size() >= mol.getNumBonds()),
               "bad bond invariants size");
  auto *morganArguments = dynamic_cast<MorganArguments *>(arguments);
  PRECONDITION(morganArguments, "bad arguments type");

  const unsigned int nAtoms = mol.getNumAtoms();
  const unsigned int nBonds = mol.getNumBonds();
  const unsigned int maxNumResults = (morganArguments->d_radius + 1) * nAtoms;

  std::vector<AtomEnvironment<OutputType> *> result;
  result.reserve(maxNumResults);

  // CIPLabeler is currently ROMol-only; bridge once if the labels are not yet
  // computed. See EDGE_CASES.md (CIPLabeler) for the deferred native port.
  if (morganArguments->df_includeChirality &&
      !Chirality::getUseLegacyStereoPerception() &&
      !mol.asROMol().hasProp(common_properties::_CIPComputed)) {
    CIPLabeler::assignCIPLabels(const_cast<ROMol &>(mol.asROMol()));
  }

  std::vector<OutputType> currentInvariants(atomInvariants->size());
  std::copy(atomInvariants->begin(), atomInvariants->end(),
            currentInvariants.begin());
  std::vector<OutputType> nextLayerInvariants(nAtoms);

  std::vector<std::pair<std::int32_t, std::uint32_t>> neighborhoodInvariants;
  neighborhoodInvariants.reserve(8);

  boost::dynamic_bitset<> includeAtoms(nAtoms);
  if (fromAtoms) {
    for (auto idx : *fromAtoms) {
      includeAtoms.set(idx, 1);
    }
  } else {
    includeAtoms.set();
  }

  boost::dynamic_bitset<> chiralAtoms(nAtoms);

  std::unordered_set<boost::dynamic_bitset<>> neighborhoods;
  neighborhoods.reserve(maxNumResults);
  std::vector<boost::dynamic_bitset<>> atomNeighborhoods(
      nAtoms, boost::dynamic_bitset<>(nBonds));
  std::vector<boost::dynamic_bitset<>> roundAtomNeighborhoods =
      atomNeighborhoods;
  boost::dynamic_bitset<> deadAtoms(nAtoms);

  // if df_onlyNonzeroInvariants is set order the atoms to make sure atoms
  // with zero invariants are processed last so that in case of duplicate
  // environments atoms with non-zero invariants are used
  std::vector<unsigned int> atomOrder(nAtoms);
  if (morganArguments->df_onlyNonzeroInvariants) {
    std::vector<std::pair<int32_t, uint32_t>> ordering;
    for (unsigned int i = 0; i < nAtoms; ++i) {
      if (!currentInvariants[i]) {
        ordering.emplace_back(1, i);
      } else {
        ordering.emplace_back(0, i);
      }
    }
    std::sort(ordering.begin(), ordering.end());
    for (unsigned int i = 0; i < nAtoms; ++i) {
      atomOrder[i] = ordering[i].second;
    }
  } else {
    for (unsigned int i = 0; i < nAtoms; ++i) {
      atomOrder[i] = i;
    }
  }

  // add the round 0 invariants to the result
  for (unsigned int i = 0; i < nAtoms; ++i) {
    if (includeAtoms[i]) {
      if (!morganArguments->df_onlyNonzeroInvariants || currentInvariants[i]) {
        result.push_back(
            new MorganAtomEnv<OutputType>(currentInvariants[i], i, 0, &mol));
      }
    }
  }

  for (unsigned int layer = 0; layer < morganArguments->d_radius; ++layer) {
    std::vector<AccumTuple> allNeighborhoodsThisRound;
    for (auto atomIdx : atomOrder) {
      if (!deadAtoms[atomIdx]) {
        const std::uint32_t degree = mol.getAtomDegree(atomIdx);
        if (!degree) {
          deadAtoms.set(atomIdx, 1);
          continue;
        }
        const AtomData &tAtom = mol.getAtom(atomIdx);

        auto [bondsBegin, bondsEnd] = mol.getAtomBonds(atomIdx);
        auto [nbrsBegin, nbrsEnd] = mol.getAtomNeighbors(atomIdx);

        neighborhoodInvariants.clear();
        for (auto bondIt = bondsBegin, nbrIt = nbrsBegin; bondIt != bondsEnd;
             ++bondIt, ++nbrIt) {
          const std::uint32_t bondIdx = *bondIt;
          const std::uint32_t oIdx = *nbrIt;
          roundAtomNeighborhoods[atomIdx][bondIdx] = 1;
          roundAtomNeighborhoods[atomIdx] |= atomNeighborhoods[oIdx];

          const auto bt =
              static_cast<std::int32_t>((*bondInvariants)[bondIdx]);
          neighborhoodInvariants.emplace_back(bt, currentInvariants[oIdx]);
        }

        std::sort(neighborhoodInvariants.begin(), neighborhoodInvariants.end());

        std::uint32_t invar = layer;
        gboost::hash_combine(invar, currentInvariants[atomIdx]);
        bool looksChiral =
            (tAtom.getChiralTag() != AtomEnums::ChiralType::CHI_UNSPECIFIED);
        for (auto it = neighborhoodInvariants.begin();
             it != neighborhoodInvariants.end(); ++it) {
          gboost::hash_combine(invar, *it);
          if (morganArguments->df_includeChirality && looksChiral &&
              !chiralAtoms[atomIdx]) {
            if (it->first !=
                static_cast<std::int32_t>(BondEnums::BondType::SINGLE)) {
              looksChiral = false;
            } else if (it != neighborhoodInvariants.begin() &&
                       it->second == (it - 1)->second) {
              looksChiral = false;
            }
          }
        }

        if (morganArguments->df_includeChirality && looksChiral) {
          chiralAtoms[atomIdx] = 1;
          std::string cip;
          mol.getAtomPropIfPresent(common_properties::_CIPCodeToken, atomIdx,
                                   cip);
          if (cip == "R") {
            gboost::hash_combine(invar, 3);
          } else if (cip == "S") {
            gboost::hash_combine(invar, 2);
          } else {
            gboost::hash_combine(invar, 1);
          }
        }

        nextLayerInvariants[atomIdx] = static_cast<OutputType>(invar);
        allNeighborhoodsThisRound.emplace_back(roundAtomNeighborhoods[atomIdx],
                                               static_cast<OutputType>(invar),
                                               atomIdx);
      }
    }

    std::sort(allNeighborhoodsThisRound.begin(),
              allNeighborhoodsThisRound.end());
    for (auto iter = allNeighborhoodsThisRound.begin();
         iter != allNeighborhoodsThisRound.end(); ++iter) {
      if (morganArguments->df_includeRedundantEnvironments ||
          neighborhoods.count(std::get<0>(*iter)) == 0) {
        if (!morganArguments->df_onlyNonzeroInvariants ||
            (*atomInvariants)[std::get<2>(*iter)]) {
          if (includeAtoms[std::get<2>(*iter)]) {
            result.push_back(new MorganAtomEnv<OutputType>(
                std::get<1>(*iter), std::get<2>(*iter), layer + 1, &mol));
            neighborhoods.insert(std::get<0>(*iter));
          }
        }
      } else {
        deadAtoms[std::get<2>(*iter)] = 1;
      }
    }

    // the invariants from this round become the next round invariants:
    currentInvariants.swap(nextLayerInvariants);
    std::fill(nextLayerInvariants.begin(), nextLayerInvariants.end(), 0);

    // this rounds calculated neighbors will be next rounds initial neighbors,
    // so the radius can grow every iteration
    atomNeighborhoods = roundAtomNeighborhoods;
  }

  return result;
}

template <typename OutputType>
std::string MorganEnvGenerator<OutputType>::infoString() const {
  return "MorganEnvironmentGenerator";
}
template <typename OutputType>
void MorganEnvGenerator<OutputType>::toJSON(
    boost::property_tree::ptree &pt) const {
  pt.put("type", "MorganEnvGenerator");
  AtomEnvironmentGenerator<OutputType>::toJSON(pt);
}

template <typename OutputType>
void MorganEnvGenerator<OutputType>::fromJSON(
    const boost::property_tree::ptree &pt) {
  AtomEnvironmentGenerator<OutputType>::fromJSON(pt);
}

template <typename OutputType>
FingerprintGenerator<OutputType> *getMorganGenerator(
    const MorganArguments &args,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, bool ownsAtomInvGen,
    bool ownsBondInvGen) {
  AtomEnvironmentGenerator<OutputType> *morganEnvGenerator =
      new MorganEnvGenerator<OutputType>();

  bool ownsAtomInvGenerator = ownsAtomInvGen;
  if (!atomInvariantsGenerator) {
    atomInvariantsGenerator = new MorganAtomInvGenerator();
    ownsAtomInvGenerator = true;
  }

  bool ownsBondInvGenerator = ownsBondInvGen;
  if (!bondInvariantsGenerator) {
    bondInvariantsGenerator = new MorganBondInvGenerator(
        args.df_useBondTypes, args.df_includeChirality);
    ownsBondInvGenerator = true;
  }

  return new FingerprintGenerator<OutputType>(
      morganEnvGenerator, new MorganArguments(args), atomInvariantsGenerator,
      bondInvariantsGenerator, ownsAtomInvGenerator, ownsBondInvGenerator);
}

template <typename OutputType>
FingerprintGenerator<OutputType> *getMorganGenerator(
    unsigned int radius, bool countSimulation, bool includeChirality,
    bool useBondTypes, bool onlyNonzeroInvariants,
    bool includeRedundantEnvironments,
    AtomInvariantsGenerator *atomInvariantsGenerator,
    BondInvariantsGenerator *bondInvariantsGenerator, std::uint32_t fpSize,
    std::vector<std::uint32_t> countBounds, bool ownsAtomInvGen,
    bool ownsBondInvGen) {
  MorganArguments arguments(radius, countSimulation, includeChirality,
                            onlyNonzeroInvariants, countBounds, fpSize,
                            includeRedundantEnvironments, useBondTypes);

  return getMorganGenerator<OutputType>(arguments, atomInvariantsGenerator,
                                        bondInvariantsGenerator, ownsAtomInvGen,
                                        ownsBondInvGen);
}

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t>
    *getMorganGenerator(const MorganArguments &, AtomInvariantsGenerator *,
                        BondInvariantsGenerator *, bool, bool);
template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t>
    *getMorganGenerator(const MorganArguments &, AtomInvariantsGenerator *,
                        BondInvariantsGenerator *, bool, bool);

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint32_t> *
getMorganGenerator(unsigned int radius, bool countSimulation,
                   bool includeChirality, bool useBondTypes,
                   bool onlyNonzeroInvariants,
                   bool includeRedundantEnvironments,
                   AtomInvariantsGenerator *atomInvariantsGenerator,
                   BondInvariantsGenerator *bondInvariantsGenerator,
                   std::uint32_t fpSize, std::vector<std::uint32_t> countBounds,
                   bool ownsAtomInvGen, bool ownsBondInvGen);

template RDKIT_FINGERPRINTS_EXPORT FingerprintGenerator<std::uint64_t> *
getMorganGenerator(unsigned int radius, bool countSimulation,
                   bool includeChirality, bool useBondTypes,
                   bool onlyNonzeroInvariants,
                   bool includeRedundantEnvironments,
                   AtomInvariantsGenerator *atomInvariantsGenerator,
                   BondInvariantsGenerator *bondInvariantsGenerator,
                   std::uint32_t fpSize, std::vector<std::uint32_t> countBounds,
                   bool ownsAtomInvGen, bool ownsBondInvGen);

}  // namespace MorganFingerprint
}  // namespace RDKit
