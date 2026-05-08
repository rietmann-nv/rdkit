#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/test_fixtures.h>

using namespace RDKit;

using bench_common::Dataset;

TEST_CASE("Chirality::findPotentialStereo", "[stereo]") {
  auto samples = bench_common::load_samples();

  BENCHMARK("Chirality::findPotentialStereo") {
    auto total = 0;

    for (auto &mol : samples) {
      auto stereo_infos = Chirality::findPotentialStereo(mol);

      // workaround for https://github.com/rdkit/rdkit/issues/8880
      mol.clearComputedProps();

      for (auto &info : stereo_infos) {
        total += info.controllingAtoms.size();
      }
    }

    return total;
  };
}

TEST_CASE("CIPLabeler::assignCIPLabels", "[stereo]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("CIPLabeler::assignCIPLabels") {
    for (auto &mol : samples) {
      CIPLabeler::assignCIPLabels(mol);
    }
  };
}

TEST_CASE("MolOps::clearSingleBondDirFlags", "[stereo]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::clearSingleBondDirFlags")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<ROMol> work;
    work.reserve(meter.runs() * samples.size());
    for (int run = 0; run < meter.runs(); ++run) {
      for (const auto &mol : samples) {
        work.emplace_back(mol);
      }
    }
    meter.measure([&](int i) {
      uint64_t total = 0;
      for (size_t s = 0; s < samples.size(); ++s) {
        auto &mol = work[i * samples.size() + s];
        MolOps::clearSingleBondDirFlags(mol);
        total += mol.getNumBonds();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::clearSingleBondDirFlags RDMol", "[stereo][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("MolOps::clearSingleBondDirFlags RDMol")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RDMol> work;
    work.reserve(meter.runs() * samples.size());
    for (int run = 0; run < meter.runs(); ++run) {
      for (const auto &mol : samples) {
        work.emplace_back(mol);
      }
    }
    meter.measure([&](int i) {
      uint64_t total = 0;
      for (size_t s = 0; s < samples.size(); ++s) {
        auto &mol = work[i * samples.size() + s];
        MolOps::clearSingleBondDirFlags(mol);
        total += mol.getNumBonds();
      }
      return total;
    });
  };
}

TEST_CASE("MolOps::assignStereochemistry", "[stereo]") {
  const auto cleanIt = true;
  const auto force = true;
  const auto flagPossibleStereoCenters = true;

  auto samples = bench_common::load_samples();

  const auto legacy = GENERATE(true, false);
  UseLegacyStereoPerceptionFixture fx(legacy);
  auto str_legacy = std::string(legacy ? "true" : "false");

  BENCHMARK("MolOps::assignStereochemistry legacy=" + str_legacy) {
    auto total = 0;

    for (auto &mol : samples) {
      MolOps::assignStereochemistry(mol, cleanIt, force,
                                    flagPossibleStereoCenters);
      for (auto atom : mol.atoms()) {
        total += atom->getChiralTag();
      }

      // workaround for https://github.com/rdkit/rdkit/issues/8880
      mol.clearComputedProps();
    }

    return total;
  };
}

// ---------------------------------------------------------------------------
// Per-bucket size + ring-count variants.

#define BENCH_FIND_POT_STEREO(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("Chirality::findPotentialStereo " SUFFIX, "[stereo]" TAG) {        \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("Chirality::findPotentialStereo " SUFFIX) {                      \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        auto stereo_infos = Chirality::findPotentialStereo(mol);               \
        mol.clearComputedProps();                                              \
        for (auto &info : stereo_infos) {                                      \
          total += info.controllingAtoms.size();                               \
        }                                                                      \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_CIP_LABELS(DATASET, SUFFIX, TAG)                                 \
  TEST_CASE("CIPLabeler::assignCIPLabels " SUFFIX, "[stereo]" TAG) {           \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("CIPLabeler::assignCIPLabels " SUFFIX) {                         \
      for (auto &mol : samples) {                                              \
        CIPLabeler::assignCIPLabels(mol);                                      \
      }                                                                        \
    };                                                                         \
  }

#define BENCH_CLEAR_BOND_DIR(DATASET, SUFFIX, TAG)                             \
  TEST_CASE("MolOps::clearSingleBondDirFlags " SUFFIX, "[stereo]" TAG) {       \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK_ADVANCED("MolOps::clearSingleBondDirFlags " SUFFIX)(             \
        Catch::Benchmark::Chronometer meter) {                                 \
      std::vector<ROMol> work;                                                 \
      work.reserve(meter.runs() * samples.size());                             \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (const auto &mol : samples) {                                      \
          work.emplace_back(mol);                                              \
        }                                                                      \
      }                                                                        \
      meter.measure([&](int i) {                                               \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          auto &mol = work[i * samples.size() + s];                            \
          MolOps::clearSingleBondDirFlags(mol);                                \
          total += mol.getNumBonds();                                          \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_CLEAR_BOND_DIR_RD(DATASET, SUFFIX, TAG)                          \
  TEST_CASE("MolOps::clearSingleBondDirFlags RDMol " SUFFIX,                   \
            "[stereo][rdmol]" TAG) {                                           \
    auto samples = bench_common::load_rdmol_samples(DATASET);                  \
    BENCHMARK_ADVANCED("MolOps::clearSingleBondDirFlags RDMol " SUFFIX)(       \
        Catch::Benchmark::Chronometer meter) {                                 \
      std::vector<RDMol> work;                                                 \
      work.reserve(meter.runs() * samples.size());                             \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (const auto &mol : samples) {                                      \
          work.emplace_back(mol);                                              \
        }                                                                      \
      }                                                                        \
      meter.measure([&](int i) {                                               \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          auto &mol = work[i * samples.size() + s];                            \
          MolOps::clearSingleBondDirFlags(mol);                                \
          total += mol.getNumBonds();                                          \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_ASSIGN_STEREO(DATASET, SUFFIX, TAG)                              \
  TEST_CASE("MolOps::assignStereochemistry " SUFFIX, "[stereo]" TAG) {         \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolOps::assignStereochemistry " SUFFIX) {                       \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        MolOps::assignStereochemistry(mol, true, true, true);                  \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->getChiralTag();                                       \
        }                                                                      \
        mol.clearComputedProps();                                              \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_STEREO_FOR(DATASET, SUFFIX, TAG)                                 \
  BENCH_FIND_POT_STEREO(DATASET, SUFFIX, TAG)                                  \
  BENCH_CIP_LABELS(DATASET, SUFFIX, TAG)                                       \
  BENCH_CLEAR_BOND_DIR(DATASET, SUFFIX, TAG)                                   \
  BENCH_CLEAR_BOND_DIR_RD(DATASET, SUFFIX, TAG)                                \
  BENCH_ASSIGN_STEREO(DATASET, SUFFIX, TAG)

BENCH_STEREO_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_STEREO_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_STEREO_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_STEREO_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_STEREO_FOR(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_STEREO_FOR(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_STEREO_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_STEREO_FOR(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_STEREO_FOR(Dataset::Rings_6, "rings 6", "[rings_6]")

// Atropisomers carry the wedge-bond / atrop stereo markers that
// findPotentialStereo and assignStereochemistry have to track.
BENCH_FIND_POT_STEREO(Dataset::Atropisomers, "atropisomers", "[atropisomers]")
BENCH_ASSIGN_STEREO(Dataset::Atropisomers, "atropisomers", "[atropisomers]")
