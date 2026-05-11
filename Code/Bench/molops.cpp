#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/MolOps.h>
#include <GraphMol/RDMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

namespace {

using bench_common::Dataset;
using bench_common::dataset_name;

template <class Mol, class Op>
auto run_per_sample(const std::vector<Mol> &samples,
                    Catch::Benchmark::Chronometer meter, Op op) {
  std::vector<Mol> work;
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
      total += op(mol);
    }
    return total;
  });
}

// Sanitize step index.  Each bench primes its input by running steps
// [cleanUp .. (stage-1)] on raw parser output; the target stage then runs
// on input in the same state it would see inside sanitizeMol.
enum class SanitizeStage {
  cleanUp = 0,                  // no priming
  cleanUpOrganometallics,        // after cleanUp
  updatePropertyCache_first,     // after cleanUpOrganometallics
  symmetrizeSSSR,                // after updatePropertyCache(1)
  Kekulize,                      // after symmetrizeSSSR
  assignRadicals,                // after Kekulize
  setAromaticity,                // after assignRadicals
  setConjugation,                // after setAromaticity
  setHybridization,              // after setConjugation
  cleanupAtropisomers,           // after setHybridization
  cleanupChirality,              // after cleanupAtropisomers
  adjustHs,                      // after cleanupChirality
};

template <class Mol>
void prime_to(Mol &mol, SanitizeStage stage) {
  using S = SanitizeStage;
  if (stage <= S::cleanUp) return;
  MolOps::cleanUp(mol);
  if (stage <= S::cleanUpOrganometallics) return;
  MolOps::cleanUpOrganometallics(mol);
  if (stage <= S::updatePropertyCache_first) return;
  mol.updatePropertyCache(true);
  if (stage <= S::symmetrizeSSSR) return;
  (void)MolOps::symmetrizeSSSR(mol);
  if (stage <= S::Kekulize) return;
  MolOps::Kekulize(mol, /*markAtomsBonds=*/true, /*canonical=*/false);
  if (stage <= S::assignRadicals) return;
  MolOps::assignRadicals(mol);
  if (stage <= S::setAromaticity) return;
  (void)MolOps::setAromaticity(mol);
  if (stage <= S::setConjugation) return;
  MolOps::setConjugation(mol);
  if (stage <= S::setHybridization) return;
  MolOps::setHybridization(mol);
  if (stage <= S::cleanupAtropisomers) return;
  MolOps::cleanupAtropisomers(mol);
  if (stage <= S::cleanupChirality) return;
  MolOps::cleanupChirality(mol);
  if (stage <= S::adjustHs) return;
  MolOps::adjustHs(mol);
}

template <class Sample, class Mol>
std::vector<Mol> load_and_prime(Dataset dataset, SanitizeStage stage);

template <>
std::vector<RWMol> load_and_prime<RWMol, RWMol>(Dataset dataset,
                                                 SanitizeStage stage) {
  auto romol_samples = bench_common::load_samples(dataset, /*sanitize=*/false);
  std::vector<RWMol> primed;
  primed.reserve(romol_samples.size());
  for (auto &mol : romol_samples) {
    primed.emplace_back(mol);
    prime_to(primed.back(), stage);
  }
  return primed;
}

template <>
std::vector<RDMol> load_and_prime<RDMol, RDMol>(Dataset dataset,
                                                 SanitizeStage stage) {
  auto samples = bench_common::load_rdmol_samples(dataset, /*sanitize=*/false);
  for (auto &mol : samples) {
    prime_to(mol, stage);
  }
  return samples;
}

}  // namespace

TEST_CASE("MolOps::addHs", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolOps::addHs") {
    auto total_atoms = 0;
    for (auto &mol : samples) {
      RWMol mol_copy(mol);
      MolOps::addHs(mol_copy);
      total_atoms += mol_copy.getNumAtoms();
    }
    return total_atoms;
  };
}

TEST_CASE("MolOps::FindSSR", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolOps::FindSSR") {
    auto total = 0;
    for (auto &mol : samples) {
      total += MolOps::findSSSR(mol);
    }
    return total;
  };
}

TEST_CASE("MolOps::FindSSR RDMol", "[molops][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK("MolOps::FindSSR RDMol") {
    auto total = 0;
    RingInfoCache ringInfo;
    for (auto &mol : samples) {
      total += MolOps::findSSSR(mol, ringInfo);
    }
    return total;
  };
}

TEST_CASE("MolOps::getMolFrags", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolOps::getMolFrags") {
    auto total = 0;
    for (auto &mol : samples) {
      std::vector<std::unique_ptr<ROMol>> frags;
      MolOps::getMolFrags(mol, frags);
      for (auto &frag : frags) {
        total += frag->getNumAtoms();
      }
    }
    return total;
  };
}

// Unified sub-step bench: load raw parser output, prime to STAGE, then
// measure OP on input that matches what sanitizeMol's caller would see at
// that point.
#define BENCH_STAGE(NAME, STAGE, OP_RO, OP_RD, COUNT, DATASET, SUFFIX, TAG)    \
  TEST_CASE("MolOps::" NAME " " SUFFIX, "[molops]" TAG) {                       \
    auto rw_samples = load_and_prime<RWMol, RWMol>(DATASET, STAGE);             \
    BENCHMARK_ADVANCED("MolOps::" NAME " " SUFFIX)(                             \
        Catch::Benchmark::Chronometer meter) {                                  \
      run_per_sample(rw_samples, meter, [](RWMol &mol) {                        \
        OP_RO;                                                                  \
        return COUNT;                                                           \
      });                                                                       \
    };                                                                          \
  }                                                                             \
  TEST_CASE("MolOps::" NAME " RDMol " SUFFIX, "[molops][rdmol]" TAG) {          \
    auto samples = load_and_prime<RDMol, RDMol>(DATASET, STAGE);                \
    BENCHMARK_ADVANCED("MolOps::" NAME " RDMol " SUFFIX)(                       \
        Catch::Benchmark::Chronometer meter) {                                  \
      run_per_sample(samples, meter, [](RDMol &mol) {                           \
        OP_RD;                                                                  \
        return COUNT;                                                           \
      });                                                                       \
    };                                                                          \
  }

// Canonical (default sample set) sub-step benches with proper priming.
BENCH_STAGE("cleanUp", SanitizeStage::cleanUp, MolOps::cleanUp(mol),
            MolOps::cleanUp(mol), mol.getNumAtoms(), Dataset::Canonical, "",
            "[canonical]")
BENCH_STAGE("cleanUpOrganometallics", SanitizeStage::cleanUpOrganometallics,
            MolOps::cleanUpOrganometallics(mol),
            MolOps::cleanUpOrganometallics(mol), mol.getNumAtoms(),
            Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("updatePropertyCache_first",
            SanitizeStage::updatePropertyCache_first,
            mol.updatePropertyCache(true), mol.updatePropertyCache(true),
            mol.getNumAtoms(), Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("symmetrizeSSSR", SanitizeStage::symmetrizeSSSR,
            (void)MolOps::symmetrizeSSSR(mol),
            (void)MolOps::symmetrizeSSSR(mol), mol.getNumAtoms(),
            Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("Kekulize", SanitizeStage::Kekulize,
            MolOps::Kekulize(mol, true, false),
            MolOps::Kekulize(mol, true, false), mol.getNumBonds(),
            Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("assignRadicals", SanitizeStage::assignRadicals,
            MolOps::assignRadicals(mol), MolOps::assignRadicals(mol),
            mol.getNumAtoms(), Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("setAromaticity", SanitizeStage::setAromaticity,
            (void)MolOps::setAromaticity(mol),
            (void)MolOps::setAromaticity(mol), mol.getNumAtoms(),
            Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("setConjugation", SanitizeStage::setConjugation,
            MolOps::setConjugation(mol), MolOps::setConjugation(mol),
            mol.getNumBonds(), Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("setHybridization", SanitizeStage::setHybridization,
            MolOps::setHybridization(mol), MolOps::setHybridization(mol),
            mol.getNumAtoms(), Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("cleanupAtropisomers", SanitizeStage::cleanupAtropisomers,
            MolOps::cleanupAtropisomers(mol),
            MolOps::cleanupAtropisomers(mol), mol.getNumAtoms(),
            Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("cleanupChirality", SanitizeStage::cleanupChirality,
            MolOps::cleanupChirality(mol), MolOps::cleanupChirality(mol),
            mol.getNumAtoms(), Dataset::Canonical, "", "[canonical]")
BENCH_STAGE("adjustHs", SanitizeStage::adjustHs, MolOps::adjustHs(mol),
            MolOps::adjustHs(mol), mol.getNumAtoms(), Dataset::Canonical, "",
            "[canonical]")

TEST_CASE("MolOps::removeHs", "[molops]") {
  auto baseSamples = bench_common::load_samples();
  std::vector<RWMol> samplesWithHs;
  samplesWithHs.reserve(baseSamples.size());
  for (auto &mol : baseSamples) {
    RWMol withHs(mol);
    MolOps::addHs(withHs);
    samplesWithHs.push_back(std::move(withHs));
  }

  BENCHMARK_ADVANCED("MolOps::removeHs")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samplesWithHs, meter, [](RWMol &mol) {
      MolOps::removeHs(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::removeHs RDMol", "[molops][rdmol]") {
  auto baseSamples = bench_common::load_rdmol_samples();
  std::vector<RDMol> samplesWithHs;
  samplesWithHs.reserve(baseSamples.size());
  for (auto &mol : baseSamples) {
    RWMol withHs(mol.asROMol());
    MolOps::addHs(withHs);
    samplesWithHs.emplace_back(withHs.asRDMol());
    samplesWithHs.back().clearCompatibilityData();
  }

  BENCHMARK_ADVANCED("MolOps::removeHs RDMol")(
      Catch::Benchmark::Chronometer meter) {
    MolOps::SanitizeTemp temp;
    MolOps::RemoveHsParameters ps;
    run_per_sample(samplesWithHs, meter, [&](RDMol &mol) {
      MolOps::removeHs(mol, ps, temp, /*sanitize=*/false);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::cleanupChirality RDMol", "[molops][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("MolOps::cleanupChirality RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::cleanupChirality(mol);
      return mol.getNumAtoms();
    });
  };
}

// ---------------------------------------------------------------------------
// Umbrella sanitize / removeHs benches.  These reflect what SmilesToMol
// actually calls on the default parse path (removeHs(mol, ps, temp,
// sanitize=true) which internally invokes sanitizeMol).  The standalone
// sub-step benches above are not additive with the umbrella benches: they
// measure the same code paths in isolation, but in real pipeline use they
// run inside sanitizeMol.

namespace {

std::vector<RDMol> load_unsanitized_rdmol_samples(Dataset dataset) {
  auto samples = bench_common::load_rdmol_samples(dataset, /*sanitize=*/false);
  // For the umbrella sanitize/removeHs bench we want an input shape that
  // matches what SmilesToMol passes to removeHs/sanitizeMol on the default
  // path -- i.e. fresh parser output, no sanitize done yet.
  return samples;
}

std::vector<ROMol> load_unsanitized_romol_samples(Dataset dataset) {
  return bench_common::load_samples(dataset, /*sanitize=*/false);
}

}  // namespace

TEST_CASE("MolOps::sanitizeMol", "[molops]") {
  auto samples = load_unsanitized_romol_samples(Dataset::Canonical);
  std::vector<RWMol> rw_samples;
  rw_samples.reserve(samples.size());
  for (auto &mol : samples) {
    rw_samples.emplace_back(mol);
  }
  BENCHMARK_ADVANCED("MolOps::sanitizeMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::sanitizeMol(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::sanitizeMol RDMol", "[molops][rdmol]") {
  auto samples = load_unsanitized_rdmol_samples(Dataset::Canonical);
  BENCHMARK_ADVANCED("MolOps::sanitizeMol RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::sanitizeMol(mol);
      return mol.getNumAtoms();
    });
  };
}

#define BENCH_SANITIZE(DATASET, SUFFIX, TAG)                                   \
  TEST_CASE("MolOps::sanitizeMol " SUFFIX, "[molops]" TAG) {                    \
    auto samples = load_unsanitized_romol_samples(DATASET);                     \
    std::vector<RWMol> rw_samples;                                              \
    rw_samples.reserve(samples.size());                                         \
    for (auto &mol : samples) {                                                 \
      rw_samples.emplace_back(mol);                                             \
    }                                                                           \
    BENCHMARK_ADVANCED("MolOps::sanitizeMol " SUFFIX)(                          \
        Catch::Benchmark::Chronometer meter) {                                  \
      run_per_sample(rw_samples, meter, [](RWMol &mol) {                        \
        MolOps::sanitizeMol(mol);                                               \
        return mol.getNumAtoms();                                               \
      });                                                                       \
    };                                                                          \
  }                                                                             \
  TEST_CASE("MolOps::sanitizeMol RDMol " SUFFIX, "[molops][rdmol]" TAG) {       \
    auto samples = load_unsanitized_rdmol_samples(DATASET);                     \
    BENCHMARK_ADVANCED("MolOps::sanitizeMol RDMol " SUFFIX)(                    \
        Catch::Benchmark::Chronometer meter) {                                  \
      run_per_sample(samples, meter, [](RDMol &mol) {                           \
        MolOps::sanitizeMol(mol);                                               \
        return mol.getNumAtoms();                                               \
      });                                                                       \
    };                                                                          \
  }

BENCH_SANITIZE(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_SANITIZE(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_SANITIZE(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_SANITIZE(Dataset::Size_60_80, "size 60-80", "[size_60_80]")
BENCH_SANITIZE(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_SANITIZE(Dataset::KekulizeHard, "kekulize_hard", "[kekulize_hard]")

// removeHs core (sanitize=false): the H-removal pass without the recursive
// sanitizeMol that SmilesToMol's default path triggers.  This is what runs
// in the pipeline BEFORE sanitize is invoked.

#define BENCH_REMOVEHS_CORE(DATASET, SUFFIX, TAG)                              \
  TEST_CASE("MolOps::removeHs core " SUFFIX, "[molops]" TAG) {                  \
    auto romol_samples = load_unsanitized_romol_samples(DATASET);               \
    std::vector<RWMol> rw_samples;                                              \
    rw_samples.reserve(romol_samples.size());                                   \
    for (auto &mol : romol_samples) {                                           \
      rw_samples.emplace_back(mol);                                             \
    }                                                                           \
    MolOps::RemoveHsParameters ps;                                              \
    BENCHMARK_ADVANCED("MolOps::removeHs core " SUFFIX)(                        \
        Catch::Benchmark::Chronometer meter) {                                  \
      run_per_sample(rw_samples, meter, [&ps](RWMol &mol) {                     \
        MolOps::removeHs(mol, ps, /*sanitize=*/false);                          \
        return mol.getNumAtoms();                                               \
      });                                                                       \
    };                                                                          \
  }                                                                             \
  TEST_CASE("MolOps::removeHs core RDMol " SUFFIX, "[molops][rdmol]" TAG) {     \
    auto samples = load_unsanitized_rdmol_samples(DATASET);                     \
    MolOps::RemoveHsParameters ps;                                              \
    BENCHMARK_ADVANCED("MolOps::removeHs core RDMol " SUFFIX)(                  \
        Catch::Benchmark::Chronometer meter) {                                  \
      MolOps::SanitizeTemp temp;                                                \
      run_per_sample(samples, meter, [&ps, &temp](RDMol &mol) {                 \
        MolOps::removeHs(mol, ps, temp, /*sanitize=*/false);                    \
        return mol.getNumAtoms();                                               \
      });                                                                       \
    };                                                                          \
  }

BENCH_REMOVEHS_CORE(Dataset::Canonical, "", "[canonical]")
BENCH_REMOVEHS_CORE(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_REMOVEHS_CORE(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_REMOVEHS_CORE(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_REMOVEHS_CORE(Dataset::Size_60_80, "size 60-80", "[size_60_80]")
BENCH_REMOVEHS_CORE(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_REMOVEHS_CORE(Dataset::KekulizeHard, "kekulize_hard", "[kekulize_hard]")

// ---------------------------------------------------------------------------
// Per-bucket size + ring-count variants. Each existing benchmark above gets a
// matching TEST_CASE for every size bucket (4) and ring-count bucket (5),
// tagged with `[size_*]` / `[rings_*]` so they can be selected individually.
//
// The bodies are macro-generated to keep this file scannable. The macro
// expands into:
//   TEST_CASE("OP <suffix>", "[molops]<tag>") { ... }
// using `LoaderRO` for ROMol-side ops and `LoaderRD` for RDMol-side ops.

#define BENCH_ROMOL_ITER_BODY(OP_NAME, OP_EXPR, COUNT_EXPR, DATASET, SUFFIX,   \
                              TAG)                                             \
  TEST_CASE("MolOps::" OP_NAME " " SUFFIX, "[molops]" TAG) {                   \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolOps::" OP_NAME " " SUFFIX) {                                 \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        OP_EXPR;                                                               \
        total += COUNT_EXPR;                                                   \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

// addHs (operates on RWMol copy; bench iterates ROMol samples)
#define BENCH_ADDHS(DATASET, SUFFIX, TAG)                                      \
  TEST_CASE("MolOps::addHs " SUFFIX, "[molops]" TAG) {                         \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolOps::addHs " SUFFIX) {                                       \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        RWMol mol_copy(mol);                                                   \
        MolOps::addHs(mol_copy);                                               \
        total += mol_copy.getNumAtoms();                                       \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_FINDSSR(DATASET, SUFFIX, TAG)                                    \
  TEST_CASE("MolOps::FindSSR " SUFFIX, "[molops]" TAG) {                       \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolOps::FindSSR " SUFFIX) {                                     \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        total += MolOps::findSSSR(mol);                                        \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_FINDSSR_RD(DATASET, SUFFIX, TAG)                                 \
  TEST_CASE("MolOps::FindSSR RDMol " SUFFIX, "[molops][rdmol]" TAG) {          \
    auto samples = bench_common::load_rdmol_samples(DATASET);                  \
    BENCHMARK("MolOps::FindSSR RDMol " SUFFIX) {                               \
      auto total = 0;                                                          \
      RingInfoCache ringInfo;                                                  \
      for (auto &mol : samples) {                                              \
        total += MolOps::findSSSR(mol, ringInfo);                              \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_GETMOLFRAGS(DATASET, SUFFIX, TAG)                                \
  TEST_CASE("MolOps::getMolFrags " SUFFIX, "[molops]" TAG) {                   \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolOps::getMolFrags " SUFFIX) {                                 \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        std::vector<std::unique_ptr<ROMol>> frags;                             \
        MolOps::getMolFrags(mol, frags);                                       \
        for (auto &frag : frags) {                                             \
          total += frag->getNumAtoms();                                        \
        }                                                                      \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

// Macro for ROMol-template ops that take Chronometer + run_per_sample.
// The body always copies the ROMol samples into RWMol so OP can mutate.
#define BENCH_RO_ADV(NAME, OP, COUNT, DATASET, SUFFIX, TAG)                    \
  TEST_CASE("MolOps::" NAME " " SUFFIX, "[molops]" TAG) {                      \
    auto romol_samples = bench_common::load_samples(DATASET);                  \
    std::vector<RWMol> rw_samples;                                             \
    rw_samples.reserve(romol_samples.size());                                  \
    for (auto &mol : romol_samples) {                                          \
      rw_samples.emplace_back(mol);                                            \
    }                                                                          \
    BENCHMARK_ADVANCED("MolOps::" NAME " " SUFFIX)(                            \
        Catch::Benchmark::Chronometer meter) {                                 \
      run_per_sample(rw_samples, meter, [](RWMol &mol) {                       \
        OP;                                                                    \
        return COUNT;                                                          \
      });                                                                      \
    };                                                                         \
  }

// Macro for RDMol-template ops.
#define BENCH_RD_ADV(NAME, OP, COUNT, DATASET, SUFFIX, TAG)                    \
  TEST_CASE("MolOps::" NAME " RDMol " SUFFIX, "[molops][rdmol]" TAG) {         \
    auto samples = bench_common::load_rdmol_samples(DATASET);                  \
    BENCHMARK_ADVANCED("MolOps::" NAME " RDMol " SUFFIX)(                      \
        Catch::Benchmark::Chronometer meter) {                                 \
      run_per_sample(samples, meter, [](RDMol &mol) {                          \
        OP;                                                                    \
        return COUNT;                                                          \
      });                                                                      \
    };                                                                         \
  }

// Per-bucket variants for non-sub-step utilities (addHs / FindSSR /
// getMolFrags) only. Sub-step benches (cleanUp, Kekulize, etc.) come from
// BENCH_STAGE_FOR_BUCKETS below, which primes inputs to the correct state.
#define BENCH_OPS_FOR(DATASET, SUFFIX, TAG)                                    \
  BENCH_ADDHS(DATASET, SUFFIX, TAG)                                            \
  BENCH_FINDSSR(DATASET, SUFFIX, TAG)                                          \
  BENCH_FINDSSR_RD(DATASET, SUFFIX, TAG)                                       \
  BENCH_GETMOLFRAGS(DATASET, SUFFIX, TAG)

BENCH_OPS_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_OPS_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_OPS_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_OPS_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_OPS_FOR(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_OPS_FOR(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_OPS_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_OPS_FOR(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_OPS_FOR(Dataset::Rings_6, "rings 6", "[rings_6]")

// Per-bucket sub-step benches using BENCH_STAGE for correct priming.
#define BENCH_STAGE_FOR_BUCKETS(DATASET, SUFFIX, TAG)                          \
  BENCH_STAGE("cleanUp", SanitizeStage::cleanUp, MolOps::cleanUp(mol),         \
              MolOps::cleanUp(mol), mol.getNumAtoms(), DATASET, SUFFIX, TAG)   \
  BENCH_STAGE("cleanUpOrganometallics",                                       \
              SanitizeStage::cleanUpOrganometallics,                           \
              MolOps::cleanUpOrganometallics(mol),                             \
              MolOps::cleanUpOrganometallics(mol), mol.getNumAtoms(), DATASET, \
              SUFFIX, TAG)                                                     \
  BENCH_STAGE("updatePropertyCache_first",                                    \
              SanitizeStage::updatePropertyCache_first,                        \
              mol.updatePropertyCache(true), mol.updatePropertyCache(true),    \
              mol.getNumAtoms(), DATASET, SUFFIX, TAG)                         \
  BENCH_STAGE("symmetrizeSSSR", SanitizeStage::symmetrizeSSSR,                \
              (void)MolOps::symmetrizeSSSR(mol),                               \
              (void)MolOps::symmetrizeSSSR(mol), mol.getNumAtoms(), DATASET,   \
              SUFFIX, TAG)                                                     \
  BENCH_STAGE("Kekulize", SanitizeStage::Kekulize,                            \
              MolOps::Kekulize(mol, true, false),                              \
              MolOps::Kekulize(mol, true, false), mol.getNumBonds(), DATASET,  \
              SUFFIX, TAG)                                                     \
  BENCH_STAGE("assignRadicals", SanitizeStage::assignRadicals,                \
              MolOps::assignRadicals(mol), MolOps::assignRadicals(mol),        \
              mol.getNumAtoms(), DATASET, SUFFIX, TAG)                         \
  BENCH_STAGE("setAromaticity", SanitizeStage::setAromaticity,                \
              (void)MolOps::setAromaticity(mol),                               \
              (void)MolOps::setAromaticity(mol), mol.getNumAtoms(), DATASET,   \
              SUFFIX, TAG)                                                     \
  BENCH_STAGE("setConjugation", SanitizeStage::setConjugation,                \
              MolOps::setConjugation(mol), MolOps::setConjugation(mol),        \
              mol.getNumBonds(), DATASET, SUFFIX, TAG)                         \
  BENCH_STAGE("setHybridization", SanitizeStage::setHybridization,            \
              MolOps::setHybridization(mol), MolOps::setHybridization(mol),    \
              mol.getNumAtoms(), DATASET, SUFFIX, TAG)                         \
  BENCH_STAGE("cleanupAtropisomers", SanitizeStage::cleanupAtropisomers,      \
              MolOps::cleanupAtropisomers(mol),                                \
              MolOps::cleanupAtropisomers(mol), mol.getNumAtoms(), DATASET,    \
              SUFFIX, TAG)                                                     \
  BENCH_STAGE("cleanupChirality", SanitizeStage::cleanupChirality,            \
              MolOps::cleanupChirality(mol), MolOps::cleanupChirality(mol),    \
              mol.getNumAtoms(), DATASET, SUFFIX, TAG)                         \
  BENCH_STAGE("adjustHs", SanitizeStage::adjustHs, MolOps::adjustHs(mol),     \
              MolOps::adjustHs(mol), mol.getNumAtoms(), DATASET, SUFFIX, TAG)

BENCH_STAGE_FOR_BUCKETS(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_STAGE_FOR_BUCKETS(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_STAGE_FOR_BUCKETS(Dataset::Size_60_80, "size 60-80", "[size_60_80]")
BENCH_STAGE_FOR_BUCKETS(Dataset::Rings_4, "rings 4", "[rings_4]")

// ---------------------------------------------------------------------------
// Edge-case benches. These target inner branches in cleanUp / assignRadicals /
// Kekulize that the canonical samples don't reach. The pre-canonical inputs
// (organometallics, hypervalent_*, pre_canonical_no2_azide) are intentionally
// not sanitize-clean and load with sanitize=false so cleanUp itself is what
// rewrites the affected atoms.

namespace {

template <bool Sanitize>
void run_cleanup_pre_canonical(Dataset dataset, const char *label) {
  auto romol_samples = bench_common::load_samples(dataset, Sanitize);
  std::vector<RWMol> rw_samples;
  rw_samples.reserve(romol_samples.size());
  for (auto &mol : romol_samples) {
    rw_samples.emplace_back(mol);
  }
  BENCHMARK_ADVANCED(label)(Catch::Benchmark::Chronometer meter) {
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::cleanUp(mol);
      return mol.getNumAtoms();
    });
  };
}

template <bool Sanitize>
void run_cleanup_pre_canonical_rdmol(Dataset dataset, const char *label) {
  auto samples = bench_common::load_rdmol_samples(dataset, Sanitize);
  BENCHMARK_ADVANCED(label)(Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::cleanUp(mol);
      return mol.getNumAtoms();
    });
  };
}

}  // namespace

TEST_CASE("MolOps::cleanUp organometallics", "[molops][organometallics]") {
  run_cleanup_pre_canonical<false>(Dataset::Organometallics,
                                   "MolOps::cleanUp organometallics");
}
TEST_CASE("MolOps::cleanUp organometallics RDMol",
          "[molops][rdmol][organometallics]") {
  run_cleanup_pre_canonical_rdmol<false>(Dataset::Organometallics,
                                         "MolOps::cleanUp organometallics RDMol");
}

TEST_CASE("MolOps::cleanUp hypervalent_halogens",
          "[molops][hypervalent_halogens]") {
  run_cleanup_pre_canonical<false>(Dataset::HypervalentHalogens,
                                   "MolOps::cleanUp hypervalent_halogens");
}
TEST_CASE("MolOps::cleanUp hypervalent_halogens RDMol",
          "[molops][rdmol][hypervalent_halogens]") {
  run_cleanup_pre_canonical_rdmol<false>(
      Dataset::HypervalentHalogens, "MolOps::cleanUp hypervalent_halogens RDMol");
}

TEST_CASE("MolOps::cleanUp hypervalent_p", "[molops][hypervalent_p]") {
  run_cleanup_pre_canonical<false>(Dataset::HypervalentP,
                                   "MolOps::cleanUp hypervalent_p");
}
TEST_CASE("MolOps::cleanUp hypervalent_p RDMol",
          "[molops][rdmol][hypervalent_p]") {
  run_cleanup_pre_canonical_rdmol<false>(Dataset::HypervalentP,
                                         "MolOps::cleanUp hypervalent_p RDMol");
}

TEST_CASE("MolOps::cleanUp pre_canonical_no2_azide",
          "[molops][pre_canonical_no2_azide]") {
  run_cleanup_pre_canonical<false>(Dataset::PreCanonicalNO2Azide,
                                   "MolOps::cleanUp pre_canonical_no2_azide");
}
TEST_CASE("MolOps::cleanUp pre_canonical_no2_azide RDMol",
          "[molops][rdmol][pre_canonical_no2_azide]") {
  run_cleanup_pre_canonical_rdmol<false>(
      Dataset::PreCanonicalNO2Azide,
      "MolOps::cleanUp pre_canonical_no2_azide RDMol");
}

TEST_CASE("MolOps::assignRadicals radicals", "[molops][radicals]") {
  auto samples = bench_common::load_samples(Dataset::Radicals);
  std::vector<RWMol> rw_samples;
  rw_samples.reserve(samples.size());
  for (auto &mol : samples) {
    rw_samples.emplace_back(mol);
  }
  BENCHMARK_ADVANCED("MolOps::assignRadicals radicals")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::assignRadicals(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::assignRadicals radicals RDMol",
          "[molops][rdmol][radicals]") {
  auto samples = bench_common::load_rdmol_samples(Dataset::Radicals);
  BENCHMARK_ADVANCED("MolOps::assignRadicals radicals RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::assignRadicals(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::Kekulize kekulize_hard", "[molops][kekulize_hard]") {
  auto samples = bench_common::load_samples(Dataset::KekulizeHard);
  std::vector<RWMol> rw_samples;
  rw_samples.reserve(samples.size());
  for (auto &mol : samples) {
    rw_samples.emplace_back(mol);
  }
  BENCHMARK_ADVANCED("MolOps::Kekulize kekulize_hard")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::Kekulize(mol, /*markAtomsBonds=*/true, /*canonical=*/false);
      return mol.getNumBonds();
    });
  };
}

TEST_CASE("MolOps::Kekulize kekulize_hard RDMol",
          "[molops][rdmol][kekulize_hard]") {
  auto samples = bench_common::load_rdmol_samples(Dataset::KekulizeHard);
  BENCHMARK_ADVANCED("MolOps::Kekulize kekulize_hard RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::Kekulize(mol, /*markAtomsBonds=*/true, /*canonical=*/false);
      return mol.getNumBonds();
    });
  };
}
