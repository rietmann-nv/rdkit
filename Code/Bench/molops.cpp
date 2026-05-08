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

TEST_CASE("MolOps::setConjugation", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::setConjugation")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](ROMol &mol) {
      MolOps::setConjugation(mol);
      return mol.getNumBonds();
    });
  };
}

TEST_CASE("MolOps::setConjugation RDMol", "[molops][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("MolOps::setConjugation RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::setConjugation(mol);
      return mol.getNumBonds();
    });
  };
}

TEST_CASE("MolOps::setHybridization", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::setHybridization")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](ROMol &mol) {
      MolOps::setHybridization(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::setHybridization RDMol", "[molops][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("MolOps::setHybridization RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::setHybridization(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::adjustHs", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::adjustHs")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> rw_samples;
    rw_samples.reserve(samples.size());
    for (auto &mol : samples) {
      rw_samples.emplace_back(mol);
    }
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::adjustHs(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::adjustHs RDMol", "[molops][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("MolOps::adjustHs RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::adjustHs(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::assignRadicals", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::assignRadicals")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> rw_samples;
    rw_samples.reserve(samples.size());
    for (auto &mol : samples) {
      rw_samples.emplace_back(mol);
    }
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::assignRadicals(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::assignRadicals RDMol", "[molops][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("MolOps::assignRadicals RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::assignRadicals(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::cleanUp", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::cleanUp")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> rw_samples;
    rw_samples.reserve(samples.size());
    for (auto &mol : samples) {
      rw_samples.emplace_back(mol);
    }
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::cleanUp(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::cleanUp RDMol", "[molops][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("MolOps::cleanUp RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::cleanUp(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::cleanupChirality", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::cleanupChirality")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> rw_samples;
    rw_samples.reserve(samples.size());
    for (auto &mol : samples) {
      rw_samples.emplace_back(mol);
    }
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::cleanupChirality(mol);
      return mol.getNumAtoms();
    });
  };
}

TEST_CASE("MolOps::Kekulize", "[molops]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("MolOps::Kekulize")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<RWMol> rw_samples;
    rw_samples.reserve(samples.size());
    for (auto &mol : samples) {
      rw_samples.emplace_back(mol);
    }
    run_per_sample(rw_samples, meter, [](RWMol &mol) {
      MolOps::Kekulize(mol, /*markAtomsBonds=*/true,
                       /*canonical=*/false);
      return mol.getNumBonds();
    });
  };
}

TEST_CASE("MolOps::Kekulize RDMol", "[molops][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("MolOps::Kekulize RDMol")(
      Catch::Benchmark::Chronometer meter) {
    run_per_sample(samples, meter, [](RDMol &mol) {
      MolOps::Kekulize(mol, /*markAtomsBonds=*/true,
                       /*canonical=*/false);
      return mol.getNumBonds();
    });
  };
}

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

#define BENCH_OPS_FOR(DATASET, SUFFIX, TAG)                                    \
  BENCH_ADDHS(DATASET, SUFFIX, TAG)                                            \
  BENCH_FINDSSR(DATASET, SUFFIX, TAG)                                          \
  BENCH_FINDSSR_RD(DATASET, SUFFIX, TAG)                                       \
  BENCH_GETMOLFRAGS(DATASET, SUFFIX, TAG)                                      \
  BENCH_RO_ADV("setConjugation", MolOps::setConjugation(mol),                  \
               mol.getNumBonds(), DATASET, SUFFIX, TAG)                        \
  BENCH_RD_ADV("setConjugation", MolOps::setConjugation(mol),                  \
               mol.getNumBonds(), DATASET, SUFFIX, TAG)                        \
  BENCH_RO_ADV("setHybridization", MolOps::setHybridization(mol),              \
               mol.getNumAtoms(), DATASET, SUFFIX, TAG)                        \
  BENCH_RD_ADV("setHybridization", MolOps::setHybridization(mol),              \
               mol.getNumAtoms(), DATASET, SUFFIX, TAG)                        \
  BENCH_RO_ADV("adjustHs", MolOps::adjustHs(mol), mol.getNumAtoms(), DATASET,  \
               SUFFIX, TAG)                                                    \
  BENCH_RD_ADV("adjustHs", MolOps::adjustHs(mol), mol.getNumAtoms(), DATASET,  \
               SUFFIX, TAG)                                                    \
  BENCH_RO_ADV("assignRadicals", MolOps::assignRadicals(mol),                  \
               mol.getNumAtoms(), DATASET, SUFFIX, TAG)                        \
  BENCH_RD_ADV("assignRadicals", MolOps::assignRadicals(mol),                  \
               mol.getNumAtoms(), DATASET, SUFFIX, TAG)                        \
  BENCH_RO_ADV("cleanupChirality", MolOps::cleanupChirality(mol),              \
               mol.getNumAtoms(), DATASET, SUFFIX, TAG)                        \
  BENCH_RD_ADV("cleanupChirality", MolOps::cleanupChirality(mol),              \
               mol.getNumAtoms(), DATASET, SUFFIX, TAG)                        \
  BENCH_RO_ADV("Kekulize", MolOps::Kekulize(mol, true, false),                 \
               mol.getNumBonds(), DATASET, SUFFIX, TAG)                        \
  BENCH_RD_ADV("Kekulize", MolOps::Kekulize(mol, true, false),                 \
               mol.getNumBonds(), DATASET, SUFFIX, TAG)

BENCH_OPS_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_OPS_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_OPS_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_OPS_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_OPS_FOR(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_OPS_FOR(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_OPS_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_OPS_FOR(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_OPS_FOR(Dataset::Rings_6, "rings 6", "[rings_6]")

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
