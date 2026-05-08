#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/RDMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

using bench_common::Dataset;

TEST_CASE("ROMol copy constructor", "[mol]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("ROMol copy constructor")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<Catch::Benchmark::storage_for<ROMol>> storage(meter.runs() *
                                                              samples.size());
    meter.measure([&](int i) {
      for (size_t sample = 0; sample < samples.size(); ++sample) {
        storage[i * samples.size() + sample].construct(samples[sample]);
      }
    });
  };
}

TEST_CASE("ROMol destructor", "[mol]") {
  auto samples = bench_common::load_samples();
  BENCHMARK_ADVANCED("ROMol destructor")(Catch::Benchmark::Chronometer meter) {
    std::vector<Catch::Benchmark::destructable_object<ROMol>> storage(
        meter.runs() * samples.size());
    for (size_t i = 0; i < storage.size(); ++i) {
      storage[i].construct(samples[i % samples.size()]);
    }
    meter.measure([&](int i) {
      for (size_t sample = 0; sample < samples.size(); ++sample) {
        storage[i * samples.size() + sample].destruct();
      }
    });
  };
}

TEST_CASE("memory pressure test", "[mol][size]") {
  auto cases = bench_common::load_samples();
  REQUIRE(!cases.empty());

  const size_t N = 10400;
  std::vector<ROMol> mols;
  mols.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    mols.emplace_back(cases[i % cases.size()]);
  }
  REQUIRE(mols.size() == N);

  BENCHMARK("memory pressure test", i) {
    // copy from one random location to another
    auto a = bench_common::nth_random(i);
    auto src_idx = a % mols.size();
    auto dst_idx = (a / mols.size()) % mols.size();
    ROMol temp(mols[src_idx]);
    mols[dst_idx] = std::move(temp);
  };
}

TEST_CASE("ROMol::getNumHeavyAtoms", "[mol]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("ROMol::getNumHeavyAtoms") {
    auto sum = 0;
    for (auto &mol : samples) {
      sum += mol.getNumHeavyAtoms();
    }
    return sum;
  };
}

TEST_CASE("RDMol copy constructor", "[mol][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("RDMol copy constructor")(
      Catch::Benchmark::Chronometer meter) {
    std::vector<Catch::Benchmark::storage_for<RDMol>> storage(meter.runs() *
                                                              samples.size());
    meter.measure([&](int i) {
      for (size_t sample = 0; sample < samples.size(); ++sample) {
        storage[i * samples.size() + sample].construct(samples[sample]);
      }
    });
  };
}

TEST_CASE("RDMol destructor", "[mol][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK_ADVANCED("RDMol destructor")(Catch::Benchmark::Chronometer meter) {
    std::vector<Catch::Benchmark::destructable_object<RDMol>> storage(
        meter.runs() * samples.size());
    for (size_t i = 0; i < storage.size(); ++i) {
      storage[i].construct(samples[i % samples.size()]);
    }
    meter.measure([&](int i) {
      for (size_t sample = 0; sample < samples.size(); ++sample) {
        storage[i * samples.size() + sample].destruct();
      }
    });
  };
}

TEST_CASE("RDMol memory pressure test", "[mol][rdmol][size]") {
  auto cases = bench_common::load_rdmol_samples();
  REQUIRE(!cases.empty());

  const size_t N = 10400;
  std::vector<RDMol> mols;
  mols.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    mols.emplace_back(cases[i % cases.size()]);
  }
  REQUIRE(mols.size() == N);

  BENCHMARK("RDMol memory pressure test", i) {
    auto a = bench_common::nth_random(i);
    auto src_idx = a % mols.size();
    auto dst_idx = (a / mols.size()) % mols.size();
    RDMol temp(mols[src_idx]);
    mols[dst_idx] = std::move(temp);
  };
}

TEST_CASE("RDMol::getNumHeavyAtoms", "[mol][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  BENCHMARK("RDMol::getNumHeavyAtoms") {
    auto sum = 0;
    for (auto &mol : samples) {
      sum += mol.getNumHeavyAtoms();
    }
    return sum;
  };
}

// ---------------------------------------------------------------------------
// Per-bucket size + ring-count variants. Memory pressure tests are NOT
// duplicated -- they're parameter-tuned for the canonical set's mol size and
// would need separate sizing per bucket; if/when needed, they can be added
// later.

#define BENCH_ROMOL_COPY(DATASET, SUFFIX, TAG)                                 \
  TEST_CASE("ROMol copy constructor " SUFFIX, "[mol]" TAG) {                   \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK_ADVANCED("ROMol copy constructor " SUFFIX)(                      \
        Catch::Benchmark::Chronometer meter) {                                 \
      std::vector<Catch::Benchmark::storage_for<ROMol>> storage(               \
          meter.runs() * samples.size());                                      \
      meter.measure([&](int i) {                                               \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          storage[i * samples.size() + sample].construct(samples[sample]);     \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_ROMOL_DTOR(DATASET, SUFFIX, TAG)                                 \
  TEST_CASE("ROMol destructor " SUFFIX, "[mol]" TAG) {                         \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK_ADVANCED("ROMol destructor " SUFFIX)(                            \
        Catch::Benchmark::Chronometer meter) {                                 \
      std::vector<Catch::Benchmark::destructable_object<ROMol>> storage(       \
          meter.runs() * samples.size());                                      \
      for (size_t i = 0; i < storage.size(); ++i) {                            \
        storage[i].construct(samples[i % samples.size()]);                     \
      }                                                                        \
      meter.measure([&](int i) {                                               \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          storage[i * samples.size() + sample].destruct();                     \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_RDMOL_COPY(DATASET, SUFFIX, TAG)                                 \
  TEST_CASE("RDMol copy constructor " SUFFIX, "[mol][rdmol]" TAG) {            \
    auto samples = bench_common::load_rdmol_samples(DATASET);                  \
    BENCHMARK_ADVANCED("RDMol copy constructor " SUFFIX)(                      \
        Catch::Benchmark::Chronometer meter) {                                 \
      std::vector<Catch::Benchmark::storage_for<RDMol>> storage(               \
          meter.runs() * samples.size());                                      \
      meter.measure([&](int i) {                                               \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          storage[i * samples.size() + sample].construct(samples[sample]);     \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_RDMOL_DTOR(DATASET, SUFFIX, TAG)                                 \
  TEST_CASE("RDMol destructor " SUFFIX, "[mol][rdmol]" TAG) {                  \
    auto samples = bench_common::load_rdmol_samples(DATASET);                  \
    BENCHMARK_ADVANCED("RDMol destructor " SUFFIX)(                            \
        Catch::Benchmark::Chronometer meter) {                                 \
      std::vector<Catch::Benchmark::destructable_object<RDMol>> storage(       \
          meter.runs() * samples.size());                                      \
      for (size_t i = 0; i < storage.size(); ++i) {                            \
        storage[i].construct(samples[i % samples.size()]);                     \
      }                                                                        \
      meter.measure([&](int i) {                                               \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          storage[i * samples.size() + sample].destruct();                     \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_ROMOL_HEAVYATOMS(DATASET, SUFFIX, TAG)                           \
  TEST_CASE("ROMol::getNumHeavyAtoms " SUFFIX, "[mol]" TAG) {                  \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("ROMol::getNumHeavyAtoms " SUFFIX) {                             \
      auto sum = 0;                                                            \
      for (auto &mol : samples) {                                              \
        sum += mol.getNumHeavyAtoms();                                         \
      }                                                                        \
      return sum;                                                              \
    };                                                                         \
  }

#define BENCH_RDMOL_HEAVYATOMS(DATASET, SUFFIX, TAG)                           \
  TEST_CASE("RDMol::getNumHeavyAtoms " SUFFIX, "[mol][rdmol]" TAG) {           \
    auto samples = bench_common::load_rdmol_samples(DATASET);                  \
    BENCHMARK("RDMol::getNumHeavyAtoms " SUFFIX) {                             \
      auto sum = 0;                                                            \
      for (auto &mol : samples) {                                              \
        sum += mol.getNumHeavyAtoms();                                         \
      }                                                                        \
      return sum;                                                              \
    };                                                                         \
  }

#define BENCH_MOL_FOR(DATASET, SUFFIX, TAG)                                    \
  BENCH_ROMOL_COPY(DATASET, SUFFIX, TAG)                                       \
  BENCH_ROMOL_DTOR(DATASET, SUFFIX, TAG)                                       \
  BENCH_RDMOL_COPY(DATASET, SUFFIX, TAG)                                       \
  BENCH_RDMOL_DTOR(DATASET, SUFFIX, TAG)                                       \
  BENCH_ROMOL_HEAVYATOMS(DATASET, SUFFIX, TAG)                                 \
  BENCH_RDMOL_HEAVYATOMS(DATASET, SUFFIX, TAG)

BENCH_MOL_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_MOL_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_MOL_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_MOL_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_MOL_FOR(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_MOL_FOR(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_MOL_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_MOL_FOR(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_MOL_FOR(Dataset::Rings_6, "rings 6", "[rings_6]")
