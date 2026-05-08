#include <catch2/catch_all.hpp>
#include <sstream>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

using bench_common::Dataset;

TEST_CASE("MolPickler::pickleMol", "[pickle]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolPickler::pickleMol") {
    std::stringstream buf;
    for (auto &mol : samples) {
      MolPickler::pickleMol(mol, buf);
    }
    return buf.str().size();
  };
}

TEST_CASE("MolPickler::molFromPickle", "[pickle]") {
  auto samples = bench_common::load_samples();
  std::vector<std::string> pickles;
  pickles.reserve(samples.size());
  for (auto &mol : samples) {
    std::string pickled;
    MolPickler::pickleMol(mol, pickled);
    pickles.push_back(std::move(pickled));
  }
  BENCHMARK("MolPickler::molFromPickle") {
    auto total_atoms = 0;
    for (auto &pickled : pickles) {
      ROMol res(pickled);
      total_atoms += res.getNumAtoms();
    }
    REQUIRE(total_atoms > 0);
    return total_atoms;
  };
}

#define BENCH_PICKLE(DATASET, SUFFIX, TAG)                                     \
  TEST_CASE("MolPickler::pickleMol " SUFFIX, "[pickle]" TAG) {                 \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolPickler::pickleMol " SUFFIX) {                               \
      std::stringstream buf;                                                   \
      for (auto &mol : samples) {                                              \
        MolPickler::pickleMol(mol, buf);                                       \
      }                                                                        \
      return buf.str().size();                                                 \
    };                                                                         \
  }

#define BENCH_UNPICKLE(DATASET, SUFFIX, TAG)                                   \
  TEST_CASE("MolPickler::molFromPickle " SUFFIX, "[pickle]" TAG) {             \
    auto samples = bench_common::load_samples(DATASET);                        \
    std::vector<std::string> pickles;                                          \
    pickles.reserve(samples.size());                                           \
    for (auto &mol : samples) {                                                \
      std::string pickled;                                                     \
      MolPickler::pickleMol(mol, pickled);                                     \
      pickles.push_back(std::move(pickled));                                   \
    }                                                                          \
    BENCHMARK("MolPickler::molFromPickle " SUFFIX) {                           \
      auto total_atoms = 0;                                                    \
      for (auto &pickled : pickles) {                                          \
        ROMol res(pickled);                                                    \
        total_atoms += res.getNumAtoms();                                      \
      }                                                                        \
      REQUIRE(total_atoms > 0);                                                \
      return total_atoms;                                                      \
    };                                                                         \
  }

#define BENCH_PICKLE_FOR(DATASET, SUFFIX, TAG)                                 \
  BENCH_PICKLE(DATASET, SUFFIX, TAG)                                           \
  BENCH_UNPICKLE(DATASET, SUFFIX, TAG)

BENCH_PICKLE_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_PICKLE_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_PICKLE_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_PICKLE_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_PICKLE_FOR(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_PICKLE_FOR(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_PICKLE_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_PICKLE_FOR(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_PICKLE_FOR(Dataset::Rings_6, "rings 6", "[rings_6]")
