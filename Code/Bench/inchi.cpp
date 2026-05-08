#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <INCHI-API/inchi.h>

using namespace RDKit;

using bench_common::Dataset;

TEST_CASE("MolToInchi", "[inchi]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolToInchi") {
    std::vector<std::string> inchis;
    for (auto &mol : samples) {
      ExtraInchiReturnValues rv;
      inchis.push_back(MolToInchi(mol, rv));
    }
    return inchis;
  };
}

TEST_CASE("InchiToInchiKey", "[inchi]") {
  auto samples = bench_common::load_samples();
  std::vector<std::string> inchis;
  for (auto &mol : samples) {
    ExtraInchiReturnValues rv;
    inchis.push_back(MolToInchi(mol, rv));
  }
  BENCHMARK("InchiToInchiKey") {
    std::vector<std::string> inchikeys;
    for (auto &inchi : inchis) {
      inchikeys.push_back(InchiToInchiKey(inchi));
    }
    return inchikeys;
  };
}

TEST_CASE("InchiToMol", "[inchi]") {
  auto samples = bench_common::load_samples();
  std::vector<std::string> inchis;
  for (auto &mol : samples) {
    ExtraInchiReturnValues rv;
    inchis.push_back(MolToInchi(mol, rv));
  }
  BENCHMARK("InchiToMol") {
    std::vector<std::unique_ptr<ROMol>> mols;
    for (auto &inchi : inchis) {
      ExtraInchiReturnValues rv_inner;
      mols.emplace_back(InchiToMol(inchi, rv_inner));
    }
    return mols;
  };
}

#define BENCH_MOL_TO_INCHI(DATASET, SUFFIX, TAG)                               \
  TEST_CASE("MolToInchi " SUFFIX, "[inchi]" TAG) {                             \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("MolToInchi " SUFFIX) {                                          \
      std::vector<std::string> inchis;                                         \
      for (auto &mol : samples) {                                              \
        ExtraInchiReturnValues rv;                                             \
        inchis.push_back(MolToInchi(mol, rv));                                 \
      }                                                                        \
      return inchis;                                                           \
    };                                                                         \
  }

#define BENCH_INCHIKEY(DATASET, SUFFIX, TAG)                                   \
  TEST_CASE("InchiToInchiKey " SUFFIX, "[inchi]" TAG) {                        \
    auto samples = bench_common::load_samples(DATASET);                        \
    std::vector<std::string> inchis;                                           \
    for (auto &mol : samples) {                                                \
      ExtraInchiReturnValues rv;                                               \
      inchis.push_back(MolToInchi(mol, rv));                                   \
    }                                                                          \
    BENCHMARK("InchiToInchiKey " SUFFIX) {                                     \
      std::vector<std::string> inchikeys;                                      \
      for (auto &inchi : inchis) {                                             \
        inchikeys.push_back(InchiToInchiKey(inchi));                           \
      }                                                                        \
      return inchikeys;                                                        \
    };                                                                         \
  }

#define BENCH_INCHI_TO_MOL(DATASET, SUFFIX, TAG)                               \
  TEST_CASE("InchiToMol " SUFFIX, "[inchi]" TAG) {                             \
    auto samples = bench_common::load_samples(DATASET);                        \
    std::vector<std::string> inchis;                                           \
    for (auto &mol : samples) {                                                \
      ExtraInchiReturnValues rv;                                               \
      inchis.push_back(MolToInchi(mol, rv));                                   \
    }                                                                          \
    BENCHMARK("InchiToMol " SUFFIX) {                                          \
      std::vector<std::unique_ptr<ROMol>> mols;                                \
      for (auto &inchi : inchis) {                                             \
        ExtraInchiReturnValues rv_inner;                                       \
        mols.emplace_back(InchiToMol(inchi, rv_inner));                        \
      }                                                                        \
      return mols;                                                             \
    };                                                                         \
  }

#define BENCH_INCHI_FOR(DATASET, SUFFIX, TAG)                                  \
  BENCH_MOL_TO_INCHI(DATASET, SUFFIX, TAG)                                     \
  BENCH_INCHIKEY(DATASET, SUFFIX, TAG)                                         \
  BENCH_INCHI_TO_MOL(DATASET, SUFFIX, TAG)

BENCH_INCHI_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_INCHI_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_INCHI_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_INCHI_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_INCHI_FOR(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_INCHI_FOR(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_INCHI_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_INCHI_FOR(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_INCHI_FOR(Dataset::Rings_6, "rings 6", "[rings_6]")
