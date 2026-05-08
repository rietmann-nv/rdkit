#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/Descriptors/Lipinski.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

using bench_common::Dataset;

TEST_CASE("Descriptors::calcNumSpiroAtoms", "[descriptors]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("Descriptors::calcNumSpiroAtoms") {
    auto sum = 0;
    for (auto &mol : samples) {
      sum += Descriptors::calcNumSpiroAtoms(mol);
    }
    return sum;
  };
}

TEST_CASE("Descriptors::calcNumBridgeheadAtoms", "[descriptors]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("Descriptors::calcNumBridgeheadAtoms") {
    auto sum = 0;
    for (auto &mol : samples) {
      sum += Descriptors::calcNumBridgeheadAtoms(mol);
    }
    return sum;
  };
}

#define BENCH_SPIRO(DATASET, SUFFIX, TAG)                                      \
  TEST_CASE("Descriptors::calcNumSpiroAtoms " SUFFIX, "[descriptors]" TAG) {   \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("Descriptors::calcNumSpiroAtoms " SUFFIX) {                      \
      auto sum = 0;                                                            \
      for (auto &mol : samples) {                                              \
        sum += Descriptors::calcNumSpiroAtoms(mol);                            \
      }                                                                        \
      return sum;                                                              \
    };                                                                         \
  }

#define BENCH_BRIDGEHEAD(DATASET, SUFFIX, TAG)                                 \
  TEST_CASE("Descriptors::calcNumBridgeheadAtoms " SUFFIX,                     \
            "[descriptors]" TAG) {                                             \
    auto samples = bench_common::load_samples(DATASET);                        \
    BENCHMARK("Descriptors::calcNumBridgeheadAtoms " SUFFIX) {                 \
      auto sum = 0;                                                            \
      for (auto &mol : samples) {                                              \
        sum += Descriptors::calcNumBridgeheadAtoms(mol);                       \
      }                                                                        \
      return sum;                                                              \
    };                                                                         \
  }

#define BENCH_DESC_FOR(DATASET, SUFFIX, TAG)                                   \
  BENCH_SPIRO(DATASET, SUFFIX, TAG)                                            \
  BENCH_BRIDGEHEAD(DATASET, SUFFIX, TAG)

BENCH_DESC_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_DESC_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_DESC_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_DESC_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_DESC_FOR(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_DESC_FOR(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_DESC_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_DESC_FOR(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_DESC_FOR(Dataset::Rings_6, "rings 6", "[rings_6]")
