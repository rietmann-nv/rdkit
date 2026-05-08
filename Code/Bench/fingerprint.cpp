#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/Fingerprints/MorganGenerator.h>
#include <GraphMol/RDMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

using bench_common::Dataset;

TEST_CASE("MorganFingerprints::getFingerprint", "[fingerprint]") {
  auto samples = bench_common::load_samples();
  const auto radius = 2;
  std::unique_ptr<FingerprintGenerator<uint64_t>> gen(
      MorganFingerprint::getMorganGenerator<uint64_t>(radius));

  BENCHMARK("MorganFingerprints::getFingerprint") {
    auto sum = 0;
    for (auto &mol : samples) {
      std::unique_ptr<ExplicitBitVect> fp(gen->getFingerprint(mol));
      sum += fp->getNumOnBits();
    }
    return sum;
  };
}

TEST_CASE("MorganFingerprints::getFingerprint RDMol",
          "[fingerprint][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();
  const auto radius = 2;
  std::unique_ptr<FingerprintGenerator<uint64_t>> gen(
      MorganFingerprint::getMorganGenerator<uint64_t>(radius));

  BENCHMARK("MorganFingerprints::getFingerprint RDMol") {
    auto sum = 0;
    for (const auto &mol : samples) {
      std::unique_ptr<ExplicitBitVect> fp(gen->getFingerprint(mol));
      sum += fp->getNumOnBits();
    }
    return sum;
  };
}

#define BENCH_MORGAN(DATASET, SUFFIX, TAG)                                     \
  TEST_CASE("MorganFingerprints::getFingerprint " SUFFIX,                      \
            "[fingerprint]" TAG) {                                             \
    auto samples = bench_common::load_samples(DATASET);                        \
    const auto radius = 2;                                                     \
    std::unique_ptr<FingerprintGenerator<uint64_t>> gen(                       \
        MorganFingerprint::getMorganGenerator<uint64_t>(radius));              \
    BENCHMARK("MorganFingerprints::getFingerprint " SUFFIX) {                  \
      auto sum = 0;                                                            \
      for (auto &mol : samples) {                                              \
        std::unique_ptr<ExplicitBitVect> fp(gen->getFingerprint(mol));         \
        sum += fp->getNumOnBits();                                             \
      }                                                                        \
      return sum;                                                              \
    };                                                                         \
  }                                                                            \
  TEST_CASE("MorganFingerprints::getFingerprint RDMol " SUFFIX,                \
            "[fingerprint][rdmol]" TAG) {                                      \
    auto samples = bench_common::load_rdmol_samples(DATASET);                  \
    const auto radius = 2;                                                     \
    std::unique_ptr<FingerprintGenerator<uint64_t>> gen(                       \
        MorganFingerprint::getMorganGenerator<uint64_t>(radius));              \
    BENCHMARK("MorganFingerprints::getFingerprint RDMol " SUFFIX) {            \
      auto sum = 0;                                                            \
      for (const auto &mol : samples) {                                        \
        std::unique_ptr<ExplicitBitVect> fp(gen->getFingerprint(mol));         \
        sum += fp->getNumOnBits();                                             \
      }                                                                        \
      return sum;                                                              \
    };                                                                         \
  }

BENCH_MORGAN(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_MORGAN(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_MORGAN(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_MORGAN(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_MORGAN(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_MORGAN(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_MORGAN(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_MORGAN(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_MORGAN(Dataset::Rings_6, "rings 6", "[rings_6]")
