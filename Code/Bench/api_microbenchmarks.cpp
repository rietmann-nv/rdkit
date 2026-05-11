// Micro-benchmarks for the ROMol vs RDMol public API surface itself,
// independent of the higher-level MolOps / SmilesToMol / Fingerprints
// pipelines.  Each test case has a ROMol leg ("[mol_api]") and an RDMol
// leg ("[mol_api][rdmol]") that exercises the same logical operation.
//
// Property benches come in pairs.
//   *_bulk           - allowed to use the RDMol bulk path (single
//                      addAtomProp call, raw pointer write, etc.) on the
//                      RDMol leg; ROMol leg uses its natural per-atom
//                      setProp loop.  Measures the win when the API call
//                      shape is allowed to differ.
//   *_random_access  - both legs do per-element get/set on a pre-existing
//                      property in a deterministic non-sequential index
//                      order.  No allocation, no bulk shortcut.  Measures
//                      the steady-state mutate/read-one-prop cost.
//
// Attribute extraction comes in two flavours where applicable.
//   *_loop  - both sides walk atoms/bonds and copy a single attribute
//             into a destination vector.  Apples-to-apples loop pattern.
//   *_bulk  - RDMol-only; uses getAtomPropArrayIfPresent or equivalent
//             to extract the whole array via memcpy/std::copy.  No
//             ROMol leg because ROMol does not expose contiguous prop
//             storage.
//
// Iteration benches keep the per-element op minimal (sum of one
// integer-typed attribute) so the measurement reflects access pattern,
// not work.

#include <catch2/catch_all.hpp>

#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

#include <GraphMol/MolOps.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RDMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RingInfo.h>
#include <RDGeneral/types.h>

#include "bench_common.hpp"

using namespace RDKit;

namespace {

using bench_common::Dataset;
using bench_common::load_rdmol_samples;
using bench_common::load_samples;
using bench_common::nth_random;

// Run OP once per sample, with a deep copy of each input made up-front
// (one copy per sample per Catch2 run) so OPs that mutate the molecule
// see a fresh state on every iteration.  Returns the accumulator from
// the last meter call, kept live to defeat dead-store elimination.
template <class Mol, class Op>
auto run_per_sample_mutating(const std::vector<Mol> &samples,
                             Catch::Benchmark::Chronometer meter, Op op) {
  std::vector<Mol> work;
  work.reserve(meter.runs() * samples.size());
  for (int run = 0; run < meter.runs(); ++run) {
    for (const auto &mol : samples) {
      work.emplace_back(mol);
    }
  }
  meter.measure([&](int run) {
    uint64_t total = 0;
    for (size_t sample = 0; sample < samples.size(); ++sample) {
      auto &mol = work[run * samples.size() + sample];
      total += op(mol);
    }
    return total;
  });
}

// Read-only variant: each meter iteration walks the same const samples,
// no copy needed.  Used for getters / extractors / iteration benches.
template <class Mol, class Op>
auto run_per_sample_readonly(const std::vector<Mol> &samples,
                             Catch::Benchmark::Chronometer meter, Op op) {
  meter.measure([&](int /*run*/) {
    uint64_t total = 0;
    for (const auto &mol : samples) {
      total += op(mol);
    }
    return total;
  });
}

// Build a primed copy of the samples with `prime` applied once per mol.
// Used for benches that need ring info, an existing property, etc.
template <class Mol, class Prime>
std::vector<Mol> prime_samples(std::vector<Mol> samples, Prime prime) {
  for (auto &mol : samples) {
    prime(mol);
  }
  return samples;
}

// Deterministic non-sequential index sequence over [0, n).  Same sequence
// for ROMol and RDMol legs so they walk the same access pattern.
std::vector<uint32_t> shuffled_indices(uint32_t n) {
  std::vector<uint32_t> indices(n);
  std::iota(indices.begin(), indices.end(), 0u);
  for (uint32_t i = 0; i + 1 < n; ++i) {
    uint32_t j = i + uint32_t(nth_random(i) % (n - i));
    std::swap(indices[i], indices[j]);
  }
  return indices;
}

// Pre-compute a shuffled index sequence per sample so the random-access
// benches don't pay for the shuffle inside the timed region.
template <class Mol>
std::vector<std::vector<uint32_t>> per_sample_shuffled_atom_indices(
    const std::vector<Mol> &samples) {
  std::vector<std::vector<uint32_t>> result;
  result.reserve(samples.size());
  for (const auto &mol : samples) {
    result.emplace_back(shuffled_indices(mol.getNumAtoms()));
  }
  return result;
}

}  // namespace

// ---------------------------------------------------------------------------
// Group 1: lifecycle (move-only; copy/dtor live in mol.cpp)
// ---------------------------------------------------------------------------

#define BENCH_MOVE_CTOR(DATASET, SUFFIX, TAG)                                  \
  TEST_CASE("ROMol move constructor " SUFFIX, "[mol_api]" TAG) {               \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol move constructor " SUFFIX)                       \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<ROMol> sources;                                              \
      sources.reserve(meter.runs() * samples.size());                          \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (const auto &mol : samples) {                                      \
          sources.emplace_back(mol);                                           \
        }                                                                      \
      }                                                                        \
      std::vector<Catch::Benchmark::storage_for<ROMol>> storage(               \
          meter.runs() * samples.size());                                      \
      meter.measure([&](int run) {                                             \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          storage[run * samples.size() + sample].construct(                    \
              std::move(sources[run * samples.size() + sample]));              \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol move constructor " SUFFIX, "[mol_api][rdmol]" TAG) {        \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol move constructor " SUFFIX)                       \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<RDMol> sources;                                              \
      sources.reserve(meter.runs() * samples.size());                          \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (const auto &mol : samples) {                                      \
          sources.emplace_back(mol);                                           \
        }                                                                      \
      }                                                                        \
      std::vector<Catch::Benchmark::storage_for<RDMol>> storage(               \
          meter.runs() * samples.size());                                      \
      meter.measure([&](int run) {                                             \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          storage[run * samples.size() + sample].construct(                    \
              std::move(sources[run * samples.size() + sample]));              \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_MOVE_ASSIGN(DATASET, SUFFIX, TAG)                                \
  TEST_CASE("ROMol move assign " SUFFIX, "[mol_api]" TAG) {                    \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol move assign " SUFFIX)                            \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<ROMol> sources;                                              \
      std::vector<ROMol> dests;                                                \
      sources.reserve(meter.runs() * samples.size());                          \
      dests.reserve(meter.runs() * samples.size());                            \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (const auto &mol : samples) {                                      \
          sources.emplace_back(mol);                                           \
          dests.emplace_back();                                                \
        }                                                                      \
      }                                                                        \
      meter.measure([&](int run) {                                             \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          dests[run * samples.size() + sample] =                               \
              std::move(sources[run * samples.size() + sample]);               \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol move assign " SUFFIX, "[mol_api][rdmol]" TAG) {             \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol move assign " SUFFIX)                            \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<RDMol> sources;                                              \
      std::vector<RDMol> dests;                                                \
      sources.reserve(meter.runs() * samples.size());                          \
      dests.reserve(meter.runs() * samples.size());                            \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (const auto &mol : samples) {                                      \
          sources.emplace_back(mol);                                           \
          dests.emplace_back();                                                \
        }                                                                      \
      }                                                                        \
      meter.measure([&](int run) {                                             \
        for (size_t sample = 0; sample < samples.size(); ++sample) {           \
          dests[run * samples.size() + sample] =                               \
              std::move(sources[run * samples.size() + sample]);               \
        }                                                                      \
      });                                                                      \
    };                                                                         \
  }

// ---------------------------------------------------------------------------
// Group 2: atom & bond iteration with minimal op
// ---------------------------------------------------------------------------

#define BENCH_ATOM_ITER_SUM_ATOMICNUM(DATASET, SUFFIX, TAG)                    \
  TEST_CASE("ROMol atoms() sum atomicNum " SUFFIX, "[mol_api]" TAG) {          \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol atoms() sum atomicNum " SUFFIX)                  \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->getAtomicNum();                                       \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol atoms loop sum atomicNum " SUFFIX,                          \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol atoms loop sum atomicNum " SUFFIX)               \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          total += mol.getAtom(i).getAtomicNum();                              \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_BOND_ITER_SUM_TWICEORDER(DATASET, SUFFIX, TAG)                   \
  TEST_CASE("ROMol bonds() sum twiceBondType " SUFFIX, "[mol_api]" TAG) {      \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol bonds() sum twiceBondType " SUFFIX)              \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto bond : mol.bonds()) {                                        \
          total += getTwiceBondType(bond->getBondType());                      \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol bonds loop sum twiceBondType " SUFFIX,                      \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol bonds loop sum twiceBondType " SUFFIX)           \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const uint32_t numBonds = mol.getNumBonds();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numBonds; ++i) {                              \
          total += getTwiceBondType(mol.getBond(i).getBondType());             \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_NEIGHBOR_WALK_SUM_DEGREE(DATASET, SUFFIX, TAG)                   \
  TEST_CASE("ROMol atomNeighbors sum degree " SUFFIX, "[mol_api]" TAG) {       \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol atomNeighbors sum degree " SUFFIX)               \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          for (auto nbr : mol.atomNeighbors(atom)) {                           \
            (void)nbr;                                                         \
            ++total;                                                           \
          }                                                                    \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol atomNeighbors range sum degree " SUFFIX,                    \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol atomNeighbors range sum degree " SUFFIX)         \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          for (auto nbr : mol.atomNeighbors(i)) {                              \
            (void)nbr;                                                         \
            ++total;                                                           \
          }                                                                    \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol getAtomNeighbors raw sum degree " SUFFIX,                   \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol getAtomNeighbors raw sum degree " SUFFIX)        \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          auto [begin, end] = mol.getAtomNeighbors(i);                         \
          total += uint32_t(end - begin);                                      \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_BONDS_FROM_ATOM_SUM_ORDER(DATASET, SUFFIX, TAG)                  \
  TEST_CASE("ROMol atomBonds sum twiceBondType " SUFFIX, "[mol_api]" TAG) {    \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol atomBonds sum twiceBondType " SUFFIX)            \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          for (auto bond : mol.atomBonds(atom)) {                              \
            total += getTwiceBondType(bond->getBondType());                    \
          }                                                                    \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol atomBonds range sum twiceBondType " SUFFIX,                 \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol atomBonds range sum twiceBondType " SUFFIX)      \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          for (auto &&bond : mol.atomBonds(i)) {                               \
            total += getTwiceBondType(bond.data().getBondType());              \
          }                                                                    \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol getAtomBonds raw sum twiceBondType " SUFFIX,                \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol getAtomBonds raw sum twiceBondType " SUFFIX)     \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          auto [begin, end] = mol.getAtomBonds(i);                             \
          for (const uint32_t *p = begin; p != end; ++p) {                     \
            total += getTwiceBondType(mol.getBond(*p).getBondType());          \
          }                                                                    \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

// ---------------------------------------------------------------------------
// Group 3: intrinsic-attribute extraction (push_back into a vector)
// ---------------------------------------------------------------------------

#define BENCH_EXTRACT_ATOMICNUMS(DATASET, SUFFIX, TAG)                         \
  TEST_CASE("ROMol extract atomicNums " SUFFIX, "[mol_api]" TAG) {             \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol extract atomicNums " SUFFIX)                     \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<int> sink;                                                   \
      run_per_sample_readonly(samples, meter, [&sink](const ROMol &mol) {     \
        sink.clear();                                                          \
        sink.reserve(mol.getNumAtoms());                                       \
        for (auto atom : mol.atoms()) {                                        \
          sink.push_back(atom->getAtomicNum());                                \
        }                                                                      \
        return uint64_t(sink.size());                                          \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol extract atomicNums " SUFFIX, "[mol_api][rdmol]" TAG) {      \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol extract atomicNums " SUFFIX)                     \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<uint32_t> sink;                                              \
      run_per_sample_readonly(samples, meter, [&sink](const RDMol &mol) {     \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        sink.clear();                                                          \
        sink.reserve(numAtoms);                                                \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          sink.push_back(mol.getAtom(i).getAtomicNum());                       \
        }                                                                      \
        return uint64_t(sink.size());                                          \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_EXTRACT_FORMAL_CHARGES(DATASET, SUFFIX, TAG)                     \
  TEST_CASE("ROMol extract formalCharges " SUFFIX, "[mol_api]" TAG) {          \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol extract formalCharges " SUFFIX)                  \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<int> sink;                                                   \
      run_per_sample_readonly(samples, meter, [&sink](const ROMol &mol) {     \
        sink.clear();                                                          \
        sink.reserve(mol.getNumAtoms());                                       \
        for (auto atom : mol.atoms()) {                                        \
          sink.push_back(atom->getFormalCharge());                             \
        }                                                                      \
        return uint64_t(sink.size());                                          \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol extract formalCharges " SUFFIX, "[mol_api][rdmol]" TAG) {   \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol extract formalCharges " SUFFIX)                  \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<int8_t> sink;                                                \
      run_per_sample_readonly(samples, meter, [&sink](const RDMol &mol) {     \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        sink.clear();                                                          \
        sink.reserve(numAtoms);                                                \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          sink.push_back(mol.getAtom(i).getFormalCharge());                    \
        }                                                                      \
        return uint64_t(sink.size());                                          \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_EXTRACT_MASSES(DATASET, SUFFIX, TAG)                             \
  TEST_CASE("ROMol extract masses " SUFFIX, "[mol_api]" TAG) {                 \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol extract masses " SUFFIX)                         \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<double> sink;                                                \
      run_per_sample_readonly(samples, meter, [&sink](const ROMol &mol) {     \
        sink.clear();                                                          \
        sink.reserve(mol.getNumAtoms());                                       \
        for (auto atom : mol.atoms()) {                                        \
          sink.push_back(atom->getMass());                                     \
        }                                                                      \
        return uint64_t(sink.size());                                          \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol extract masses " SUFFIX, "[mol_api][rdmol]" TAG) {          \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol extract masses " SUFFIX)                         \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<double> sink;                                                \
      run_per_sample_readonly(samples, meter, [&sink](const RDMol &mol) {     \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        sink.clear();                                                          \
        sink.reserve(numAtoms);                                                \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          sink.push_back(mol.getAtom(i).getMass());                            \
        }                                                                      \
        return uint64_t(sink.size());                                          \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_EXTRACT_BOND_TYPES(DATASET, SUFFIX, TAG)                         \
  TEST_CASE("ROMol extract bondTypes " SUFFIX, "[mol_api]" TAG) {              \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol extract bondTypes " SUFFIX)                      \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<uint8_t> sink;                                               \
      run_per_sample_readonly(samples, meter, [&sink](const ROMol &mol) {     \
        sink.clear();                                                          \
        sink.reserve(mol.getNumBonds());                                       \
        for (auto bond : mol.bonds()) {                                        \
          sink.push_back(uint8_t(bond->getBondType()));                        \
        }                                                                      \
        return uint64_t(sink.size());                                          \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol extract bondTypes " SUFFIX, "[mol_api][rdmol]" TAG) {       \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol extract bondTypes " SUFFIX)                      \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<uint8_t> sink;                                               \
      run_per_sample_readonly(samples, meter, [&sink](const RDMol &mol) {     \
        const uint32_t numBonds = mol.getNumBonds();                           \
        sink.clear();                                                          \
        sink.reserve(numBonds);                                                \
        for (uint32_t i = 0; i < numBonds; ++i) {                              \
          sink.push_back(uint8_t(mol.getBond(i).getBondType()));               \
        }                                                                      \
        return uint64_t(sink.size());                                          \
      });                                                                      \
    };                                                                         \
  }

// ---------------------------------------------------------------------------
// Group 4: user property -- bulk and random-access
//
// The ROMol leg stores per-atom int props via Atom::setProp(key, value).
// The RDMol leg has three flavours:
//   _loop : per-atom setSingleAtomProp(token, idx, value) -- matches the
//           ROMol access pattern but uses RDMol's PropToken-keyed array.
//   _bulk : addAtomProp(token, default) returns a T* into the contiguous
//           array, then a tight loop fills it.  This is the RDMol-native
//           bulk path.
//   _broadcast : addAtomProp(token, value) -- single call sets every
//           atom's value in one shot.  Only meaningful for uniform-value
//           writes.
// ---------------------------------------------------------------------------

#define BENCH_PROP_SET_VARIED(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("ROMol setProp int per-atom varied " SUFFIX, "[mol_api]" TAG) {    \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol setProp int per-atom varied " SUFFIX)            \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_mutating(samples, meter, [](ROMol &mol) {                 \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          mol.getAtomWithIdx(i)->setProp("bench_int", int(i));                 \
        }                                                                      \
        return uint64_t(numAtoms);                                             \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol setSingleAtomProp int per-atom varied " SUFFIX,             \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    PropToken token("bench_int");                                              \
    BENCHMARK_ADVANCED("RDMol setSingleAtomProp int per-atom varied " SUFFIX)  \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_mutating(samples, meter, [&token](RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          mol.setSingleAtomProp(token, i, int(i));                             \
        }                                                                      \
        return uint64_t(numAtoms);                                             \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol addAtomProp int bulk-pointer varied " SUFFIX,               \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    PropToken token("bench_int");                                              \
    BENCHMARK_ADVANCED("RDMol addAtomProp int bulk-pointer varied " SUFFIX)    \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_mutating(samples, meter, [&token](RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        int *data = mol.addAtomProp<int>(token, 0);                            \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          data[i] = int(i);                                                    \
        }                                                                      \
        return uint64_t(numAtoms);                                             \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_PROP_SET_UNIFORM(DATASET, SUFFIX, TAG)                           \
  TEST_CASE("ROMol setProp int per-atom uniform " SUFFIX, "[mol_api]" TAG) {   \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol setProp int per-atom uniform " SUFFIX)           \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_mutating(samples, meter, [](ROMol &mol) {                 \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          mol.getAtomWithIdx(i)->setProp("bench_uniform", int(42));            \
        }                                                                      \
        return uint64_t(numAtoms);                                             \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol addAtomProp int broadcast uniform " SUFFIX,                 \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    PropToken token("bench_uniform");                                          \
    BENCHMARK_ADVANCED("RDMol addAtomProp int broadcast uniform " SUFFIX)      \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_mutating(samples, meter, [&token](RDMol &mol) {           \
        mol.addAtomProp<int>(token, 42);                                       \
        return uint64_t(mol.getNumAtoms());                                    \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol setSingleAtomProp int per-atom uniform " SUFFIX,            \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    PropToken token("bench_uniform");                                          \
    BENCHMARK_ADVANCED("RDMol setSingleAtomProp int per-atom uniform "         \
                       SUFFIX)                                                 \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_mutating(samples, meter, [&token](RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          mol.setSingleAtomProp(token, i, int(42));                            \
        }                                                                      \
        return uint64_t(numAtoms);                                             \
      });                                                                      \
    };                                                                         \
  }

// Random-access mutate of a property that already exists.  Both legs
// pre-populate before the timed region so neither pays for allocation.
#define BENCH_PROP_SET_RANDOM_ACCESS(DATASET, SUFFIX, TAG)                     \
  TEST_CASE("ROMol setProp int random-access " SUFFIX, "[mol_api]" TAG) {      \
    auto raw = load_samples(DATASET);                                          \
    auto samples = prime_samples(std::move(raw), [](ROMol &mol) {              \
      const uint32_t numAtoms = mol.getNumAtoms();                             \
      for (uint32_t i = 0; i < numAtoms; ++i) {                                \
        mol.getAtomWithIdx(i)->setProp("bench_int", int(0));                   \
      }                                                                        \
    });                                                                        \
    auto orderings = per_sample_shuffled_atom_indices(samples);                \
    BENCHMARK_ADVANCED("ROMol setProp int random-access " SUFFIX)              \
    (Catch::Benchmark::Chronometer meter) {                                    \
      meter.measure([&](int /*run*/) {                                         \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          auto &mol = const_cast<ROMol &>(samples[s]);                         \
          for (uint32_t i : orderings[s]) {                                    \
            mol.getAtomWithIdx(i)->setProp("bench_int", int(i + 1));           \
            ++total;                                                           \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol setSingleAtomProp int random-access " SUFFIX,               \
            "[mol_api][rdmol]" TAG) {                                          \
    PropToken token("bench_int");                                              \
    auto raw = load_rdmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [&token](RDMol &mol) {        \
      mol.addAtomProp<int>(token, 0);                                          \
    });                                                                        \
    auto orderings = per_sample_shuffled_atom_indices(samples);                \
    BENCHMARK_ADVANCED("RDMol setSingleAtomProp int random-access " SUFFIX)    \
    (Catch::Benchmark::Chronometer meter) {                                    \
      meter.measure([&](int /*run*/) {                                         \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          auto &mol = const_cast<RDMol &>(samples[s]);                         \
          for (uint32_t i : orderings[s]) {                                    \
            mol.setSingleAtomProp(token, i, int(i + 1));                       \
            ++total;                                                           \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol getAtomPropArrayIfPresent int random-access " SUFFIX,       \
            "[mol_api][rdmol]" TAG) {                                          \
    PropToken token("bench_int");                                              \
    auto raw = load_rdmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [&token](RDMol &mol) {        \
      mol.addAtomProp<int>(token, 0);                                          \
    });                                                                        \
    auto orderings = per_sample_shuffled_atom_indices(samples);                \
    BENCHMARK_ADVANCED(                                                        \
        "RDMol getAtomPropArrayIfPresent int random-access " SUFFIX)           \
    (Catch::Benchmark::Chronometer meter) {                                    \
      meter.measure([&](int /*run*/) {                                         \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          auto &mol = const_cast<RDMol &>(samples[s]);                         \
          int *data = mol.getAtomPropArrayIfPresent<int>(token);               \
          for (uint32_t i : orderings[s]) {                                    \
            data[i] = int(i + 1);                                              \
            ++total;                                                           \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }

// Random-access read of a property that already exists.
#define BENCH_PROP_GET_RANDOM_ACCESS(DATASET, SUFFIX, TAG)                     \
  TEST_CASE("ROMol getProp int random-access " SUFFIX, "[mol_api]" TAG) {      \
    auto raw = load_samples(DATASET);                                          \
    auto samples = prime_samples(std::move(raw), [](ROMol &mol) {              \
      const uint32_t numAtoms = mol.getNumAtoms();                             \
      for (uint32_t i = 0; i < numAtoms; ++i) {                                \
        mol.getAtomWithIdx(i)->setProp("bench_int", int(i));                   \
      }                                                                        \
    });                                                                        \
    auto orderings = per_sample_shuffled_atom_indices(samples);                \
    BENCHMARK_ADVANCED("ROMol getProp int random-access " SUFFIX)              \
    (Catch::Benchmark::Chronometer meter) {                                    \
      meter.measure([&](int /*run*/) {                                         \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          const auto &mol = samples[s];                                        \
          for (uint32_t i : orderings[s]) {                                    \
            total += uint32_t(mol.getAtomWithIdx(i)->getProp<int>(             \
                "bench_int"));                                                 \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol getAtomProp int random-access " SUFFIX,                     \
            "[mol_api][rdmol]" TAG) {                                          \
    PropToken token("bench_int");                                              \
    auto raw = load_rdmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [&token](RDMol &mol) {        \
      const uint32_t numAtoms = mol.getNumAtoms();                             \
      int *data = mol.addAtomProp<int>(token, 0);                              \
      for (uint32_t i = 0; i < numAtoms; ++i) {                                \
        data[i] = int(i);                                                      \
      }                                                                        \
    });                                                                        \
    auto orderings = per_sample_shuffled_atom_indices(samples);                \
    BENCHMARK_ADVANCED("RDMol getAtomProp int random-access " SUFFIX)          \
    (Catch::Benchmark::Chronometer meter) {                                    \
      meter.measure([&](int /*run*/) {                                         \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          const auto &mol = samples[s];                                        \
          for (uint32_t i : orderings[s]) {                                    \
            total += uint32_t(mol.getAtomProp<int>(token, i));                 \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol getAtomPropArrayIfPresent int random-access read "          \
            SUFFIX, "[mol_api][rdmol]" TAG) {                                  \
    PropToken token("bench_int");                                              \
    auto raw = load_rdmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [&token](RDMol &mol) {        \
      const uint32_t numAtoms = mol.getNumAtoms();                             \
      int *data = mol.addAtomProp<int>(token, 0);                              \
      for (uint32_t i = 0; i < numAtoms; ++i) {                                \
        data[i] = int(i);                                                      \
      }                                                                        \
    });                                                                        \
    auto orderings = per_sample_shuffled_atom_indices(samples);                \
    BENCHMARK_ADVANCED(                                                        \
        "RDMol getAtomPropArrayIfPresent int random-access read " SUFFIX)      \
    (Catch::Benchmark::Chronometer meter) {                                    \
      meter.measure([&](int /*run*/) {                                         \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < samples.size(); ++s) {                          \
          const auto &mol = samples[s];                                        \
          const int *data = mol.getAtomPropArrayIfPresent<int>(token);         \
          for (uint32_t i : orderings[s]) {                                    \
            total += uint32_t(data[i]);                                        \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }

// hasProp scan over all atoms for a property that IS present (set by
// the prime step).
#define BENCH_PROP_HAS_PRESENT(DATASET, SUFFIX, TAG)                           \
  TEST_CASE("ROMol hasProp present " SUFFIX, "[mol_api]" TAG) {                \
    auto raw = load_samples(DATASET);                                          \
    auto samples = prime_samples(std::move(raw), [](ROMol &mol) {              \
      for (auto atom : mol.atoms()) {                                          \
        atom->setProp("bench_int", int(0));                                    \
      }                                                                        \
    });                                                                        \
    BENCHMARK_ADVANCED("ROMol hasProp present " SUFFIX)                        \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->hasProp("bench_int") ? 1u : 0u;                       \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol hasAtomProp present " SUFFIX, "[mol_api][rdmol]" TAG) {     \
    PropToken token("bench_int");                                              \
    auto raw = load_rdmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [&token](RDMol &mol) {        \
      mol.addAtomProp<int>(token, 0);                                          \
    });                                                                        \
    BENCHMARK_ADVANCED("RDMol hasAtomProp present " SUFFIX)                    \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [&token](const RDMol &mol) {     \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          total += mol.hasAtomProp(token, i) ? 1u : 0u;                        \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

// hasProp scan for a property that does NOT exist (cold lookup).
#define BENCH_PROP_HAS_ABSENT(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("ROMol hasProp absent " SUFFIX, "[mol_api]" TAG) {                 \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol hasProp absent " SUFFIX)                         \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->hasProp("bench_absent") ? 1u : 0u;                    \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol hasAtomProp absent " SUFFIX, "[mol_api][rdmol]" TAG) {      \
    PropToken token("bench_absent");                                           \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol hasAtomProp absent " SUFFIX)                     \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [&token](const RDMol &mol) {     \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          total += mol.hasAtomProp(token, i) ? 1u : 0u;                        \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

// Bulk clear of an atom-property array.  ROMol clears each atom's prop;
// RDMol drops the contiguous array in one call.
#define BENCH_PROP_CLEAR_BULK(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("ROMol clearProp per-atom " SUFFIX, "[mol_api]" TAG) {             \
    auto raw = load_samples(DATASET);                                          \
    BENCHMARK_ADVANCED("ROMol clearProp per-atom " SUFFIX)                     \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<ROMol> work;                                                 \
      work.reserve(meter.runs() * raw.size());                                 \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (auto &mol : raw) {                                                \
          work.emplace_back(mol);                                              \
          for (auto atom : work.back().atoms()) {                              \
            atom->setProp("bench_int", int(0));                                \
          }                                                                    \
        }                                                                      \
      }                                                                        \
      meter.measure([&](int run) {                                             \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < raw.size(); ++s) {                              \
          auto &mol = work[run * raw.size() + s];                              \
          for (auto atom : mol.atoms()) {                                      \
            atom->clearProp("bench_int");                                      \
            ++total;                                                           \
          }                                                                    \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol clearAtomPropIfPresent " SUFFIX,                            \
            "[mol_api][rdmol]" TAG) {                                          \
    PropToken token("bench_int");                                              \
    auto raw = load_rdmol_samples(DATASET);                                    \
    BENCHMARK_ADVANCED("RDMol clearAtomPropIfPresent " SUFFIX)                 \
    (Catch::Benchmark::Chronometer meter) {                                    \
      std::vector<RDMol> work;                                                 \
      work.reserve(meter.runs() * raw.size());                                 \
      for (int run = 0; run < meter.runs(); ++run) {                           \
        for (auto &mol : raw) {                                                \
          work.emplace_back(mol);                                              \
          work.back().addAtomProp<int>(token, 0);                              \
        }                                                                      \
      }                                                                        \
      meter.measure([&](int run) {                                             \
        uint64_t total = 0;                                                    \
        for (size_t s = 0; s < raw.size(); ++s) {                              \
          auto &mol = work[run * raw.size() + s];                              \
          mol.clearAtomPropIfPresent(token);                                   \
          total += mol.getNumAtoms();                                          \
        }                                                                      \
        return total;                                                          \
      });                                                                      \
    };                                                                         \
  }

// ---------------------------------------------------------------------------
// Group 5: simple cached accessors (degree, totalDegree, numImplicitHs)
// All depend on the property cache being populated, which load_*_samples
// guarantees by passing sanitize=true.
// ---------------------------------------------------------------------------

#define BENCH_ACCESSOR_DEGREE(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("ROMol getDegree all atoms " SUFFIX, "[mol_api]" TAG) {            \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol getDegree all atoms " SUFFIX)                    \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->getDegree();                                          \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol getAtomDegree all atoms " SUFFIX,                           \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol getAtomDegree all atoms " SUFFIX)                \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          total += mol.getAtomDegree(i);                                       \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_ACCESSOR_NUM_IMPLICIT_HS(DATASET, SUFFIX, TAG)                   \
  TEST_CASE("ROMol getNumImplicitHs all atoms " SUFFIX, "[mol_api]" TAG) {     \
    auto samples = load_samples(DATASET);                                      \
    BENCHMARK_ADVANCED("ROMol getNumImplicitHs all atoms " SUFFIX)             \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += atom->getNumImplicitHs();                                   \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol getNumImplicitHs all atoms " SUFFIX,                        \
            "[mol_api][rdmol]" TAG) {                                          \
    auto samples = load_rdmol_samples(DATASET);                                \
    BENCHMARK_ADVANCED("RDMol getNumImplicitHs all atoms " SUFFIX)             \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          total += mol.getAtom(i).getNumImplicitHs();                          \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

// ---------------------------------------------------------------------------
// Group 6: ring info (cacheable -- prime via findSSSR before measure)
// ---------------------------------------------------------------------------

#define BENCH_RING_NUM_ATOM_RINGS(DATASET, SUFFIX, TAG)                        \
  TEST_CASE("ROMol numAtomRings all atoms " SUFFIX, "[mol_api]" TAG) {         \
    auto raw = load_samples(DATASET);                                          \
    auto samples = prime_samples(std::move(raw), [](ROMol &mol) {              \
      MolOps::findSSSR(mol);                                                   \
    });                                                                        \
    BENCHMARK_ADVANCED("ROMol numAtomRings all atoms " SUFFIX)                 \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        const auto *info = mol.getRingInfo();                                  \
        uint32_t total = 0;                                                    \
        for (auto atom : mol.atoms()) {                                        \
          total += info->numAtomRings(atom->getIdx());                         \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol numAtomRings all atoms " SUFFIX,                            \
            "[mol_api][rdmol]" TAG) {                                          \
    auto raw = load_rdmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [](RDMol &mol) {              \
      RingInfoCache &cache = mol.getRingInfo();                                \
      MolOps::findSSSR(mol, cache);                                            \
    });                                                                        \
    BENCHMARK_ADVANCED("RDMol numAtomRings all atoms " SUFFIX)                 \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const auto &info = mol.getRingInfo();                                  \
        const uint32_t numAtoms = mol.getNumAtoms();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numAtoms; ++i) {                              \
          total += info.numAtomRings(i);                                       \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

#define BENCH_RING_NUM_BOND_RINGS(DATASET, SUFFIX, TAG)                        \
  TEST_CASE("ROMol numBondRings all bonds " SUFFIX, "[mol_api]" TAG) {         \
    auto raw = load_samples(DATASET);                                          \
    auto samples = prime_samples(std::move(raw), [](ROMol &mol) {              \
      MolOps::findSSSR(mol);                                                   \
    });                                                                        \
    BENCHMARK_ADVANCED("ROMol numBondRings all bonds " SUFFIX)                 \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const ROMol &mol) {           \
        const auto *info = mol.getRingInfo();                                  \
        uint32_t total = 0;                                                    \
        for (auto bond : mol.bonds()) {                                        \
          total += info->numBondRings(bond->getIdx());                         \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }                                                                            \
  TEST_CASE("RDMol numBondRings all bonds " SUFFIX,                            \
            "[mol_api][rdmol]" TAG) {                                          \
    auto raw = load_rdmol_samples(DATASET);                                    \
    auto samples = prime_samples(std::move(raw), [](RDMol &mol) {              \
      RingInfoCache &cache = mol.getRingInfo();                                \
      MolOps::findSSSR(mol, cache);                                            \
    });                                                                        \
    BENCHMARK_ADVANCED("RDMol numBondRings all bonds " SUFFIX)                 \
    (Catch::Benchmark::Chronometer meter) {                                    \
      run_per_sample_readonly(samples, meter, [](const RDMol &mol) {           \
        const auto &info = mol.getRingInfo();                                  \
        const uint32_t numBonds = mol.getNumBonds();                           \
        uint32_t total = 0;                                                    \
        for (uint32_t i = 0; i < numBonds; ++i) {                              \
          total += info.numBondRings(i);                                       \
        }                                                                      \
        return uint64_t(total);                                                \
      });                                                                      \
    };                                                                         \
  }

// ---------------------------------------------------------------------------
// Bench instantiations across datasets
// ---------------------------------------------------------------------------

#define BENCH_API_FOR(DATASET, SUFFIX, TAG)                                    \
  BENCH_MOVE_CTOR(DATASET, SUFFIX, TAG)                                        \
  BENCH_MOVE_ASSIGN(DATASET, SUFFIX, TAG)                                      \
  BENCH_ATOM_ITER_SUM_ATOMICNUM(DATASET, SUFFIX, TAG)                          \
  BENCH_BOND_ITER_SUM_TWICEORDER(DATASET, SUFFIX, TAG)                         \
  BENCH_NEIGHBOR_WALK_SUM_DEGREE(DATASET, SUFFIX, TAG)                         \
  BENCH_BONDS_FROM_ATOM_SUM_ORDER(DATASET, SUFFIX, TAG)                        \
  BENCH_EXTRACT_ATOMICNUMS(DATASET, SUFFIX, TAG)                               \
  BENCH_EXTRACT_FORMAL_CHARGES(DATASET, SUFFIX, TAG)                           \
  BENCH_EXTRACT_MASSES(DATASET, SUFFIX, TAG)                                   \
  BENCH_EXTRACT_BOND_TYPES(DATASET, SUFFIX, TAG)                               \
  BENCH_PROP_SET_VARIED(DATASET, SUFFIX, TAG)                                  \
  BENCH_PROP_SET_UNIFORM(DATASET, SUFFIX, TAG)                                 \
  BENCH_PROP_SET_RANDOM_ACCESS(DATASET, SUFFIX, TAG)                           \
  BENCH_PROP_GET_RANDOM_ACCESS(DATASET, SUFFIX, TAG)                           \
  BENCH_PROP_HAS_PRESENT(DATASET, SUFFIX, TAG)                                 \
  BENCH_PROP_HAS_ABSENT(DATASET, SUFFIX, TAG)                                  \
  BENCH_PROP_CLEAR_BULK(DATASET, SUFFIX, TAG)                                  \
  BENCH_ACCESSOR_DEGREE(DATASET, SUFFIX, TAG)                                  \
  BENCH_ACCESSOR_NUM_IMPLICIT_HS(DATASET, SUFFIX, TAG)                         \
  BENCH_RING_NUM_ATOM_RINGS(DATASET, SUFFIX, TAG)                              \
  BENCH_RING_NUM_BOND_RINGS(DATASET, SUFFIX, TAG)

BENCH_API_FOR(Dataset::Canonical, "", "[canonical]")
BENCH_API_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_API_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_API_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")
BENCH_API_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
