#include <catch2/catch_all.hpp>

#include <GraphMol/RDMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "bench_common.hpp"

namespace bench_common {

std::vector<RDKit::ROMol> load_samples() {
  std::vector<RDKit::ROMol> ret;
  for (auto smiles : SAMPLES) {
    auto mol = RDKit::v2::SmilesParse::MolFromSmiles(smiles);
    REQUIRE(mol);
    ret.push_back(std::move(*mol));
  }
  return ret;
}

std::vector<RDKit::RDMol> load_rdmol_samples() {
  RDKit::SmilesParserParams params;
  RDKit::SmilesParseTemp temp;

  std::vector<RDKit::RDMol> ret;
  ret.reserve(std::size(SAMPLES));
  for (auto smiles : SAMPLES) {
    auto &mol = ret.emplace_back();
    const bool success = RDKit::SmilesToMol(smiles, params, mol, temp);
    REQUIRE(success);
    // Sanitization in the native RDMol parser path goes through asRWMol(),
    // which materializes a CompatibilityData block. Discard it so the
    // benchmarks that follow measure the flat-array path only.
    mol.clearCompatibilityData();
  }
  return ret;
}

}  // namespace bench_common
