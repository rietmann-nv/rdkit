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
    // The default-params SmilesToMol(RDMol&) path is now native end to
    // end (commits 297ea7237..d9e052c7e), but the
    // SANITIZE_CLEANUP_ORGANOMETALLICS step still bridges through
    // asRWMol() to call Canon::rankMolAtoms (deferred Canon module
    // port). That bridge allocates a CompatibilityData block. Drop it
    // here so subsequent benchmarks measure the flat-array path
    // without any leftover compat overhead.
    mol.clearCompatibilityData();
  }
  return ret;
}

}  // namespace bench_common
