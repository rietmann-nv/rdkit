#include <catch2/catch_all.hpp>
#include <string>

#include "bench_common.hpp"

#include <GraphMol/RDMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace RDKit;

TEST_CASE("SmilesToMol", "[smiles]") {
  BENCHMARK("SmilesToMol") {
    auto total_atoms = 0;
    for (auto smiles : bench_common::SAMPLES) {
      auto mol = v2::SmilesParse::MolFromSmiles(smiles);
      REQUIRE(mol);
      total_atoms += mol->getNumAtoms();
    }
    return total_atoms;
  };
}

TEST_CASE("SmilesToMol no-sanitize", "[smiles]") {
  v2::SmilesParse::SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;
  BENCHMARK("SmilesToMol no-sanitize") {
    auto total_atoms = 0;
    for (auto smiles : bench_common::SAMPLES) {
      auto mol = v2::SmilesParse::MolFromSmiles(smiles, params);
      REQUIRE(mol);
      total_atoms += mol->getNumAtoms();
    }
    return total_atoms;
  };
}

TEST_CASE("SmilesToMol RDMol no-sanitize", "[smiles][rdmol]") {
  SmilesParserParams params;
  params.sanitize = false;
  params.removeHs = false;
  BENCHMARK("SmilesToMol RDMol no-sanitize") {
    SmilesParseTemp temp;
    auto total_atoms = 0;
    for (auto smiles : bench_common::SAMPLES) {
      RDMol mol;
      const bool success = SmilesToMol(smiles, params, mol, temp);
      REQUIRE(success);
      total_atoms += mol.getNumAtoms();
    }
    return total_atoms;
  };
}

TEST_CASE("SmilesToMol RDMol", "[smiles][rdmol]") {
  // Default-params (sanitize=true, removeHs=true) end-to-end; this is the
  // SMILES path the entire migration was aiming at.
  SmilesParserParams params;
  BENCHMARK("SmilesToMol RDMol") {
    SmilesParseTemp temp;
    auto total_atoms = 0;
    for (auto smiles : bench_common::SAMPLES) {
      RDMol mol;
      const bool success = SmilesToMol(smiles, params, mol, temp);
      REQUIRE(success);
      total_atoms += mol.getNumAtoms();
    }
    return total_atoms;
  };
}

TEST_CASE("MolToSmiles", "[smiles]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolToSmiles") {
    auto total_length = 0;
    for (auto &mol : samples) {
      auto smiles = MolToSmiles(mol);
      total_length += smiles.size();
    }
    return total_length;
  };
}

TEST_CASE("MolToCXSmiles", "[smiles]") {
  auto samples = bench_common::load_samples();
  BENCHMARK("MolToCXSmiles") {
    auto total_length = 0;
    for (auto &mol : samples) {
      auto smiles = MolToCXSmiles(mol);
      total_length += smiles.size();
    }
    return total_length;
  };
}
