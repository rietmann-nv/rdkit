#include <catch2/catch_all.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include <GraphMol/RDMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "bench_common.hpp"

namespace bench_common {

const char *dataset_name(Dataset dataset) {
  switch (dataset) {
    case Dataset::Canonical:
      return "canonical";
    case Dataset::Size_00_20:
      return "size_00_20";
    case Dataset::Size_20_40:
      return "size_20_40";
    case Dataset::Size_40_60:
      return "size_40_60";
    case Dataset::Size_60_80:
      return "size_60_80";
    case Dataset::Rings_2:
      return "rings_2";
    case Dataset::Rings_3:
      return "rings_3";
    case Dataset::Rings_4:
      return "rings_4";
    case Dataset::Rings_5:
      return "rings_5";
    case Dataset::Rings_6:
      return "rings_6";
    case Dataset::Radicals:
      return "radicals";
    case Dataset::Organometallics:
      return "organometallics";
    case Dataset::HypervalentHalogens:
      return "hypervalent_halogens";
    case Dataset::HypervalentP:
      return "hypervalent_p";
    case Dataset::PreCanonicalNO2Azide:
      return "pre_canonical_no2_azide";
    case Dataset::Atropisomers:
      return "atropisomers";
    case Dataset::KekulizeHard:
      return "kekulize_hard";
  }
  return "unknown";
}

namespace {

std::vector<std::string> load_dataset_file(Dataset dataset) {
  const char *rdbase_env = std::getenv("RDBASE");
  if (rdbase_env == nullptr) {
    throw std::runtime_error(
        "RDBASE environment variable not set; required to locate "
        "Code/Bench/data/*.smi");
  }
  std::filesystem::path path(rdbase_env);
  path /= "Code";
  path /= "Bench";
  path /= "data";
  path /= std::string(dataset_name(dataset)) + ".smi";

  std::ifstream file(path);
  if (!file) {
    std::ostringstream msg;
    msg << "could not open dataset file: " << path.string();
    throw std::runtime_error(msg.str());
  }

  std::vector<std::string> smiles;
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) {
      continue;
    }
    auto tab = line.find('\t');
    auto space = line.find(' ');
    auto cut = std::min(tab, space);
    if (cut == std::string::npos) {
      smiles.push_back(line);
    } else {
      smiles.emplace_back(line.substr(0, cut));
    }
  }
  return smiles;
}

const std::vector<std::string> &cached_dataset_smiles(Dataset dataset) {
  static std::mutex mutex;
  static std::unordered_map<int, std::vector<std::string>> cache;
  std::scoped_lock lock(mutex);
  const int key = static_cast<int>(dataset);
  auto it = cache.find(key);
  if (it != cache.end()) {
    return it->second;
  }
  auto inserted = cache.emplace(key, load_dataset_file(dataset));
  return inserted.first->second;
}

const std::vector<std::string> &canonical_smiles() {
  static const std::vector<std::string> result = [] {
    std::vector<std::string> ret;
    ret.reserve(std::size(SAMPLES));
    for (auto smiles : SAMPLES) {
      ret.emplace_back(smiles);
    }
    return ret;
  }();
  return result;
}

}  // namespace

const std::vector<std::string> &dataset_smiles(Dataset dataset) {
  if (dataset == Dataset::Canonical) {
    return canonical_smiles();
  }
  return cached_dataset_smiles(dataset);
}

std::vector<RDKit::ROMol> load_samples() {
  return load_samples(Dataset::Canonical);
}

std::vector<RDKit::RDMol> load_rdmol_samples() {
  return load_rdmol_samples(Dataset::Canonical);
}

std::vector<RDKit::ROMol> load_samples(Dataset dataset, bool sanitize) {
  RDKit::v2::SmilesParse::SmilesParserParams params;
  params.sanitize = sanitize;
  params.removeHs = sanitize;

  const auto &smiles_list = dataset_smiles(dataset);
  std::vector<RDKit::ROMol> ret;
  ret.reserve(smiles_list.size());
  for (const auto &smiles : smiles_list) {
    auto mol = RDKit::v2::SmilesParse::MolFromSmiles(smiles, params);
    REQUIRE(mol);
    ret.push_back(std::move(*mol));
  }
  return ret;
}

std::vector<RDKit::RDMol> load_rdmol_samples(Dataset dataset, bool sanitize) {
  RDKit::SmilesParserParams params;
  params.sanitize = sanitize;
  params.removeHs = sanitize;
  RDKit::SmilesParseTemp temp;

  const auto &smiles_list = dataset_smiles(dataset);
  std::vector<RDKit::RDMol> ret;
  ret.reserve(smiles_list.size());
  for (const auto &smiles : smiles_list) {
    auto &mol = ret.emplace_back();
    const bool success = RDKit::SmilesToMol(smiles.c_str(), params, mol, temp);
    REQUIRE(success);
    // The default-params SmilesToMol(RDMol&) path is now native end to
    // end on canonical input (commits 297ea7237..97657fc7a). The four
    // remaining bridges (Canon::rankFragmentAtoms, Canon::rankMolAtoms,
    // and the CXSMILES stereo helpers) only fire on inputs the canonical
    // bench set doesn't contain. As a defensive guarantee that the
    // benchmarks read pristine flat-array RDMols, drop any compat block
    // that an unaudited path might have allocated.
    mol.clearCompatibilityData();
  }
  return ret;
}

}  // namespace bench_common
