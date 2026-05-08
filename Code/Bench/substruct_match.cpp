#include <catch2/catch_all.hpp>
#include <filesystem>
#include <istream>
#include <string>
#include <vector>

#include "bench_common.hpp"

#include <GraphMol/RDMol.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

using bench_common::Dataset;

namespace bench_substruct_match {

std::filesystem::path relative_to_rdbase(
    const std::filesystem::path &relative) {
  char *rdbase = std::getenv("RDBASE");
  if (!rdbase) {
    throw std::runtime_error("RDBASE environment variable not set");
  }
  std::filesystem::path path(rdbase);
  path /= relative;
  return path;
}

std::vector<std::string> parse_smarts_examples(std::istream &in) {
  std::vector<std::string> examples;
  std::string line;
  while (std::getline(in, line)) {
    if (line.starts_with("#") || line.starts_with(" ") || line.empty()) {
      continue;
    }
    std::istringstream ss(line);
    std::string smarts;
    std::getline(ss, smarts, ' ');
    examples.push_back(line);
  }
  return examples;
}

std::vector<std::string> load_smarts_examples(std::filesystem::path &path) {
  std::ifstream file(path);
  return parse_smarts_examples(file);
}

std::vector<ROMol> load_smarts_queries(const std::filesystem::path &relative) {
  auto path = relative_to_rdbase(relative);
  auto examples = load_smarts_examples(path);
  std::vector<ROMol> queries;
  for (const auto &example : examples) {
    auto query = v2::SmilesParse::MolFromSmarts(example);
    REQUIRE(query);
    queries.push_back(std::move(*query));
  }
  return queries;
}

}  // namespace bench_substruct_match

TEST_CASE("ROMol::GetSubstructMatch", "[substruct_match]") {
  auto samples = bench_common::load_samples();

  auto query = v2::SmilesParse::MolFromSmarts(
      "[#6,#7]1:[#6,#7]:[#6](:[#6,#7]:[#6,#7]:[#6]:1-[#17,#35,#53,#9])-[#5](-[#8])-[#8]");
  BENCHMARK("ROMol::GetSubstructMatch") {
    auto total = 0;
    for (auto &mol : samples) {
      auto matches = SubstructMatch(mol, *query);
      total += matches.size();
    }
    return total;
  };
}

TEST_CASE("RDMol::GetSubstructMatch", "[substruct_match][rdmol]") {
  auto samples = bench_common::load_rdmol_samples();

  auto query = v2::SmilesParse::MolFromSmarts(
      "[#6,#7]1:[#6,#7]:[#6](:[#6,#7]:[#6,#7]:[#6]:1-[#17,#35,#53,#9])-[#5](-[#8])-[#8]");
  REQUIRE(query);
  const RDMol &queryRDMol = query->asRDMol();
  BENCHMARK("RDMol::GetSubstructMatch") {
    auto total = 0;
    for (const auto &mol : samples) {
      auto matches = SubstructMatch(mol, queryRDMol);
      total += matches.size();
    }
    return total;
  };
}

TEST_CASE("ROMol::GetSubstructMatch RLewis", "[substruct_match]") {
  auto mols = bench_common::load_samples();
  auto queries = bench_substruct_match::load_smarts_queries(
      "Data/SmartsLib/RLewis_smarts.txt");

  BENCHMARK("ROMol::GetSubstructMatch RLewis") {
    auto total = 0;
    for (const auto &mol : mols) {
      for (const auto &query : queries) {
        auto matches = SubstructMatch(mol, query);
        total += matches.size();
      }
    }
    return total;
  };
}

TEST_CASE("RDMol::GetSubstructMatch RLewis", "[substruct_match][rdmol]") {
  auto mols = bench_common::load_rdmol_samples();
  auto queries = bench_substruct_match::load_smarts_queries(
      "Data/SmartsLib/RLewis_smarts.txt");
  std::vector<const RDMol *> queryRDMols;
  queryRDMols.reserve(queries.size());
  for (const auto &q : queries) {
    queryRDMols.push_back(&q.asRDMol());
  }

  BENCHMARK("RDMol::GetSubstructMatch RLewis") {
    auto total = 0;
    for (const auto &mol : mols) {
      for (const auto *query : queryRDMols) {
        auto matches = SubstructMatch(mol, *query);
        total += matches.size();
      }
    }
    return total;
  };
}

TEST_CASE("ROMol::GetSubstructMatch patty", "[substruct_match]") {
  auto mols = bench_common::load_samples();
  auto queries = bench_substruct_match::load_smarts_queries(
      "Data/SmartsLib/patty_rules.txt");

  BENCHMARK("ROMol::GetSubstructMatch patty") {
    auto total = 0;
    for (const auto &mol : mols) {
      for (const auto &query : queries) {
        auto matches = SubstructMatch(mol, query);
        total += matches.size();
      }
    }
    return total;
  };
}

TEST_CASE("RDMol::GetSubstructMatch patty", "[substruct_match][rdmol]") {
  auto mols = bench_common::load_rdmol_samples();
  auto queries = bench_substruct_match::load_smarts_queries(
      "Data/SmartsLib/patty_rules.txt");
  std::vector<const RDMol *> queryRDMols;
  queryRDMols.reserve(queries.size());
  for (const auto &q : queries) {
    queryRDMols.push_back(&q.asRDMol());
  }

  BENCHMARK("RDMol::GetSubstructMatch patty") {
    auto total = 0;
    for (const auto &mol : mols) {
      for (const auto *query : queryRDMols) {
        auto matches = SubstructMatch(mol, *query);
        total += matches.size();
      }
    }
    return total;
  };
}

#define BENCH_SUBSTRUCT(DATASET, SUFFIX, TAG)                                  \
  TEST_CASE("ROMol::GetSubstructMatch " SUFFIX, "[substruct_match]" TAG) {     \
    auto samples = bench_common::load_samples(DATASET);                        \
    auto query = v2::SmilesParse::MolFromSmarts(                               \
        "[#6,#7]1:[#6,#7]:[#6](:[#6,#7]:[#6,#7]:[#6]:1-[#17,#35,#53,#9])-[#"   \
        "5](-[#8])-[#8]");                                                     \
    BENCHMARK("ROMol::GetSubstructMatch " SUFFIX) {                            \
      auto total = 0;                                                          \
      for (auto &mol : samples) {                                              \
        auto matches = SubstructMatch(mol, *query);                            \
        total += matches.size();                                               \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_SUBSTRUCT_RLEWIS(DATASET, SUFFIX, TAG)                           \
  TEST_CASE("ROMol::GetSubstructMatch RLewis " SUFFIX,                         \
            "[substruct_match]" TAG) {                                         \
    auto mols = bench_common::load_samples(DATASET);                           \
    auto queries = bench_substruct_match::load_smarts_queries(                 \
        "Data/SmartsLib/RLewis_smarts.txt");                                   \
    BENCHMARK("ROMol::GetSubstructMatch RLewis " SUFFIX) {                     \
      auto total = 0;                                                          \
      for (const auto &mol : mols) {                                           \
        for (const auto &query : queries) {                                    \
          auto matches = SubstructMatch(mol, query);                           \
          total += matches.size();                                             \
        }                                                                      \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_SUBSTRUCT_PATTY(DATASET, SUFFIX, TAG)                            \
  TEST_CASE("ROMol::GetSubstructMatch patty " SUFFIX,                          \
            "[substruct_match]" TAG) {                                         \
    auto mols = bench_common::load_samples(DATASET);                           \
    auto queries = bench_substruct_match::load_smarts_queries(                 \
        "Data/SmartsLib/patty_rules.txt");                                     \
    BENCHMARK("ROMol::GetSubstructMatch patty " SUFFIX) {                      \
      auto total = 0;                                                          \
      for (const auto &mol : mols) {                                           \
        for (const auto &query : queries) {                                    \
          auto matches = SubstructMatch(mol, query);                           \
          total += matches.size();                                             \
        }                                                                      \
      }                                                                        \
      return total;                                                            \
    };                                                                         \
  }

#define BENCH_SUBSTRUCT_FOR(DATASET, SUFFIX, TAG)                              \
  BENCH_SUBSTRUCT(DATASET, SUFFIX, TAG)                                        \
  BENCH_SUBSTRUCT_RLEWIS(DATASET, SUFFIX, TAG)                                 \
  BENCH_SUBSTRUCT_PATTY(DATASET, SUFFIX, TAG)

BENCH_SUBSTRUCT_FOR(Dataset::Size_00_20, "size 00-20", "[size_00_20]")
BENCH_SUBSTRUCT_FOR(Dataset::Size_20_40, "size 20-40", "[size_20_40]")
BENCH_SUBSTRUCT_FOR(Dataset::Size_40_60, "size 40-60", "[size_40_60]")
BENCH_SUBSTRUCT_FOR(Dataset::Size_60_80, "size 60-80", "[size_60_80]")

BENCH_SUBSTRUCT_FOR(Dataset::Rings_2, "rings 2", "[rings_2]")
BENCH_SUBSTRUCT_FOR(Dataset::Rings_3, "rings 3", "[rings_3]")
BENCH_SUBSTRUCT_FOR(Dataset::Rings_4, "rings 4", "[rings_4]")
BENCH_SUBSTRUCT_FOR(Dataset::Rings_5, "rings 5", "[rings_5]")
BENCH_SUBSTRUCT_FOR(Dataset::Rings_6, "rings 6", "[rings_6]")
