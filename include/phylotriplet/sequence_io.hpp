#pragma once
#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "sequence.hpp"

namespace phylotriplet {
namespace detail {

// -- String utilities
// ----------------------------------------------------------

inline std::string to_upper(std::string s) {
  for (char& c : s)
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

inline std::string trim(const std::string& s) {
  const size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  const size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

inline bool starts_with_icase(const std::string& s, const char* prefix) {
  const size_t plen = std::strlen(prefix);
  if (s.size() < plen) return false;
  for (size_t i = 0; i < plen; ++i)
    if (std::toupper(static_cast<unsigned char>(s[i])) !=
        std::toupper(static_cast<unsigned char>(prefix[i])))
      return false;
  return true;
}

// -- NEXUS comment stripping
// ---------------------------------------------------

/**
 * Replace NEXUS-style nested comments [...] with spaces, preserving newlines
 * so that line-based parsing still works correctly.
 */
inline std::string strip_nexus_comments(const std::string& s) {
  std::string r;
  r.reserve(s.size());
  int depth = 0;
  for (char c : s) {
    if (c == '[') {
      ++depth;
      r += ' ';
    } else if (c == ']' && depth) {
      --depth;
      r += ' ';
    } else if (depth) {
      r += (c == '\n' ? '\n' : ' ');
    } else {
      r += c;
    }
  }
  return r;
}

// -- NEXUS keyword parsers
// -----------------------------------------------------

/**
 * Parse an integer value from "KEY=<int>" in an already-upper-cased string.
 * Returns `def` if the key is absent.
 */
inline int parse_int_key(const std::string& upper, const std::string& key,
                         int def = 0) {
  const std::string needle = key + "=";
  const size_t p = upper.find(needle);
  if (p == std::string::npos) return def;
  try {
    return std::stoi(upper.substr(p + needle.size()));
  } catch (...) {
    return def;
  }
}

/**
 * Expand a NEXUS CHARSET range specification into a sorted vector of
 * 0-indexed column positions.
 *
 * Supports single positions ("5") and ranges ("1-100 201-300").
 * Step syntax (e.g. "1-100\3") is not supported.
 */
inline std::vector<size_t> parse_charset_ranges(const std::string& s) {
  std::vector<size_t> cols;
  std::istringstream ss(s);
  std::string tok;
  while (ss >> tok) {
    const size_t dash = tok.find('-');
    if (dash != std::string::npos) {
      const size_t lo = std::stoull(tok.substr(0, dash));
      const size_t hi = std::stoull(tok.substr(dash + 1));
      for (size_t i = lo; i <= hi; ++i) cols.push_back(i - 1);  // 0-indexed
    } else {
      cols.push_back(std::stoull(tok) - 1);
    }
  }
  std::sort(cols.begin(), cols.end());
  cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
  return cols;
}

// -- NEXUS name reading
// --------------------------------------------------------

/**
 * Read the next taxon name from `line` starting at `pos`, advancing `pos`
 * past the name.  Handles NEXUS single-quoted names ('' inside = literal ').
 */
inline std::string read_nexus_name(const std::string& line, size_t& pos) {
  // Skip leading whitespace
  while (pos < line.size() &&
         std::isspace(static_cast<unsigned char>(line[pos])))
    ++pos;

  if (pos >= line.size()) return "";

  if (line[pos] == '\'') {
    // Quoted name
    ++pos;
    std::string name;
    while (pos < line.size()) {
      if (line[pos] == '\'') {
        ++pos;
        if (pos < line.size() && line[pos] == '\'') {
          name += '\'';
          ++pos;  // escaped quote
        } else {
          break;
        }
      } else {
        name += line[pos++];
      }
    }
    return name;
  }

  // Unquoted: read until whitespace
  const size_t start = pos;
  while (pos < line.size() &&
         !std::isspace(static_cast<unsigned char>(line[pos])))
    ++pos;
  return line.substr(start, pos - start);
}

// -- NEXUS MATRIX parsers
// ------------------------------------------------------

/**
 * Parse MATRIX content (everything after the MATRIX keyword up to ";") for
 * sequential format: one line per taxon, name then sequence (spaces in the
 * sequence are stripped).  Handles sequences wrapped across multiple lines.
 */
inline std::vector<Sequence> parse_matrix_sequential(
    const std::vector<std::string>& lines, int ntax, int /*nchar*/) {
  std::vector<Sequence> seqs;
  std::map<std::string, size_t> name_idx;

  for (const auto& raw : lines) {
    const std::string line = trim(raw);
    if (line.empty()) continue;

    size_t pos = 0;
    const std::string name = read_nexus_name(line, pos);
    if (name.empty()) continue;

    // Determine taxon slot
    size_t idx;
    auto it = name_idx.find(name);
    if (it == name_idx.end()) {
      if (static_cast<int>(seqs.size()) >= ntax) continue;
      idx = seqs.size();
      name_idx[name] = idx;
      seqs.push_back({name, {}});
    } else {
      idx = it->second;  // continuation line (wrapped sequence)
    }

    // Append sequence data (strip any embedded spaces)
    while (pos < line.size() &&
           std::isspace(static_cast<unsigned char>(line[pos])))
      ++pos;
    for (size_t i = pos; i < line.size(); ++i)
      if (!std::isspace(static_cast<unsigned char>(line[i])))
        seqs[idx].bases.push_back(char_to_nucleotide(line[i]));
  }

  return seqs;
}

/**
 * Parse MATRIX content for interleaved format.
 *
 * The first block of ntax lines sets up taxon names.  In subsequent blocks,
 * if the first token of a line matches a known taxon name it is stripped;
 * otherwise the entire line (including the first token) is treated as
 * sequence data.  Blank lines are treated as block separators.
 */
inline std::vector<Sequence> parse_matrix_interleaved(
    const std::vector<std::string>& lines, int ntax, int /*nchar*/) {
  std::vector<Sequence> seqs;
  std::map<std::string, size_t> name_idx;
  bool first_block = true;
  int block_pos = 0;  // line index within current block

  for (const auto& raw : lines) {
    const std::string line = trim(raw);
    if (line.empty()) {
      if (block_pos > 0) {
        block_pos = 0;
        first_block = false;
      }
      continue;
    }

    size_t pos = 0;
    const std::string token = read_nexus_name(line, pos);
    if (token.empty()) continue;

    size_t idx;
    bool token_is_name;

    if (first_block) {
      if (static_cast<int>(seqs.size()) >= ntax) continue;
      idx = seqs.size();
      name_idx[token] = idx;
      seqs.push_back({token, {}});
      token_is_name = true;
    } else if (name_idx.count(token)) {
      idx = name_idx[token];
      token_is_name = true;  // program repeats names in subsequent blocks
    } else {
      // Subsequent block without repeated names: token is sequence data
      idx = static_cast<size_t>(block_pos);
      token_is_name = false;
      if (idx < seqs.size())
        for (char c : token) seqs[idx].bases.push_back(char_to_nucleotide(c));
    }

    // Append rest of sequence data on this line
    if (token_is_name) {
      while (pos < line.size() &&
             std::isspace(static_cast<unsigned char>(line[pos])))
        ++pos;
    }
    if (idx < seqs.size()) {
      for (size_t j = pos; j < line.size(); ++j)
        if (!std::isspace(static_cast<unsigned char>(line[j])))
          seqs[idx].bases.push_back(char_to_nucleotide(line[j]));
    }

    ++block_pos;
    if (block_pos >= ntax) {
      block_pos = 0;
      first_block = false;
    }
  }

  return seqs;
}

}  // namespace detail

// -- Public NEXUS API
// ----------------------------------------------------------

/**
 * Read a NEXUS file from an input stream.
 *
 * Supports:
 *  - Multiple DATA/CHARACTERS blocks (one per gene)
 *  - A single DATA block + SETS block with CHARSET commands (splits by locus)
 *  - Sequential and interleaved MATRIX formats
 *  - DNA and RNA datatypes (U treated as T)
 *  - NEXUS single-quoted taxon names
 *  - Nested comments [...]
 */
inline MultiGeneAlignment read_nexus(std::istream& stream,
                                     const std::string& source = "<stream>",
                                     bool filter_invariant = true) {
  // Read entire content and strip comments
  std::string content((std::istreambuf_iterator<char>(stream)),
                      std::istreambuf_iterator<char>());
  content = detail::strip_nexus_comments(content);

  if (detail::to_upper(content).find("#NEXUS") == std::string::npos)
    throw std::runtime_error("Not a NEXUS file (missing #NEXUS header): " +
                             source);

  // Split into semicolon-delimited statements
  std::vector<std::string> stmts;
  {
    std::string cur;
    for (char c : content) {
      if (c == ';') {
        cur = detail::trim(cur);
        if (!cur.empty()) stmts.push_back(cur);
        cur.clear();
      } else {
        cur += c;
      }
    }
  }

  MultiGeneAlignment result;
  // Ordered charset list to preserve SETS block ordering
  std::vector<std::pair<std::string, std::vector<size_t>>> charsets;
  int gene_count = 0;

  size_t i = 0;
  while (i < stmts.size()) {
    const std::string su = detail::to_upper(detail::trim(stmts[i]));
    const std::string upper_full = detail::to_upper(stmts[i]);

    if (su.find("BEGIN") != std::string::npos) {
      const std::string block_name =
          detail::trim(stmts[i].substr(upper_full.find("BEGIN") + 5));
      const std::string ublock = detail::to_upper(detail::trim(block_name));
      ++i;

      // Collect block body until END
      std::vector<std::string> body;
      while (i < stmts.size() &&
             detail::to_upper(detail::trim(stmts[i])) != "END") {
        body.push_back(stmts[i]);
        ++i;
      }

      if (ublock == "DATA" || ublock == "CHARACTERS") {
        ++gene_count;
        std::string gene_name = "gene_" + std::to_string(gene_count);
        int ntax = 0, nchar = 0;
        bool interleave = false;
        std::string matrix_content;

        for (const auto& stmt : body) {
          const std::string ss = detail::to_upper(detail::trim(stmt));
          const std::string su2 = detail::to_upper(stmt);

          if (detail::starts_with_icase(ss, "DIMENSIONS")) {
            const size_t kp = su2.find("DIMENSIONS");
            const std::string rest = su2.substr(kp + 10);
            ntax = detail::parse_int_key(rest, "NTAX");
            nchar = detail::parse_int_key(rest, "NCHAR");

          } else if (detail::starts_with_icase(ss, "FORMAT")) {
            const size_t kp = su2.find("FORMAT");
            const std::string upper_rest = su2.substr(kp + 6);
            if (upper_rest.find("INTERLEAVE") != std::string::npos)
              interleave = true;

            // Validate datatype
            const size_t dt = upper_rest.find("DATATYPE=");
            if (dt != std::string::npos) {
              std::string dtype = upper_rest.substr(dt + 9);
              dtype = dtype.substr(0, dtype.find_first_of(" \t\r\n;"));
              if (dtype != "DNA" && dtype != "RNA" && dtype != "NUCLEOTIDE")
                throw std::runtime_error("Unsupported NEXUS DATATYPE '" +
                                         dtype + "' in " + source +
                                         " (expected DNA, RNA, or NUCLEOTIDE)");
            }

          } else if (detail::starts_with_icase(ss, "TITLE")) {
            const size_t kp =
                stmt.find_first_not_of(" \t\r\n", su2.find("TITLE") + 5);
            if (kp != std::string::npos) {
              gene_name = detail::trim(stmt.substr(kp));
              if (gene_name.size() >= 2 && gene_name.front() == '\'' &&
                  gene_name.back() == '\'')
                gene_name = gene_name.substr(1, gene_name.size() - 2);
            }

          } else if (detail::starts_with_icase(ss, "MATRIX")) {
            const size_t kp = su2.find("MATRIX");
            matrix_content = stmt.substr(kp + 6);
          }
        }

        if (matrix_content.empty())
          throw std::runtime_error("DATA block with no MATRIX statement in " +
                                   source);
        if (ntax == 0 || nchar == 0)
          throw std::runtime_error("DATA block missing DIMENSIONS in " +
                                   source);

        std::vector<std::string> mlines;
        std::istringstream mss(matrix_content);
        std::string mline;
        while (std::getline(mss, mline)) mlines.push_back(mline);

        GeneAlignment gene;
        gene.gene_name = gene_name;
        gene.sequences =
            interleave ? detail::parse_matrix_interleaved(mlines, ntax, nchar)
                       : detail::parse_matrix_sequential(mlines, ntax, nchar);
        if (filter_invariant) gene.filter_invariant_sites();
        result.genes.push_back(std::move(gene));

      } else if (ublock == "SETS") {
        for (const auto& stmt : body) {
          const std::string ss = detail::to_upper(detail::trim(stmt));
          if (detail::starts_with_icase(ss, "CHARSET")) {
            const size_t kp = ss.find("CHARSET");
            std::string rest = detail::trim(stmt.substr(stmt.find_first_not_of(
                " \t\r\n", stmt.find("CHARSET") != std::string::npos
                               ? stmt.find("CHARSET") + 7
                               : 7)));
            const size_t eq = rest.find('=');
            if (eq == std::string::npos) continue;
            std::string csname = detail::trim(rest.substr(0, eq));
            // Unquote if needed
            if (csname.size() >= 2 && csname.front() == '\'' &&
                csname.back() == '\'')
              csname = csname.substr(1, csname.size() - 2);
            const std::string ranges = detail::trim(rest.substr(eq + 1));
            charsets.push_back({csname, detail::parse_charset_ranges(ranges)});
            (void)kp;
          }
        }
      }
      // Unknown blocks are skipped.
    }
    ++i;
  }

  // If there is exactly one data block AND a SETS block, split by charset.
  if (!charsets.empty() && result.genes.size() == 1) {
    const GeneAlignment& full = result.genes[0];
    MultiGeneAlignment split;
    for (const auto& cs : charsets) {
      GeneAlignment gene;
      gene.gene_name = cs.first;
      gene.sequences.resize(full.sequences.size());
      for (size_t t = 0; t < full.sequences.size(); ++t) {
        gene.sequences[t].taxon = full.sequences[t].taxon;
        for (size_t col : cs.second) {
          if (col < full.sequences[t].bases.size())
            gene.sequences[t].bases.push_back(full.sequences[t].bases[col]);
        }
      }
      if (filter_invariant) gene.filter_invariant_sites();
      split.genes.push_back(std::move(gene));
    }
    return split;
  }

  return result;
}

/** Read a NEXUS file by path. */
inline MultiGeneAlignment read_nexus_file(const std::string& filepath,
                                           bool filter_invariant = true) {
  std::ifstream f(filepath);
  if (!f.is_open()) throw std::runtime_error("Cannot open file: " + filepath);
  return read_nexus(f, filepath, filter_invariant);
}

// -- Public PHYLIP API
// ---------------------------------------------------------

/**
 * Read a single alignment from a PHYLIP-format stream.
 *
 * Supports relaxed PHYLIP (taxon name is the first whitespace-delimited token
 * on each data line; names may be any length).
 *
 * @param interleaved  Set true for interleaved format.  In the first block,
 *                     each line begins with the taxon name.  In subsequent
 *                     blocks the name is absent (strict PHYLIP) or present
 *                     (some programs); if the first token matches a known name
 *                     it is stripped, otherwise the whole line is sequence
 * data.
 */
inline GeneAlignment read_phylip(std::istream& stream,
                                 const std::string& gene_name = "gene",
                                 bool interleaved = false,
                                 const std::string& source = "<stream>",
                                 bool filter_invariant = true) {
  // Skip blank lines to find the header
  std::string header;
  while (std::getline(stream, header)) {
    header = detail::trim(header);
    if (!header.empty()) break;
  }
  if (header.empty()) throw std::runtime_error("Empty PHYLIP file: " + source);

  int ntax = 0, nchar = 0;
  {
    std::istringstream hss(header);
    if (!(hss >> ntax >> nchar) || ntax <= 0 || nchar <= 0)
      throw std::runtime_error("Invalid PHYLIP header '" + header + "' in " +
                               source);
  }

  GeneAlignment gene;
  gene.gene_name = gene_name;
  gene.sequences.resize(ntax);

  std::vector<std::string> lines;
  {
    std::string ln;
    while (std::getline(stream, ln)) lines.push_back(ln);
  }

  if (!interleaved) {
    // Sequential: one line per taxon
    int ti = 0;
    for (const auto& raw : lines) {
      const std::string line = detail::trim(raw);
      if (line.empty()) continue;
      if (ti >= ntax) break;
      std::istringstream lss(line);
      lss >> gene.sequences[ti].taxon;
      std::string chunk;
      while (lss >> chunk)
        for (char c : chunk)
          gene.sequences[ti].bases.push_back(char_to_nucleotide(c));
      ++ti;
    }
  } else {
    // Interleaved: first block has names; subsequent blocks have sequence only
    // (or repeated names — detected and stripped)
    std::map<std::string, int> name_idx;
    bool first_block = true;
    int block_pos = 0;

    for (const auto& raw : lines) {
      const std::string line = detail::trim(raw);
      if (line.empty()) {
        if (block_pos > 0) {
          block_pos = 0;
          first_block = false;
        }
        continue;
      }

      std::istringstream lss(line);
      std::string first_tok;
      lss >> first_tok;

      if (first_block) {
        gene.sequences[block_pos].taxon = first_tok;
        name_idx[first_tok] = block_pos;
      } else if (name_idx.count(first_tok)) {
        // Repeated name: strip it
        block_pos = name_idx[first_tok];
      } else {
        // No name: first_tok is sequence data
        for (char c : first_tok)
          gene.sequences[block_pos].bases.push_back(char_to_nucleotide(c));
      }

      // Append remaining tokens as sequence
      std::string chunk;
      while (lss >> chunk)
        for (char c : chunk)
          gene.sequences[block_pos].bases.push_back(char_to_nucleotide(c));

      ++block_pos;
      if (block_pos >= ntax) {
        block_pos = 0;
        first_block = false;
      }
    }
  }

  if (filter_invariant) gene.filter_invariant_sites();
  return gene;
}

/** Read a single PHYLIP alignment from a file.
 *  gene_name defaults to the filename stem (no extension). */
inline GeneAlignment read_phylip_file(const std::string& filepath,
                                      const std::string& gene_name = "",
                                      bool interleaved = false,
                                      bool filter_invariant = true) {
  std::ifstream f(filepath);
  if (!f.is_open()) throw std::runtime_error("Cannot open file: " + filepath);

  std::string name = gene_name;
  if (name.empty()) {
    size_t slash = filepath.find_last_of("/\\");
    name = (slash != std::string::npos) ? filepath.substr(slash + 1) : filepath;
    size_t dot = name.find_last_of('.');
    if (dot != std::string::npos) name = name.substr(0, dot);
  }
  return read_phylip(f, name, interleaved, filepath, filter_invariant);
}

/**
 * Read multiple PHYLIP files as a multi-gene alignment.
 *
 * @param filepaths   One file per gene.
 * @param gene_names  Optional names for each gene; defaults to filename stem.
 * @param interleaved Applied uniformly to all files.
 */
inline MultiGeneAlignment read_phylip_files(
    const std::vector<std::string>& filepaths,
    const std::vector<std::string>& gene_names = {}, bool interleaved = false,
    bool filter_invariant = true) {
  MultiGeneAlignment result;
  result.genes.reserve(filepaths.size());
  for (size_t i = 0; i < filepaths.size(); ++i) {
    const std::string& name = (i < gene_names.size()) ? gene_names[i] : "";
    result.genes.push_back(
        read_phylip_file(filepaths[i], name, interleaved, filter_invariant));
  }
  return result;
}

}  // namespace phylotriplet