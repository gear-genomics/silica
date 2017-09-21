/*
============================================================================
Indigo: InDel Discovery in Sanger Chromatograms
============================================================================
Copyright (C) 2017 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/


#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <sdsl/suffix_arrays.hpp>
#include <htslib/faidx.h>

using namespace sdsl;

struct Config {
  std::size_t pre_context;
  std::size_t post_context;
  std::size_t max_locations;
  boost::filesystem::path outfile;
  boost::filesystem::path infile;
  boost::filesystem::path genome;
};

int main(int argc, char** argv) {
  Config c;
  
  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "Genome file")
    ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.fa"), "output file")
    ("prefix,p", boost::program_options::value<std::size_t>(&c.pre_context)->default_value(10), "prefix length")
    ("suffix,s", boost::program_options::value<std::size_t>(&c.post_context)->default_value(10), "suffix length")
    ("maxmatches,m", boost::program_options::value<std::size_t>(&c.max_locations)->default_value(5), "max. number of matches")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "seq.fasta")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -g <ref.fa.gz> sequences.fasta" << std::endl;
    std::cout << visible_options << "\n";
    return -1;
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parse chromosomes" << std::endl;

  // Parse chromosome lengths
  std::vector<uint32_t> seqlen;
  if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
    std::cerr << "Input reference file is missing: " << c.genome.string() << std::endl;
    return 1;
  }
  faidx_t* fai = fai_load(c.genome.string().c_str());
  if (fai == NULL) {
    if (fai_build(c.genome.string().c_str()) == -1) {
      std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
      return 1;
    } else fai = fai_load(c.genome.string().c_str());
  }
  seqlen.resize(faidx_nseq(fai));
  for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
    std::string seqname(faidx_iseq(fai, refIndex));
    seqlen[refIndex] = faidx_seq_len(fai, seqname.c_str()) + 1;
  }
  
  // Reference index
  csa_wt<> fm_index;  
  
  // Load FM index
  boost::filesystem::path op = c.genome.parent_path() / c.genome.stem();
  std::string index_file = op.string() + ".fm9";
  
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load FM-Index" << std::endl;
  if (!load_from_file(fm_index, index_file)) {
    std::cerr << "Index cannot be loaded!" << std::endl;
    return 1;
  }

  // Parse input fasta
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Parse input FASTA" << std::endl;

  if (!(boost::filesystem::exists(c.infile) && boost::filesystem::is_regular_file(c.infile) && boost::filesystem::file_size(c.infile))) {
    std::cerr << "Input fasta file is missing: " << c.infile.string() << std::endl;
    return 1;
  }
  typedef std::map<std::string, std::string> TInputFasta;
  TInputFasta infa;
  std::ifstream fafile(c.infile.string().c_str());
  if (fafile.good()) {
    std::string fan;
    std::string tmpfasta;
    std::string line;
    while(std::getline(fafile, line)) {
      if (!line.empty()) {
	if (line[0] == '>') {
	  if ((!fan.empty()) && (!tmpfasta.empty())) {
	    infa.insert(std::make_pair(fan, tmpfasta));
	    tmpfasta = "";
	  }
	  fan = line.substr(1);
	} else {
	  tmpfasta += boost::to_upper_copy(line);
	}
      }
    }
    if ((!fan.empty()) && (!tmpfasta.empty())) {
      infa.insert(std::make_pair(fan, tmpfasta));
    }
  }

  // Query FM-Index
  std::ofstream ofile(c.outfile.string().c_str());
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Query FM-Index" << std::endl;
  for(TInputFasta::const_iterator itFa = infa.begin(); itFa != infa.end(); ++itFa) {
    std::string query = itFa->second;
    std::size_t m = query.size();
    std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
    if (occs > 0) {
      auto locations = locate(fm_index, query.begin(), query.begin() + m);
      std::sort(locations.begin(), locations.end());
      for(std::size_t i = 0, pre_extract = c.pre_context, post_extract = c.post_context; i < std::min(occs, c.max_locations); ++i) {
	int64_t bestPos = locations[i];
	int64_t cumsum = 0;
	uint32_t refIndex = 0;
	for(; bestPos >= cumsum + seqlen[refIndex]; ++refIndex) cumsum += seqlen[refIndex];
	uint32_t chrpos = bestPos - cumsum;
	if (pre_extract > locations[i]) {
	  pre_extract = locations[i];
	}
	if (locations[i]+m+post_extract > fm_index.size()) {
	  post_extract = fm_index.size() - locations[i] - m;
	}
	auto s = extract(fm_index, locations[i]-pre_extract, locations[i]+m+post_extract-1);
	std::string pre = s.substr(0, pre_extract);
	s = s.substr(pre_extract);
	if (pre.find_last_of('\n') != std::string::npos) {
	  pre = pre.substr(pre.find_last_of('\n')+1);
	}
	std::string post = s.substr(m);
	post = post.substr(0, post.find_first_of('\n'));
	std::string fullseq = pre + s.substr(0, m) + post;
	ofile << ">" << std::string(faidx_iseq(fai, refIndex)) << ":" << chrpos << " hit:" << (i+1) << " " << itFa->first << std::endl;
	ofile << fullseq << std::endl;
      }
    }
  }
  ofile.close();

  // Clean-up
  if (fai != NULL) fai_destroy(fai);
  
  // Done
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;

  return 0;
}
