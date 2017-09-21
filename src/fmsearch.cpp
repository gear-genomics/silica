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
  uint16_t prefix;
  uint16_t suffix;
  uint16_t maxmatch;
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
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "FM-index")
    ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.fa"), "output file")
    ("prefix,p", boost::program_options::value<uint16_t>(&c.prefix)->default_value(60), "prefix length")
    ("suffix,s", boost::program_options::value<uint16_t>(&c.suffix)->default_value(60), "suffix length")
    ("maxmatches,m", boost::program_options::value<uint16_t>(&c.maxmatch)->default_value(10), "max. number of matches")
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
  if ((vm.count("help")) || (!vm.count("input-file"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] sequences.fasta" << std::endl;
    std::cout << visible_options << "\n";
    return -1;
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Reference index
  csa_wt<> fm_index;
    
  // Load FM index
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load FM-Index" << std::endl;
  if (!load_from_file(fm_index, c.genome.string())) {
    std::cerr << "Index cannot be loaded!" << std::endl;
    return 1;
  }

  // Find reference matches
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Find Reference Match" << std::endl;

  // Done
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;

  return 0;
}
