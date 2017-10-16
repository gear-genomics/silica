/*
============================================================================
FMsearch: FM-Index Search
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
#include <vector>

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
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

#include "neighbors.h"
#include "align.h"
#include "needle.h"
#include "thal.h"

using namespace sdsl;
using namespace fmsearch;

struct Config {
  bool indel;
  bool align;
  uint32_t cutTemp;
  uint32_t maxProdSize;
  uint32_t kmer;
  uint32_t distance;
  std::size_t pre_context;
  std::size_t post_context;
  std::size_t max_locations;
  boost::filesystem::path outfile;
  boost::filesystem::path primfile;
  boost::filesystem::path resfile;
  boost::filesystem::path infile;
  boost::filesystem::path genome;
};

struct primerBind {
  std::string chrom;
  uint32_t pos;
  bool onFor;
  double temp;
  std::string primer;
  uint32_t leng;
  std::string genSeq;
};

struct sortPrimer
{
    inline bool operator() (const primerBind& a, const primerBind& b)
    {
        return (a.temp > b.temp);
    }
};

struct pcrProduct {
  std::string chrom;
  uint32_t leng;
  std::string seq;
  uint32_t forPos;
  double forTemp;
  std::string forPrimer;
  uint32_t revPos;
  double revTemp;
  std::string revPrimer;
};

struct sortProducts
{
    inline bool operator() (const pcrProduct& a, const pcrProduct& b)
    {
        return (a.leng < b.leng);
    }
};



void addUniqe(std::vector<primerBind>* coll, primerBind* prim) {
  int found = 0;
  for(std::vector<primerBind>::iterator it = coll->begin(); it != coll->end(); ++it) {
    if ((prim->chrom.compare(it->chrom) == 0) &&
        ((prim->pos - it->pos == 2) ||
         (prim->pos - it->pos == 1) ||
         (prim->pos - it->pos == 0) || // Better keep track during edit distance generation and compensate
         (it->pos - prim->pos == 1) ||
         (it->pos - prim->pos == 2)) && 
        (prim->primer.compare(it->primer) == 0)){
      found = 1;
    }
  }
  if (found == 0) {
    coll->push_back(*prim);
  }
}


int main(int argc, char** argv) {
  Config c;
  // CMD Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
    ;

  boost::program_options::options_description appr("Approximate Search Options");
  appr.add_options()
    ("cutTemp,x", boost::program_options::value<uint32_t>(&c.cutTemp)->default_value(35), "minimal primer meting temperature to consider")
    ("maxProdSize,u", boost::program_options::value<uint32_t>(&c.maxProdSize)->default_value(15000), "maximal PCR Product size")
    ("kmer,k", boost::program_options::value<uint32_t>(&c.kmer)->default_value(15), "k-mer size")
    ("maxmatches,m", boost::program_options::value<std::size_t>(&c.max_locations)->default_value(10000), "max. number of matches per k-mer")
    ("distance,d", boost::program_options::value<uint32_t>(&c.distance)->default_value(1), "neighborhood distance")
    ("hamming,n", "use hamming neighborhood instead of edit distance")
    ;

  boost::program_options::options_description outp("Output Options");
  outp.add_options()
    ("prefix,p", boost::program_options::value<std::size_t>(&c.pre_context)->default_value(3), "prefix length")
    ("suffix,s", boost::program_options::value<std::size_t>(&c.post_context)->default_value(3), "suffix length")
    ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.fa"), "output file")
    ("result,r", boost::program_options::value<boost::filesystem::path>(&c.resfile)->default_value("res.txt"), "result file")
    ("primer,p", boost::program_options::value<boost::filesystem::path>(&c.primfile)->default_value("primer.fa"), "primer locations file")
    ("align,a", "write alignments to stderr")
    ;
    
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "seq.fasta")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(appr).add(outp).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(appr).add(outp);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -g <ref.fa.gz> sequences.fasta" << std::endl;
    std::cout << visible_options << "\n";
    return -1;
  }

  // Cmd switches
  if (!vm.count("align")) c.align = false;
  else c.align = true;
  if (!vm.count("hamming")) c.indel = true;
  else c.indel = false;

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

  // Define alphabet
  typedef std::set<char> TAlphabet;
  char tmp[] = {'A', 'C', 'G', 'T'};
  TAlphabet alphabet(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));

  // Initialize thal arguments
  boost::filesystem::path exepath = boost::filesystem::system_complete(argv[0]).parent_path();
  std::string cfgpath = exepath.string() + "/primer3_config/";
  primer3thal::thal_args a;
  primer3thal::set_thal_default_args(&a);
  a.temponly=1;
  a.type = primer3thal::thal_end1;
  a.mv = 50.0;
  a.dv =  1.5;
  a.dna_conc = 50.0;
  a.dntp = 0.6;
  a.temp = 50.0 + primer3thal::ABSOLUTE_ZERO;
  primer3thal::get_thermodynamic_values(cfgpath.c_str());

  // Query FM-Index
  std::ofstream ofile(c.outfile.string().c_str());
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Query FM-Index" << std::endl;
  std::vector<primerBind> forBind;
  std::vector<primerBind> revBind;  
  for(TInputFasta::const_iterator itFa = infa.begin(); itFa != infa.end(); ++itFa) {
    std::string qr = itFa->second;
    if (qr.size() < c.kmer) continue;
    int32_t koffset = qr.size() - c.kmer;
    qr = qr.substr(qr.size() - c.kmer);
    typedef std::set<std::string> TStringSet;
    TStringSet fwdset;
    neighbors(qr, alphabet, c.distance, c.indel, fwdset);
    // Debug
    //for(TStringSet::iterator it = fwdset.begin(); it != fwdset.end(); ++it) std::cerr << *it << std::endl;
    TStringSet revset;
    reverseComplement(qr);
    neighbors(qr, alphabet, c.distance, c.indel, revset);
    int32_t qhits = 0;
    for(int32_t fwdrev = 0; fwdrev < 2; ++fwdrev) {
      TStringSet::iterator its;
      TStringSet::iterator itsEnd;
      if (fwdrev == 0) {
	its = fwdset.begin();
	itsEnd = fwdset.end();
      } else {
	its = revset.begin();
	itsEnd = revset.end();
      }
      for(; its != itsEnd; ++its, ++qhits) {
	std::string query(*its);
	std::size_t m = query.size();
	std::size_t occs = sdsl::count(fm_index, query.begin(), query.end());
	if (occs > 0) {
	  auto locations = locate(fm_index, query.begin(), query.begin() + m);
	  std::sort(locations.begin(), locations.end());
	  std::size_t pre_extract = c.pre_context;
	  std::size_t post_extract = c.post_context;
	  if (fwdrev == 0) pre_extract += koffset;
	  else post_extract += koffset;
	  for(std::size_t i = 0; i < std::min(occs, c.max_locations); ++i) {
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

	    // thermodynamical alignemnt
	    std::string genomicseq = pre + s.substr(0, m) + post;
	    std::string primer = itFa->second;
	    if (fwdrev == 0) reverseComplement(primer);
	    primer3thal::oligo1 = (unsigned char*) primer.c_str();
	    primer3thal::oligo2 = (unsigned char*) genomicseq.c_str();
	    primer3thal::thal_results o;
	    bool thalsuccess = primer3thal::thal(primer3thal::oligo1, primer3thal::oligo2, &a, &o);
	    if ((!thalsuccess) || (o.temp == primer3thal::THAL_ERROR_SCORE)) {
	      std::cerr << "Error during thermodynamical calculation!" << std::endl;
	      return -1;
	    }

            // Score suitable primers
            primerBind prim;
            prim.chrom = std::string(faidx_iseq(fai, refIndex));
            prim.pos = chrpos;
            prim.temp = o.temp;
            prim.primer = itFa->first;
            prim.leng = 0;
            prim.genSeq = genomicseq;
            if (o.temp > c.cutTemp) {
              if (fwdrev == 0) {
                 prim.onFor = true;
                //forBind.push_back(prim);
                addUniqe(&forBind, &prim);
              } else {
                prim.onFor = false;
                //revBind.push_back(prim);
                addUniqe(&revBind, &prim);
              }
            }
                           
	    // Output fasta
	    ofile << ">" << std::string(faidx_iseq(fai, refIndex)) << ":" << chrpos;
	    if (fwdrev == 0) ofile << " fwd";
	    else ofile << " rev";
	    ofile << " temp:" << o.temp;
	    ofile << " hit:(" << qhits << "," << i << ") " << itFa->first << std::endl;
	    ofile << genomicseq << std::endl;

	    // Debug alignment output
	    if (c.align) {
	      typedef boost::multi_array<char, 2> TAlign;
	      TAlign align;
	      AlignConfig<true, false> semiglobal;
	      DnaScore<int> sc(5, -4, -4, -4);
	      int32_t score = 0;
	      primer = itFa->second;
	      if (fwdrev == 0) score = needle(primer, genomicseq, align, semiglobal, sc);
	      else {
		reverseComplement(primer);
		score = needle(primer, genomicseq, align, semiglobal, sc);
	      }
	      std::cerr << ">" << std::string(faidx_iseq(fai, refIndex)) << ":" << chrpos;
	      if (fwdrev == 0) std::cerr << " fwd";
	      else std::cerr << " rev";
	      std::cerr << " temp:" << o.temp;
	      std::cerr << " hit:(" << qhits << "," << i << ") " << itFa->first << std::endl;
	      std::cerr << "AlignScore: " << score << std::endl;
	      typedef typename TAlign::index TAIndex;
	      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
		for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
		  std::cerr << align[i][j];
		}
		std::cerr << std::endl;
	      }
	      std::cerr << std::endl;
	    }	    
	  }
	}
      }
    }
  }
  ofile.close();

  // Find PCR products
  pcrProduct pcrProd;
  std::vector<pcrProduct> pcrColl;
  for(std::vector<primerBind>::iterator fw = forBind.begin(); fw != forBind.end(); ++fw) {
    for(std::vector<primerBind>::iterator rv = revBind.begin(); rv != revBind.end(); ++rv) {
      if ((fw->chrom.compare(rv->chrom) == 0) &&
          (rv->pos - fw->pos < c.maxProdSize)){
        pcrProd.leng = rv->pos - fw->pos;
        pcrProd.chrom = fw->chrom;
        pcrProd.forPos = fw->pos;
        pcrProd.forTemp = fw->temp;
        pcrProd.forPrimer = fw->primer;
        pcrProd.revPos = rv->pos;
        pcrProd.revTemp = rv->temp;
        pcrProd.revPrimer = rv->primer;
        pcrColl.push_back(pcrProd);
      }
    }
  }

  std::sort(pcrColl.begin(), pcrColl.end(), sortProducts());

  std::ofstream rfile(c.resfile.string().c_str());
  int count = 0;
  for(std::vector<pcrProduct>::iterator it = pcrColl.begin(); it != pcrColl.end(); ++it) {
    rfile << "PCR Product " << ++count << std::endl;
    rfile << "  Size: " << it->leng << std::endl;
    rfile << "  Chrom: " << it->chrom << std::endl;
    rfile << "  Forward Primer: " << it->forPrimer << std::endl;
    rfile << "    Melting Temperature: " << it->forTemp << " C" << std::endl;
    rfile << "    Position: " << it->forPos << std::endl;
    rfile << "  Reverse Primer: " << it->revPrimer << std::endl;
    rfile << "    Melting Temperature: " << it->revTemp << " C" << std::endl;
    rfile << "    Position: " << it->revPos << std::endl << std::endl;   
  }
  rfile.close();

  //Print Primers
  forBind.insert( forBind.end(), revBind.begin(), revBind.end() );
  std::sort(forBind.begin(), forBind.end(), sortPrimer());
  std::ofstream forfile(c.primfile.string().c_str());
  for(std::vector<primerBind>::iterator it = forBind.begin(); it != forBind.end(); ++it) {
    forfile << ">" << it->chrom << ":" << it->pos;
    if (it->onFor) forfile << " fwd";
    else forfile << " rev";
    forfile << " temp:" << it->temp;
    forfile << " primer:" << it->primer << std::endl;
    forfile << it->genSeq << std::endl;
  }
  forfile.close();

  // Clean-up
  primer3thal::destroy_thal_structures();
  if (fai != NULL) fai_destroy(fai);
  
  // Done
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;

  return 0;
}
