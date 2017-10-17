/*
============================================================================
Silica: In-silico PCR
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
using namespace silica;

struct Config {
  bool indel;
  bool align;
  double cutTemp;
  uint32_t maxProdSize;
  double targetTemp;
  double cutofPen;
  double penDiff;
  double penMis;
  double penLen;
  uint32_t kmer;
  uint32_t distance;
  std::size_t pre_context;
  std::size_t post_context;
  std::size_t max_locations;
  boost::filesystem::path outfile;
  boost::filesystem::path primfile;
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
  std::string seq;
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
  double penalty;
  std::string seq;
  uint32_t forPos;
  double forTemp;
  std::string forPrimer;
  std::string forSeq;
   uint32_t revPos;
  double revTemp;
  std::string revPrimer;
  std::string revSeq;
};

struct sortProducts
{
    inline bool operator() (const pcrProduct& a, const pcrProduct& b)
    {
        return (a.penalty < b.penalty);
    }
};

void addUniqe(std::vector<primerBind>* coll, primerBind* prim) {
  int found = 0;
  for(std::vector<primerBind>::iterator it = coll->begin(); it != coll->end(); ++it) {
    if ((prim->chrom.compare(it->chrom) == 0) &&
        (prim->pos - it->pos == 0) && 
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

  // Initialize thal arguments
  boost::filesystem::path exepath = boost::filesystem::system_complete(argv[0]).parent_path();
  std::string cfgpath = exepath.string() + "/primer3_config/";
  primer3thal::thal_args a;
  primer3thal::set_thal_default_args(&a);
  a.temponly=1;
  a.type = primer3thal::thal_end1;
  primer3thal::get_thermodynamic_values(cfgpath.c_str());

  // CMD Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome file")
    ;

  boost::program_options::options_description appr("Approximate Search Options");
  appr.add_options()
    ("kmer,k", boost::program_options::value<uint32_t>(&c.kmer)->default_value(15), "k-mer size")
    ("maxmatches,m", boost::program_options::value<std::size_t>(&c.max_locations)->default_value(10000), "max. number of matches per k-mer")
    ("distance,d", boost::program_options::value<uint32_t>(&c.distance)->default_value(1), "neighborhood distance")
    ("hamming,n", "use hamming neighborhood instead of edit distance")
    ;

  boost::program_options::options_description score("Parameters for Scoring and Penalty Calculation");
  score.add_options()
    ("cutTemp,c", boost::program_options::value<double>(&c.cutTemp)->default_value(40.0), "minimal primer meting temperature to consider")
    ("maxProdSize,l", boost::program_options::value<uint32_t>(&c.maxProdSize)->default_value(15000), "maximal PCR Product size amplified")
    ("targetTm", boost::program_options::value<double>(&c.targetTemp)->default_value(58.0), "intended target Tm for primers")
    ("CutoffPenalty", boost::program_options::value<double>(&c.cutofPen)->default_value(-1.0), "maximal penalty for products to consider, -1 = keep all")
    ("penaltyTmDiff", boost::program_options::value<double>(&c.penDiff)->default_value(0.2), "multiplication factor for deviation of primer Tm penalty")
    ("penaltyTmMismatch", boost::program_options::value<double>(&c.penMis)->default_value(0.4), "multiplication factor for Tm pair difference penalty")
    ("penaltyLength", boost::program_options::value<double>(&c.penLen)->default_value(0.001), "multiplication factor for amplicon length penalty")
    ;

  boost::program_options::options_description tmcalc("Parameters for Tm Calculation");
  tmcalc.add_options()
    ("enttemp", boost::program_options::value<double>(&a.temp)->default_value(37.0), "temperature for entropie and entalpie calculation in Celsius")
    ("monovalent", boost::program_options::value<double>(&a.mv)->default_value(50.0), "concentration of monovalent ions in mMol")
    ("divalent", boost::program_options::value<double>(&a.dv)->default_value(1.5), "concentration of divalent ions in mMol")
    ("dna", boost::program_options::value<double>(&a.dna_conc)->default_value(50.0), "concentration of annealing(!) Oligos in nMol")
    ("dntp", boost::program_options::value<double>(&a.dntp)->default_value(0.6), "the sum  of all dNTPs in mMol")
    ;

  boost::program_options::options_description outp("Output Options");
  outp.add_options()
    ("prefix,p", boost::program_options::value<std::size_t>(&c.pre_context)->default_value(3), "prefix length")
    ("suffix,s", boost::program_options::value<std::size_t>(&c.post_context)->default_value(3), "suffix length")
    ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("amplicons.txt"), "output file")
    ("primer,p", boost::program_options::value<boost::filesystem::path>(&c.primfile)->default_value("primers.fa"), "primer locations file")
    ("align,x", "write alignments to stderr")
    ;
    
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "seq.fasta")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(appr).add(score).add(tmcalc).add(outp).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(appr).add(score).add(tmcalc).add(outp);
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

  // Fix provided Temperature
  a.temp += primer3thal::ABSOLUTE_ZERO;

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

  // Query FM-Index
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Query FM-Index" << std::endl;
  std::vector<primerBind> forBind;
  std::vector<primerBind> revBind;  
  for(TInputFasta::const_iterator itFa = infa.begin(); itFa != infa.end(); ++itFa) {
    std::string qr = itFa->second;
    if (qr.size() < c.kmer) continue;
    int32_t koffset = qr.size() - c.kmer;
    qr = qr.substr(qr.size() - c.kmer);
    typedef std::set<primLoci> TStringSet;
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
	std::string query(its->seq);
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
            prim.temp = o.temp;
            prim.primer = itFa->first;
            prim.leng = 0;
            prim.seq = itFa->second;
            prim.genSeq = genomicseq;
            if (o.temp > c.cutTemp) {
              if (fwdrev == 0) {
                prim.onFor = true;
                prim.pos = chrpos - its->fivePrim;
                //forBind.push_back(prim);
                addUniqe(&forBind, &prim);
              } else {
                prim.onFor = false;
                prim.pos = chrpos + its->threePrim;
                //revBind.push_back(prim);
                addUniqe(&revBind, &prim);
              }
            }
                           
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
        pcrProd.forSeq = fw->seq;
        pcrProd.revPos = rv->pos;
        pcrProd.revTemp = rv->temp;
        pcrProd.revPrimer = rv->primer;
        pcrProd.revSeq = rv->seq;
        pcrProd.seq = "stest";
        // Calculate Penalty
        double pen = std::abs(c.targetTemp - fw->temp) * c.penDiff;
        pen += std::abs(c.targetTemp - rv->temp) * c.penDiff;
        pen += std::abs(fw->temp - rv->temp) * c.penMis;
        pen += pcrProd.leng * c.penLen;

        pcrProd.penalty = pen;
        if ((c.cutofPen < 0) || (pen < c.cutofPen)) {
          pcrColl.push_back(pcrProd);
        }
      }
    }
  }

  std::sort(pcrColl.begin(), pcrColl.end(), sortProducts());

  std::ofstream rfile(c.outfile.string().c_str());
  int count = 0;
  for(std::vector<pcrProduct>::iterator it = pcrColl.begin(); it != pcrColl.end(); ++it) {
    rfile << "Amplicon_" << count << "_Length=" << it->leng << std::endl;
    rfile << "Amplicon_" << count << "_Penalty=" << it->penalty << std::endl;
    rfile << "Amplicon_" << count << "_For_Pos=" << it->chrom << ":" << it->forPos << std::endl;
    rfile << "Amplicon_" << count << "_For_Tm=" << it->forTemp << std::endl;
    rfile << "Amplicon_" << count << "_For_Name=" << it->forPrimer << std::endl;
    rfile << "Amplicon_" << count << "_For_Seq=" << it->forSeq << std::endl;
    rfile << "Amplicon_" << count << "_Rev_Pos=" << it->chrom << ":" << it->revPos << std::endl;
    rfile << "Amplicon_" << count << "_Rev_Tm=" << it->revTemp << std::endl;
    rfile << "Amplicon_" << count << "_Rev_Name=" << it->revPrimer << std::endl;
    rfile << "Amplicon_" << count << "_Rev_Seq=" << it->revSeq << std::endl;
    rfile << "Amplicon_" << count << "_Seq=" << it->seq<< std::endl;
    count++;
  }
  rfile.close();

  //Print Primers
  forBind.insert( forBind.end(), revBind.begin(), revBind.end() );
  std::sort(forBind.begin(), forBind.end(), sortPrimer());
  std::ofstream forfile(c.primfile.string().c_str());
  count = 0;
  for(std::vector<primerBind>::iterator it = forBind.begin(); it != forBind.end(); ++it) {
    forfile << "Primer_" << count << "_Tm="  << it->temp << std::endl;
    forfile << "Primer_" << count << "_Pos="  <<  it->chrom << ":" << it->pos << std::endl;
    forfile << "Primer_" << count << "_Ori=";
    if (it->onFor) forfile << "forward" << std::endl;
    else forfile << "reverse" << std::endl;
    forfile << "Primer_" << count << "_Name"  << it->primer << std::endl;
    forfile << "Primer_" << count << "_Seq="  << it->seq << std::endl;
    forfile << "Primer_" << count << "_Genome="  << it->genSeq << std::endl;
    count++;
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
