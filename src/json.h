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

#ifndef JSON_H
#define JSON_H

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

namespace silica {

template<typename TPcrProducts, typename TPrimerName, typename TPrimerSeq>
inline void
ampliconTxtOut(std::string const& outfile, faidx_t* fai, TPcrProducts const& pcrColl, TPrimerName const& pName, TPrimerSeq const& pSeq) {
  std::ofstream rfile(outfile.c_str());
  int32_t count = 0;
  for(typename TPcrProducts::const_iterator it = pcrColl.begin(); it != pcrColl.end(); ++it, ++count) {
    std::string chrom(faidx_iseq(fai, it->refIndex));
    rfile << "Amplicon_" << count << "_Length=" << it->leng << std::endl;
    rfile << "Amplicon_" << count << "_Penalty=" << it->penalty << std::endl;
    rfile << "Amplicon_" << count << "_For_Pos=" << chrom << ":" << it->forPos << std::endl;
    rfile << "Amplicon_" << count << "_For_Tm=" << it->forTemp << std::endl;
    rfile << "Amplicon_" << count << "_For_Name=" << pName[it->forId] << std::endl;
    rfile << "Amplicon_" << count << "_For_Seq=" << pSeq[it->forId] << std::endl;
    rfile << "Amplicon_" << count << "_Rev_Pos=" << chrom << ":" << it->revPos << std::endl;
    rfile << "Amplicon_" << count << "_Rev_Tm=" << it->revTemp << std::endl;
    rfile << "Amplicon_" << count << "_Rev_Name=" << pName[it->revId] << std::endl;
    rfile << "Amplicon_" << count << "_Rev_Seq=" << pSeq[it->revId] << std::endl;
    int32_t sl = -1;
    char* seq = faidx_fetch_seq(fai, chrom.c_str(), it->forPos, it->revPos, &sl);
    std::string seqstr = boost::to_upper_copy(std::string(seq));
    rfile << "Amplicon_" << count << "_Seq=" << seqstr << std::endl;
    free(seq);
  }
  rfile.close();
}

template<typename TPcrProducts, typename TPrimerName, typename TPrimerSeq>
inline void
ampliconJsonOut(std::string const& outfile, faidx_t* fai, TPcrProducts const& pcrColl, TPrimerName const& pName, TPrimerSeq const& pSeq) {
  std::ofstream rfile(outfile.c_str());
  int32_t count = 0;
  rfile << "[" << std::endl;
  for(typename TPcrProducts::const_iterator it = pcrColl.begin(); it != pcrColl.end(); ++it, ++count) {
    std::string chrom(faidx_iseq(fai, it->refIndex));
    if (count) rfile << "," << std::endl;
    rfile << "{\"Id\": " << count << ", ";
    rfile << "\"Length\": " << it->leng << ", ";
    rfile << "\"Penalty\": " << it->penalty << ", ";
    rfile << "\"Chrom\": \"" << chrom << "\", ";
    rfile << "\"ForPos\": " << it->forPos << ", ";
    rfile << "\"ForTm\": " << it->forTemp << ", ";
    rfile << "\"ForName\": \"" << pName[it->forId] << "\", ";
    rfile << "\"ForSeq\": \"" << pSeq[it->forId] << "\", ";
    rfile << "\"RevPos\": " << it->revPos << ", ";
    rfile << "\"RevTm\": " << it->revTemp << ", ";
    rfile << "\"RevName\": \"" << pName[it->revId] << "\", ";
    rfile << "\"RevSeq\": \"" << pSeq[it->revId] << "\", ";
    int32_t sl = -1;
    char* seq = faidx_fetch_seq(fai, chrom.c_str(), it->forPos, it->revPos, &sl);
    std::string seqstr = boost::to_upper_copy(std::string(seq));
    rfile << "\"Seq\": \"" << seqstr << "\"}";
    free(seq);
  }
  rfile << std::endl;
  rfile << "]" << std::endl;
  rfile.close();
}

}

#endif
