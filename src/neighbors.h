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

#ifndef NEIGHBORS_H
#define NEIGHBORS_H

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

namespace fmsearch {

class primLoci {
  public:
  std::string seq;
  int fivePrim;
  int threePrim;
  bool operator <(const primLoci &b) const {return this->seq < b.seq;};
  int area() {return 1;}
};

inline void
reverseComplement(std::string& sequence) {
  std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
  std::size_t i = 0;
  for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
    switch (*revIt) {
    case 'A': sequence[i]='T'; break;
    case 'C': sequence[i]='G'; break;
    case 'G': sequence[i]='C'; break;
    case 'T': sequence[i]='A'; break;
    case 'N': sequence[i]='N'; break;
    default: break;
    }
  }
}

// Special insert function that makes sure strset has no superstrings
template<typename TStrSet>
inline void
_insert(TStrSet& strset, primLoci const& s) {
  bool insertS = true;
  for (typename TStrSet::iterator it = strset.begin(); it != strset.end(); ) {
    if (it->seq.find(s.seq) != std::string::npos) strset.erase(it++); // s is a substring of *it, erase *it
    else {
      if (s.seq.find(it->seq) != std::string::npos) insertS = false; // *it is a substring of s, do not insert s
      ++it;
    }
  }
  if (insertS) strset.insert(s);
}

template<typename TAlphabet, typename TStringSet>
inline void
_neighbors(primLoci const& query, TAlphabet const& alphabet, int32_t dist, bool indel, int32_t pos, TStringSet& strset, int lb) {
  for(int32_t i = pos; i < (int32_t) query.seq.size();++i) {
    if (dist > 0) {
      if (indel) {
	// Insertion
	for(typename TAlphabet::const_iterator ait = alphabet.begin(); ait != alphabet.end(); ++ait) {
	  std::string ins("N");
	  ins[0] = *ait;
          primLoci newst;
	  newst.seq = query.seq.substr(0, i) + ins + query.seq.substr(i);
          if ((lb == 1) && (i == 0)) {
            newst.fivePrim = -1;
          } else {
            newst.fivePrim = 0;
          }
          if ((lb == 1) && (i == query.seq.size() - 1)) {
            newst.threePrim = -1;
          } else {
            newst.threePrim = 0;
          }
	  _neighbors(newst, alphabet, dist - 1, indel, pos, strset, 0);
	}
	// Deletion
	primLoci newst;
        newst.seq = query.seq.substr(0, i) + query.seq.substr(i + 1);
        if ((lb == 1) && (i == 0)) {
          newst.fivePrim = 1;
        } else {
          newst.fivePrim = 0;
        }
        if ((lb == 1) && (i == query.seq.size() - 1)) {
          newst.threePrim = 1;
        } else {
          newst.threePrim = 0;
        }
        _neighbors(newst, alphabet, dist - 1, indel, pos + 1, strset, 0);
      }
      for(typename TAlphabet::const_iterator ait = alphabet.begin(); ait != alphabet.end(); ++ait) {
	if (*ait != query.seq[i]) {
	  primLoci newst;
          newst.seq = query.seq;
	  newst.seq[i] = *ait;
          newst.fivePrim = 0;
          newst.threePrim = 0;
	  _neighbors(newst, alphabet, dist - 1, indel, pos+1, strset, 0);
	}
      }
    }
  }
  if ((indel) && (dist > 0)) {
    for(typename TAlphabet::const_iterator ait = alphabet.begin(); ait != alphabet.end(); ++ait) {
      std::string ins("N");
      ins[0] = *ait;
      primLoci newst;
      newst.seq = query.seq + ins;
      //strset.insert(newst);
      _insert(strset, newst);
    }
  }
  //strset.insert(query);
  _insert(strset, query);
}
      

template<typename TAlphabet, typename TStringSet>
inline void
neighbors(std::string const& query, TAlphabet const& alphabet, int32_t dist, bool indel, TStringSet& strset) {
  primLoci loc;
  loc.seq = query;
  _neighbors(loc, alphabet, dist, indel, 0, strset, 1);
}


}

#endif
