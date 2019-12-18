/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTA_CHUNK
#define H_FASTA_CHUNK

#include "Globals.h"
#include "Common.h"
#include "Buffer.h"
#include "utils.h"

#include <vector>

namespace dsrc
{

namespace fq
{

typedef core::DataChunk FastaDataChunk;

struct FastaChunk{

	FastaDataChunk * chunk;
	uint64 start;
	uint64 end;
	uint64 nseqs;
	bool startSplit;
	bool endSplit;

};

template<class SketchRes>
struct SeqInfo{
	
	uint64 gid; //sequence global id		
	std::vector<SketchRes> sketchs; 
	bool is_complete;

};

//TODO: replace it with real sketch library
class MySketch {
public:
	MySketch(){};
	~MySketch(){};
	void merge(){};

public:
	std::vector<uint64> sketch;
	
};

typedef std::vector<SeqInfo<MySketch> > SeqInfos;

} // namespace fq

} // namespace dsrc

#endif
