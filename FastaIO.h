/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/

#ifndef H_FASTQREADER
#define H_FASTQREADER

#include <vector>
#include <string>

#include "Globals.h"
#include "Common.h"
#include "DataQueue.h"
#include "DataPool.h"
#include "FastaStream.h"
#include "FastaChunk.h"
#include "Sequence.h"

namespace dsrc
{

namespace fq
{

typedef core::TDataQueue<FastaDataChunk> FastaDataQueue;
typedef core::TDataPool<FastaDataChunk> FastaDataPool;


class FastaReader //: public IFastqIoOperator
{
public:
    FastaReader(FastaFileReader& reader_, FastaDataPool& pool_)
	    :   recordsPool(pool_)
		,	fileReader(reader_)
		,	numParts(0)
	{};

	FastaChunk* readNextChunk();
	
	int64 Read(byte* memory_, uint64 size_)
	{
		int64 n = fileReader.Read(memory_, size_);
		return n;
	}

private:

	FastaDataPool&      recordsPool;
	FastaFileReader&	fileReader;
	uint32 numParts;
};

int chunkFormat(FastaDataChunk* &chunk, std::vector<Sequence*>&, bool);


string getSequence(FastaDataChunk* &chunk, int &pos);	//addbyxxm
string getLine(FastaDataChunk* &chunk, int &pos);

} // namespace fq

} // namespace dsrc

#endif
