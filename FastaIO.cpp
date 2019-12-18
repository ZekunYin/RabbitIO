/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
  
#include <cstdio>
#include <vector>
#include <map>
#include <cstdio>

#include "FastaIO.h"
#include "Buffer.h"
#include "FastaStream.h" 
#include "Sequence.h"

namespace dsrc
{

namespace fq
{

FastaChunk* FastaReader::readNextChunk(){
	FastaDataChunk* part = NULL;
	recordsPool.Acquire(part);
	FastaChunk *dataPart = new FastaChunk;
	dataPart->chunk = part;
	if(fileReader.ReadNextChunk(dataPart, this->seqInfos))
	{
		return dataPart;
	}
	else
	{
		recordsPool.Release(part);
		return NULL;
	}
}

int chunkFormat(FastaDataChunk* &chunk, std::vector<Sequence*> &data, bool mHasQuality){
	//format a whole chunk and return number of reads
	int seq_count = 0;
	int line_count = 0;
	int pos_ = 0;

	while(true){
		//TODO rewrite to deal with part sequence
		//string name = getLine(chunk, pos_);
		//if(name.empty()) break;//dsrc guarantees that read are completed!
		//std::cout << name << std::endl;

		//string sequence = getSequence(chunk, pos_);
		//std::cout << sequence << std::endl;

		//data.push_back(new Read(name, sequence));
		//seq_count++;

	}

	return seq_count;
}

string getSequence(FastaDataChunk* &chunk, int &pos){	//addbyxxm
	int start_pos = pos;
	char * data = (char *)chunk->data.Pointer();

	while(pos <= (chunk->size + 1))
	{
		if(data[pos] == '\n' || data[pos] == '\r' || pos == (chunk->size + 1)){//the pos == chunk->size + 1 cannot be true when the final char is '\n'
			pos++;
			if(data[pos] == '>' || pos == chunk->size + 1)//so this is the pos == chunk->size + 1;
				return string(data+start_pos, pos-start_pos-1);
		}
		else{
			pos++;
		}
	}
	return "";
}

string getLine(FastaDataChunk* &chunk, int &pos){
	int start_pos = pos;
	char* data = (char *)chunk->data.Pointer();

	while(pos <= (chunk->size + 1)){
		if(data[pos] == '\n' || data[pos] == '\r' || pos == (chunk->size + 1)){
			//find a line
			pos++;
			return string(data+start_pos, pos-start_pos - 1);
		}
		else{
			pos++;
		}
	}
	return "";
}

} // namesapce fq

} // namespace dsrc
