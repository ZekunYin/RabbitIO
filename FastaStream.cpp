/*
  This file is a part of DSRC software distributed under GNU GPL 2 licence.
  The homepage of the DSRC project is http://sun.aei.polsl.pl/dsrc
  
  Authors: Lucas Roguski and Sebastian Deorowicz
  
  Version: 2.00
*/
 
#include "FastaStream.h"
#include <iostream>
namespace dsrc
{

namespace fq
{

bool FastaFileReader::ReadNextChunk(FastaChunk* dataChunk_)
{
	FastaDataChunk *chunk_ = dataChunk_->chunk;
	if (Eof())
	{
		chunk_->size = 0;
		return false;
	}

	// flush the data from previous incomplete chunk
	uchar* data = chunk_->data.Pointer();
	const uint64 cbufSize = chunk_->data.Size();
	chunk_->size = 0;
	int64 toRead = cbufSize - bufferSize;// buffersize: size left from last chunk
	
	if (bufferSize > 0)
	{
		std::copy(swapBuffer.Pointer(), swapBuffer.Pointer() + bufferSize, data);
		chunk_->size = bufferSize;
		bufferSize = 0;
	}

	// read the next chunk
	int64 r = this->Read(data + chunk_->size, toRead);
	//std::cout << "r is :" << r << std::endl;
	
	if (r > 0)
	{
		if (r == toRead)	// somewhere before end
		{
		    //uint64 chunkEnd = cbufSize - SwapBufferSize; // Swapbuffersize: 1 << 13
			////std::cout << "chunkend  cbufsize Swapbuffersize: " << chunkEnd <<" "<< cbufSize << " " << SwapBufferSize << std::endl;
			//chunkEnd = GetNextRecordPos(data, chunkEnd, cbufSize);
			//chunk_->size = chunkEnd - 1; //remove \n
			//if (usesCrlf)
			//	chunk_->size -= 1;

			//std::copy(data + chunkEnd, data + cbufSize, swapBuffer.Pointer());
			//bufferSize = cbufSize - chunkEnd;

			//dealing with halo region
			uint64 chunkEnd = findCutPos(data, cbufSize);
		}
		else				// at the end of file
		{
			chunk_->size += r - 1;	// skip the last EOF symbol
			if (usesCrlf)
				chunk_->size -= 1;

			eof = true;
		}
	}
	else
	{
		eof = true;
	}

	return true;
}

//需要向前找到倒数第2条read的开头，不仅仅是@号这么简单
uint64 FastaFileReader::GetPreviousRecordPos(uchar* data_, uint64 pos_,const uint64 size_)
{	int offset =2;

	SkipToSol(data_,pos_,size_);
	if(usesCrlf){
		offset=3;
	}
	while(data_[pos_+offset] !='@'){ //+2
		//std::cout<<"pos_"<<pos_<<std::endl;
		SkipToSol(data_,pos_,size_);
	}
	//标记一下，看是否是质量分
	uint64 pos0=pos_+offset;
	SkipToSol(data_,pos_,size_);
	if(data_[pos_+offset]== '+'){
		//说明上一个@号是质量分
		SkipToSol(data_,pos_,size_);
		SkipToSol(data_,pos_,size_);
		//此时应该是name
		if(data_[pos_+offset]!='@'){
			std::cout << "core dump is " << data_[pos_+offset] << std::endl;
		}else{
			return pos_+offset;
		}
	}
	else{
		return pos0;
	}
	return 0;
}

uint64 FastaFileReader::GetNextRecordPos(uchar* data_, uint64 pos_, const uint64 size_)
{
	SkipToEol(data_, pos_, size_);
	++pos_;

	if(1)//if(fasta) fastabyxxm
	{
		while(data_[pos_] != '>')
		{
			FastaSkipToEol(data_, pos_, size_);
			++pos_;
		}
		uint64 pos0 = pos_;
		return pos0;
	}
			
	else//if(fastq) 
	{
		// find beginning of the next record
		while (data_[pos_] != '@')
		{
			SkipToEol(data_, pos_, size_);
			++pos_;
		}
		uint64 pos0 = pos_;
			
		SkipToEol(data_, pos_, size_);
		++pos_;

		if (data_[pos_] == '@')			// previous one was a quality field
			return pos_;
		
		SkipToEol(data_, pos_, size_);
		++pos_;
		if(data_[pos_] != '+')
			std::cout << "core dump is pos: " << pos_ << " char: " << data_[pos_] << std::endl;	
		ASSERT(data_[pos_] == '+');	// pos0 was the start of tag
		return pos0;
	}
}

uint64 FastaFileReader::findCutPos(uchar* data_, const uint64 size_)
{
	int count = 0;
	uint64 pos_ = 0;
	uint64 cut_ = 0;
	if(data_[0] == '>') //start with '>'
	{
		while(pos_ < size_){
			if(data_[pos_] == '>'){
				if(FindEol(data_, pos_, size_)) //find name
				{
					++pos_;
				}else{
					cut_ = pos_; //find a cut: incomplete name
					break;
				}
			}else{
				++pos_;
			}
	
		}

	}else{ //start with {ACGT}
		while(pos_ < size_){
			if(data_[pos_] == '>'){
				if(FindEol(data_, pos_, size_)) //find name
				{
					++pos_;
				}else{
					cut_ = pos_; //find a cut: incomplete name
					break;
				}
			}else{
				++pos_;
			}
	
		}	
	}

	return cut_;	
}

} // namespace fq

} // namespace dsrc

