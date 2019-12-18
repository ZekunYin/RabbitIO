#include "FastaIO.h"
#include "FastaStream.h"
#include "FastaChunk.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdio>

using namespace std;

int main(){

	//init datapool
	string fileName = "./data/clean_total_genomic.fna";
	//string fileName = "./data/clean_fungi.58.1.genomic.fna";
	//string fileName = "test.fna";
	
	//vector<dsrc::fq::SeqInfo<MySketch> > seqInfos;

	dsrc::fq::FastaDataPool *fastaPool = new dsrc::fq::FastaDataPool(16, 1<<26);
	dsrc::fq::FastaFileReader *fileReader = new dsrc::fq::FastaFileReader(fileName);
	dsrc::fq::FastaReader *fastaReader = new dsrc::fq::FastaReader(*fileReader, *fastaPool);
	dsrc::fq::FastaChunk *chunk;

	int nChunks = 0;
	while((chunk = fastaReader->readNextChunk()) != NULL)
	{
		nChunks++;	
		chunk->print();
		fastaPool->Release(chunk->chunk);
		delete chunk;
	}
	cout << "nChunks: " << nChunks << endl;
	cout << "totalSeqs: " << fileReader->totalSeqs << endl;
	cout << "totalSeqs info: " << fastaReader->seqInfos.size() << endl;
	delete fastaReader;
	delete fileReader;	
	delete fastaPool;
	return 0;
}
