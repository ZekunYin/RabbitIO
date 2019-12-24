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

	mash::fa::FastaDataPool *fastaPool    = new mash::fa::FastaDataPool(16, 1<<26);
	mash::fa::FastaFileReader *fileReader = new mash::fa::FastaFileReader(fileName);
	mash::fa::FastaReader *fastaReader    = new mash::fa::FastaReader(*fileReader, *fastaPool);
	mash::fa::FastaChunk *chunk;

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
