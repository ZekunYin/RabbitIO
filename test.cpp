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
	string fileName = "test.fna";
	dsrc::fq::FastaDataPool *fastaPool = new dsrc::fq::FastaDataPool(128, 1 << 26);
	dsrc::fq::FastaFileReader *fileReader = new dsrc::fq::FastaFileReader(fileName);
	dsrc::fq::FastaReader *fastaReader = new dsrc::fq::FastaReader(*fileReader, *fastaPool);
	dsrc::fq::FastaChunk *chunk;

	int nChunks = 0;
	while((chunk = fastaReader->readNextChunk()) != NULL)
	{
		nChunks++;	
		printf("datapart->%p\n", chunk);	
		printf("part->%p\n", chunk->chunk);	
	}
	cout << "nChunks: " << nChunks << endl;
	return 0;
}
