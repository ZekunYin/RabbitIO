#ifndef H_SEQUENCE
#define H_SEQUENCE

#include <string>

namespace dsrc
{
namespace fq
{

using namespace std;

class Sequence{
public:
	Sequence(bool complete = true):mComplete(complete){};
	Sequence(string name, string seq, bool complete = true):
			mName(name),
			mSeq(seq),
			mComplete(complete){};
	Sequence(string seq, bool complete = true):
			mSeq(seq), mComplete(complete){};
	
	~Sequence(){};

public:
	string mName;
	string mSeq;
	bool mComplete;	
};

} //dsrc
} //fq
#endif // H_SEQUENCE
