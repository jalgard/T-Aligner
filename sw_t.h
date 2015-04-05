
#include <string>
#include <vector>
#include <fstream>
using namespace std;

const int Smatch=15;	//15 //25
const int Smis=-14;	//-2 //-14
const int Sgap=-17;	//-4 //-17
const int SgapT=5;	//5 //5

string DNAupperCase(const string& s);

void SWFillAlnMatrix(int** H,const string& Rd,const string& Rf);
void SWBacktrack(int** H,const string& Rd,const string& Rf, string& Rda, string& Rfa);
bool SWBacktrackT(int** H,const string& Rd,const string& Rf, string& Rda, string& Rfa);
bool SWBacktrackT2(int** H,const string& Rd,const string& Rf, string& Rda, string& Rfa);

void AlignWaterReadnRefRegion(int** H,string& Rd, string& Rf, int RdFrom, int RfFrom);
bool AlignWaterReadnRefRegionT(int** H,string& Rd, string& Rf, int RdFrom, int RfFrom);
bool AlignWaterReadnRefRegionT2(int** H,string& Rd, string& Rf, int RdFrom, int RfFrom);