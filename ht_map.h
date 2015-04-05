#include <map>
#include <string>
#include <set>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <algorithm>

using namespace std;

/***
	Creates hash map for reads
	1 - AGTC...C
	2 - ATCA...A
	3 - GGGC...C
	...
	
	Selecting
fw	5'-AGCTCAGAACAAGTTTACGAGGCAGTTACC...ATTACAGGACTGGATGGACTTAGGT-3'
    key1                                ATTACAGGAC
	key2								     AGGACTGGAT
	key3									      TGGATGGACT
	key4										       GGACTTAGGT

rc	5'-ACCTAAGTCCATCCAGTCCTGTAAT...GGTAACTGCCTCGTAAACTTGTTCTGAGCT-3'
	key-1								CTGCCTCGTA
	key-2								     TCGTAAACTT
	key-3									      AACTTGTTCT
	key-4										       GTTCTGAGCT
													   
	1. Hash rightmost 35 nt of fw and of rc (total 7+7 hash entries for 1 read)
		This version: hash 30 nt starting at end-31 of each read (1 hash from 1 read)
	
***/

string IntToStr(int tmp);

vector<string> &Split(const string &s, char delim, vector<string> &elems);
vector<string> Split(const string &s, char delim);

class HRD
{
public:
	int n;
	int offset;
};
string RCDNA(const string& DNA);
void CreateHashMapReadsFromFASTQ(const char* fqfile, map<int,string>& reads, map<string, vector<HRD> >& hash_fw, map<string, vector<HRD> >& hash_rc);
void LoadReadsFromFASTQ(const char* fqfile, map<int,string>& reads);
/***
	Having fw and rc hashmaps, map reads on the reference
***/

void ReadWindowFile(const char* fafile,string& ref, string& hap0, int& _ap, string& mainpath);
void ReadReferenceFromFile(const char* fafile,string& ref);
void MapReadsOnReference(string& ref, map<int,string>& reads, map<string, vector<HRD> >& hash_fw, map<string, vector<HRD> >& hash_rc, map<int,int>& fw_map, map<int,int>& rc_map);
void MapReadsOnReferenceMultiSeed(string& ref,int ap, map<int,string>& reads, map<int,int>& fw_map, map<int,int>& rc_map);
void MapReadsOnReferenceMultiSeedUltra(string& ref,int ap, map<int,string>& reads, map<int,int>& fw_map, map<int,int>& rc_map);


/***
	Function builds T-profile of read-ref alignment:
	-----GTTTGGATTGAGTGAAG--TGGATTTGGGCCC
	-----G---GGA--GAGTGAAG--TGGATTTGGGCCC
	---TGGTTTGGATTGAGTGAAG--TGGATTTGGGCCC
	ACATGG---GGA--GAG-GAAGTTTGGA---GGGCCC
	
	Will have following profiles:
	1TTT3TT3T4T3TTT6
	7T4T3TTT6
	T2TTT3TT3T4T3TTT6
	3T12TTT9
***/
class JPair
{
public:
	int x;
	int y;
};

bool SortByY(JPair a, JPair b);

void ClassifyEditingSitesOnGivenPath(string& ref, string& aln_ref, string& path, map<int,int>& sites);
bool CompareReadWithPath(string& read, string& path, string& ref, map<int, int>& sites, ofstream& ofile);
void CompareReadWithMainPath(string& read, string& path, string& ref, map<int, int>& sites, ofstream& ofile);

void AlignReadsWithRef(int alnpos, set<int>& mapped, map<int,int>& fw, map<int,int>& rc, map<int,string>& rd, string ref, map<int,string>& aln);

void TruncateRead(int alnpos, string& _rd, string& _ref);