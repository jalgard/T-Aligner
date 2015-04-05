#include "sw_t.h"
#include <iostream>
using namespace std;

extern int mismatch_max_g;

string DNAupperCase(const string& s)
{
	string us="";
	for(int i=0;i<s.size();i++)
	{
		if(s[i]=='a') us+="A";
		else if(s[i]=='t') us+="T";
		else if(s[i]=='g') us+="G";
		else if(s[i]=='c') us+="C";
		else us+=s[i];
	}
	return us;
}

void SWFillAlnMatrix(int** H,const string& Rd,const string& Rf)
{
	for(int i=0; i <= Rd.length(); i++)	
		H[i][0]=0;
	for(int j=0; j <= Rf.length(); j++)
		H[0][j]=0;
	for(int i=1; i <= Rd.length(); i++)
	for(int j=1; j <= Rf.length(); j++)
	{
		int diag=0;
		Rd[i-1]==Rf[j-1]?diag=H[i-1][j-1]+Smatch:diag=H[i-1][j-1]+Smis;
		int Agap=Rf[j-1]=='T'?H[i][j-1]+SgapT:H[i][j-1]+Sgap;
		int Bgap=Rd[i-1]=='T'?H[i-1][j]+SgapT:H[i-1][j]+Sgap;
		int mS=0;
		Agap>Bgap?mS=Agap:mS=Bgap;	
		mS>diag?mS=mS:mS=diag;
		H[i][j]=mS;
	}
}

void SWBacktrack(int** H,const string& Rd,const string& Rf, string& Rda, string& Rfa)
{
	int i = Rd.length()-1;	
	int j = Rf.length()-1;	
	while(1)
	{
		if(i>=0&&j>=0&&((Rd[i]==Rf[j] && H[i][j]+Smatch==H[i+1][j+1])||(Rd[i]!=Rf[j] && H[i][j]+Smis==H[i+1][j+1])))
		{
			Rda=Rd[i]+Rda;
			Rfa=Rf[j]+Rfa;
			i--;j--;
		}
		else if(j>=0&&(((H[i+1][j]+SgapT==H[i+1][j+1])&&Rf[j]=='T')||(H[i+1][j]+Sgap==H[i+1][j+1])))
		{
			Rda="-"+Rda;
			Rfa=Rf[j]+Rfa;
			j--;
		}
		else if(i>=0&&(((H[i][j+1]+SgapT==H[i+1][j+1])&&Rd[i]=='T')||(H[i][j+1]+Sgap==H[i+1][j+1])))
		{
			Rda=Rd[i]+Rda;
			Rfa="-"+Rfa;
			i--;
		}
		if(i==0||j==0)	
		{	
			if(Rd[i]==Rf[j]) 
			{
				Rda=Rd[i]+Rda; Rfa=Rf[j]+Rfa; i--; j--;
			}
			while(i>=0) {Rda=Rd[i]+Rda; Rfa="-"+Rfa;i--;}	
			while(j>=0) {Rda="-"+Rda; Rfa=Rf[j]+Rfa;j--;}
			break;
		}
	}
}

bool SWBacktrackT(int** H,const string& Rd,const string& Rf, string& Rda, string& Rfa)
{
	int i = Rd.length()-1;	
	int j = Rf.length()-1;	
	while(1)
	{
		if(i>=0&&j>=0&&(Rd[i]==Rf[j] && H[i][j]+Smatch==H[i+1][j+1]))
		{
			Rda=Rd[i]+Rda;
			Rfa=Rf[j]+Rfa;
			i--;j--;
		}
		else if(i>=0&&j>=0&&(Rd[i]!=Rf[j] && H[i][j]+Smis==H[i+1][j+1]))
		{
			return false;
		}
		else if(j>=0&&H[i+1][j]+SgapT==H[i+1][j+1]&&Rf[j]=='T')
		{
			Rda="-"+Rda;
			Rfa=Rf[j]+Rfa;
			j--;
		}
		else if(j>=0&&(H[i+1][j]+Sgap==H[i+1][j+1]))
		{
			return false;
		}
		else if(i>=0&&H[i][j+1]+SgapT==H[i+1][j+1]&&Rd[i]=='T')
		{
			Rda=Rd[i]+Rda;
			Rfa="-"+Rfa;
			i--;
		}
		else if(i>=0&&(H[i][j+1]+Sgap==H[i+1][j+1]))
		{
			return false;
		}
		if(i==0||j==0)	
		{	
			if(j==0&&i>0) return false;
			if(Rd[i]==Rf[j]) 
			{
				Rda=Rd[i]+Rda; Rfa=Rf[j]+Rfa; i--; j--;
			}
			while(i>=0) {Rda=Rd[i]+Rda; Rfa="-"+Rfa;i--;}	
			while(j>=0) {Rda="-"+Rda; Rfa=Rf[j]+Rfa;j--;}
			return true;
		}
	}
	return true;
}

bool SWBacktrackT2(int** H,const string& Rd,const string& Rf, string& Rda, string& Rfa)
{
	int i = Rd.length()-1;	
	int j = Rf.length()-1;	
	int mmc=0;
	while(1)
	{
		if(i>=0&&j>=0&&(Rd[i]==Rf[j] && H[i][j]+Smatch==H[i+1][j+1]))
		{
			Rda=Rd[i]+Rda;
			Rfa=Rf[j]+Rfa;
			i--;j--;
		}
		else if(i>=0&&j>=0&&(Rd[i]!=Rf[j] && H[i][j]+Smis==H[i+1][j+1]))
		{
			if(mmc>=mismatch_max_g||Rf[i]=='T'||Rd[j]=='T') return false;
			mmc++;
			Rda=Rf[i]+Rda; //correct letter
			Rfa=Rf[j]+Rfa;
			i--;j--;
		}
		else if(j>=0&&H[i+1][j]+SgapT==H[i+1][j+1]&&Rf[j]=='T')
		{
			Rda="-"+Rda;
			Rfa=Rf[j]+Rfa;
			j--;
		}
		else if(j>=0&&(H[i+1][j]+Sgap==H[i+1][j+1]))
		{
			return false;
		}
		else if(i>=0&&H[i][j+1]+SgapT==H[i+1][j+1]&&Rd[i]=='T')
		{
			Rda=Rd[i]+Rda;
			Rfa="-"+Rfa;
			i--;
		}
		else if(i>=0&&(H[i][j+1]+Sgap==H[i+1][j+1]))
		{
			return false;
		}
		if(i==0||j==0)	
		{	
			if(j==0&&i>0) return false;
			if(Rd[i]==Rf[j]) 
			{
				Rda=Rd[i]+Rda; Rfa=Rf[j]+Rfa; i--; j--;
			}
			while(i>=0) {Rda=Rd[i]+Rda; Rfa="-"+Rfa;i--;}	
			while(j>=0) {Rda="-"+Rda; Rfa=Rf[j]+Rfa;j--;}
			return true;
		}
	}
	return true;
}
/* VER.OLD.
bool SWBacktrackT(int** H,const string& Rd,const string& Rf, string& Rda, string& Rfa)
{
	int i = Rd.length()-1;	
	int j = Rf.length()-1;	
	while(1)
	{
		if(i>=0&&j>=0&&(Rd[i]==Rf[j] && H[i][j]+Smatch==H[i+1][j+1]))
		{
			Rda=Rd[i]+Rda;
			Rfa=Rf[j]+Rfa;
			i--;j--;
		}
		else if(i>=0&&j>=0&&(Rd[i]!=Rf[j] && H[i][j]+Smis==H[i+1][j+1]))
		{
			return false;
		}
		else if(j>=0&&H[i+1][j]+SgapT==H[i+1][j+1]&&Rf[j]=='T')
		{
			Rda="-"+Rda;
			Rfa=Rf[j]+Rfa;
			j--;
		}
		else if(j>=0&&(H[i+1][j]+Sgap==H[i+1][j+1]))
		{
			return false;
		}
		else if(i>=0&&H[i][j+1]+SgapT==H[i+1][j+1]&&Rd[i]=='T')
		{
			Rda=Rd[i]+Rda;
			Rfa="-"+Rfa;
			i--;
		}
		else if(i>=0&&(H[i][j+1]+Sgap==H[i+1][j+1]))
		{
			return false;
		}
		if(i==0||j==0)	
		{	
			if(j==0&&i>0) return false;
			if(Rd[i]==Rf[j]) 
			{
				Rda=Rd[i]+Rda; Rfa=Rf[j]+Rfa; i--; j--;
			}
			while(i>=0) {Rda=Rd[i]+Rda; Rfa="-"+Rfa;i--;}	
			while(j>=0) {Rda="-"+Rda; Rfa=Rf[j]+Rfa;j--;}
			return true;
		}
	}
	return true;
}
*/


void AlignWaterReadnRefRegion(int** H, string& Rd, string& Rf, int RdFrom, int RfFrom)
{
	string RdPart = Rd.substr(0,RdFrom);
	string RfPart = Rf.substr(0,RfFrom);
	string alignedRdPart="";
	string alignedRfPart="";
	SWFillAlnMatrix(H,RdPart,RfPart);
	SWBacktrack(H,RdPart,RfPart,alignedRdPart,alignedRfPart);
	Rd=alignedRdPart+Rd.substr(RdFrom);
	Rf=alignedRfPart+Rf.substr(RfFrom);
}

bool AlignWaterReadnRefRegionT(int** H, string& Rd, string& Rf, int RdFrom, int RfFrom)
{
	string RdPart = Rd.substr(0,RdFrom);
	string RfPart = Rf.substr(0,RfFrom);
	string alignedRdPart="";
	string alignedRfPart="";
	SWFillAlnMatrix(H,RdPart,RfPart);
	if(SWBacktrackT(H,RdPart,RfPart,alignedRdPart,alignedRfPart))
	{
		Rd=alignedRdPart+Rd.substr(RdFrom);
		Rf=alignedRfPart+Rf.substr(RfFrom);
		return true;
	}
	else
	{
		return false;
	}
}

bool AlignWaterReadnRefRegionT2(int** H, string& Rd, string& Rf, int RdFrom, int RfFrom)
{
	string RdPart = Rd.substr(0,RdFrom);
	string RfPart = Rf.substr(0,RfFrom);
	string alignedRdPart="";
	string alignedRfPart="";
	SWFillAlnMatrix(H,RdPart,RfPart);
	if(SWBacktrackT2(H,RdPart,RfPart,alignedRdPart,alignedRfPart))
	{
		Rd=alignedRdPart+Rd.substr(RdFrom);
		Rf=alignedRfPart+Rf.substr(RfFrom);
		return true;
	}
	else
	{
		return false;
	}
}
/***
int main(int argc, char** argv)
{
	string Ref="AAGTATTATAACTTACAATTTGTTAACGTCAGTTTTTTTTATGTGTAAAGCTAGCATTAATTATGCTTCGACTAAGCTTACATACAAAATGGCGTTAACCCCGGAAAGAAGGGACAAGACCAAGAAGACCAGAGAAGCGGGAACAGCGATGGTGGGAGAACGGAGCAGTCGTAGGAAAGCCCACAGCCAGAGCTGTTTTGGTCGTAGAGGGGTTCGTTTTTCCGAGAATAAAGATTCGTTTCTGGAAGGGGAGTCAGGCCTACTGATTTTTTCGCCAACTTTTACACAGGAGAGGGCCAGTTTAAGTTTGCACCCAAAAAGGGGAGATTTTTCGAAGGGAAAATTAATTTATTAATTTGTTAGATAATATAAACGAAAACGTCACGCTAAAAATCC";
	string Read="ATTTTTTTGTTTGGCCATGTATTTATGTGCACCCAAATTTATTTTATTTGGGGAGATTTTTCGAAGGGAAAATTAATTTATTAATTTGTTAGATAATATAAACGAAAACGTCACGCTAAAAATCC";	
	int** H = new int*[10000];	
	for(int i=0;i<10000;i++)
		H[i]=new int[10000];
	cout << AlignWaterReadnRefRegionT(H,Read,Ref,Read.size(),Ref.size()) << endl;
	cout << Read << endl;
	cout << Ref << endl;
}
***/