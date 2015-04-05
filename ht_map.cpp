#include "ht_map.h"
#include "sw_t.h"
#include <iostream>
#include <fstream>
#include <deque>
#define dbgSAY cout <<
#define dbgEND << endl;

using namespace std;

extern int seed_step_g;
extern int seed_len_g;

string IntToStr(int tmp)
{
    ostringstream out;
    out << tmp;
    return out.str();
}
vector<string> &Split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
vector<string> Split(const string &s, char delim) {
    vector<string> elems;
    Split(s, delim, elems);
    return elems;
}


string RCDNA(const string& DNA)
{
	string _rcdna=DNA;
	int _szdna=DNA.size()-1;
	map<char,char> _rcmap;
	_rcmap['A']='T';	_rcmap['G']='C';	_rcmap['C']='G';	_rcmap['T']='A';
	_rcmap['N']='N'; _rcmap['Y']='A'; _rcmap['W']='G'; _rcmap['H']='C'; _rcmap['R']='A'; // fix X
	_rcmap['M']='T'; _rcmap['S']='G'; _rcmap['D']='A'; _rcmap['B']='C'; _rcmap['K']='G'; // fix X
	_rcmap['V']='C'; // fix X
	for(int i=_szdna;i>=0;i--)
	{
		_rcdna[_szdna-i]=_rcmap[DNA[i]];
	}
	return _rcdna;
}
void CreateHashMapReadsFromFASTQ(const char* fqfile, map<int,string>& reads, map<string, vector<HRD> >& hash_fw, map<string, vector<HRD> >& hash_rc)
{
	ifstream fastq(fqfile);
	
	string line;
	unsigned int line_count=4;
	int read_count=0;
	
	while (!fastq.eof()) {
		getline(fastq, line, '\n' );
		if(line_count%4==1)
		{
			line_count++;
			if(line.size()>75)	// read must be at least 75 nt long
			{
				read_count++;
				reads[read_count]=line;
			}
		}
		else
		{
			line_count++;
		}
	}
	//dbgSAY "Reads loaded " << read_count dbgEND
	//dbgSAY "Map size: " << reads.size() dbgEND
	
	//dbgSAY "Hashing reads..." dbgEND
	for(map<int,string>::iterator rd=reads.begin(); rd!=reads.end(); rd++)
	{
		
		/***
			Let's take exactly 1 30-nt seed for exact match at #end-1 position of read (3')
		***/
		int _p = rd->second.size()-31;
		for(int _i=0;_i<1;_i++)
		{
			HRD _hrd;
			_hrd.n=rd->first;	_hrd.offset=_p;
			hash_fw[rd->second.substr(_p,30)].push_back(_hrd);
			_p-=30;
		}
		string _rc = RCDNA(rd->second);
		_p = _rc.size()-31;
		for(int _i=0;_i<1;_i++)
		{
			HRD _hrd;
			_hrd.n=rd->first;	_hrd.offset=_p;
			hash_rc[_rc.substr(_p,30)].push_back(_hrd);
			_p-=30;
		}
	}
	
	//dbgSAY "Forward read hashes " << hash_fw.size() dbgEND
	//dbgSAY "Reverse-complement read hashes " << hash_rc.size() dbgEND
	
}
void LoadReadsFromFASTQ(const char* fqfile, map<int,string>& reads)
{
	ifstream fastq(fqfile);
	
	string line;
	unsigned int line_count=4;
	int read_count=0;
	
	while (!fastq.eof()) {
		getline(fastq, line, '\n' );
		if(line[line.size()-1]=='\r') line=line.substr(0,line.size()-1);
		if(line_count%4==1)
		{
			line_count++;
			if(line.size()>75)	// read must be at least 75 nt long
			{
				read_count++;
				reads[read_count]=line;
			}
		}
		else
		{
			line_count++;
		}
	}
}
void ReadWindowFile(const char* fafile,string& ref, string& hap0, int& _ap, string& mainpath)
{
	ifstream fasta(fafile);
	string line;
	string fa_title="";
	map<string,string> refs;
	while (!fasta.eof()) {
		getline(fasta, line, '\n' );
		if(line[0]=='>')
		{
			if(line[line.size()-1]=='\r') line=line.substr(0,line.size()-1);
			refs[line.substr(1)]="";
			fa_title=line.substr(1);
		}
		else
		{
			if(line[line.size()-1]=='\r') line=line.substr(0,line.size()-1);
			refs[fa_title]+=line;
		}
	}
	ref=refs["ref"];
	hap0=refs["hap0"];
	mainpath=refs["main"];
	_ap=atoi(refs["start"].c_str());
}
void ReadReferenceFromFile(const char* fafile,string& ref)
{
	ifstream fasta(fafile);
	string line;
	string fa_title="";
	map<string,string> refs;
	while (!fasta.eof()) {
		getline(fasta, line, '\n' );
		if(line[0]=='>')
		{
			refs[line.substr(1)]="";
			fa_title=line.substr(1);
		}
		else
		{
			refs[fa_title]+=line;
		}
	}
	ref=refs[fa_title];
	//dbgSAY "Reference " << fa_title << " loaded: " << ref dbgEND
}

void MapReadsOnReference(string& ref, map<int,string>& reads, map<string, vector<HRD> >& hash_fw, map<string, vector<HRD> >& hash_rc, map<int,int>& fw_map, map<int,int>& rc_map)
{
	set<int> mapped;
	for(int _p=0; _p < ref.size()-30; _p++)
	{
		string _cpat = ref.substr(_p,30);
		if(hash_fw.count(_cpat))
		{
			for(vector<HRD>::iterator hh=hash_fw[_cpat].begin(); hh!=hash_fw[_cpat].end(); hh++)
			{
				if(mapped.count((*hh).n)==0)
				{
					mapped.insert((*hh).n);
					fw_map[(*hh).n]=_p;
				}
			}
		}
		if(hash_rc.count(_cpat))
		{
			for(vector<HRD>::iterator hh=hash_rc[_cpat].begin(); hh!=hash_rc[_cpat].end(); hh++)
			{
				if(mapped.count((*hh).n)==0)
				{
					mapped.insert((*hh).n);
					rc_map[(*hh).n]=_p;
				}
			}
		}
	}
	//dbgSAY "Total reads mapped: " << mapped.size() dbgEND
	//dbgSAY "In fw-orientation: " << fw_map.size() dbgEND
	//dbgSAY "In rc-orientation: " << rc_map.size() dbgEND
}

void MapReadsOnReferenceMultiSeed(string& ref,int ap, map<int,string>& reads, map<int,int>& fw_map, map<int,int>& rc_map)
{
	string seed_pat=ref.substr(ap+1,30);
	int _p=ap;
	for(map<int,string>::iterator it=reads.begin();it!=reads.end();it++)
	{
		//cout << "IN: " << it->second << endl << "RC: ";
		string _rc = RCDNA(it->second);
		//cout << _rc << endl << "Search: " << seed_pat << endl;
		size_t fwp = it->second.find(seed_pat);
		size_t rcp = _rc.find(seed_pat);
		//cout << seed_pat << ">" << "%" << fwp << ":" << rcp << " ";
		if(fwp!=string::npos&&fwp>50)
		{
			reads[it->first]=it->second.substr(0,fwp+30);
			fw_map[it->first]=_p;
		}
		if(rcp!=string::npos&&rcp>50)
		{
			string tmp = RCDNA(_rc.substr(0,rcp+30));
			reads[it->first]=tmp;
			rc_map[it->first]=_p;
		}
	}
	//dbgSAY "Total reads mapped: " << mapped.size() dbgEND
	dbgSAY "In fw-orientation: " << fw_map.size() dbgEND
	dbgSAY "In rc-orientation: " << rc_map.size() dbgEND
}

void MapReadsOnReferenceMultiSeedUltra(string& ref,int ap, map<int,string>& reads, map<int,int>& fw_map, map<int,int>& rc_map)
{
	set<int> mapped;
	
	for(int b=ap+1;b<=ref.size()-seed_len_g+1;b+=seed_step_g)
	{
		string seed_pat=ref.substr(b,seed_len_g);
		int _p=b-1;
		for(map<int,string>::iterator it=reads.begin();it!=reads.end();it++)
		{
			if(mapped.count(it->first)==0)
			{
				//cout << "IN: " << it->second << endl << "RC: ";
				string _rc = RCDNA(it->second);
				//cout << _rc << endl << "Search: " << seed_pat << endl;
				size_t fwp = it->second.find(seed_pat);
				size_t rcp = _rc.find(seed_pat);
				//cout << seed_pat << ">" << "%" << fwp << ":" << rcp << " ";
				if(fwp!=string::npos&&fwp>50)
				{
					mapped.insert(it->first);
					//cout << reads[it->first] << " was mapped at " << fwp << " with seed " << seed_pat << endl; 
					reads[it->first]=it->second.substr(0,fwp+seed_len_g);
					fw_map[it->first]=_p;
				}
				else if(rcp!=string::npos&&rcp>50) // may be else
				{
					mapped.insert(it->first);
					string tmp = RCDNA(_rc.substr(0,rcp+seed_len_g));
					reads[it->first]=tmp;
					rc_map[it->first]=_p;
				}
			}
		}
	}
	//dbgSAY "Total reads mapped: " << mapped.size() dbgEND
	dbgSAY "In fw-orientation: " << fw_map.size() dbgEND
	dbgSAY "In rc-orientation: " << rc_map.size() dbgEND
}

bool SortByY(JPair a, JPair b)
{
	return a.y > b.y;
}


/***
	Takes full length reference ref;
	reference and path aligned with main algorithm;
	returns map<int,int>:
	maps site position - number of non-T letters after left start of the reference
	to indel counter, which is = 0 when there is no editing;
	> 0 for any deletion, < 0 for insertion
***/
void ClassifyEditingSitesOnGivenPath(string& ref, string& aln_ref, string& path, map<int,int>& sites)
{
	int current_site=0;
	sites[0]=0;bool inited=false;
	for(int i=0;i<ref.size(); i++)
		if(ref[i]=='A'||ref[i]=='C'||ref[i]=='G')
		{	sites[++current_site]=0; }
	int num_sites=current_site; current_site=0;
	int indel=0;
	for(int i=0; i<aln_ref.size(); i++)
	{
		if(aln_ref[i]=='T'&&path[i]=='T') continue;
		if(aln_ref[i]=='-'&&path[i]=='-') continue;
		
		if(path[i]!='-'&&!inited) inited=true;
		
		if(inited&&aln_ref[i]=='T'&&path[i]=='-')	indel--;
		else if(inited&&path[i]=='T'&&aln_ref[i]=='-')	indel++;
		if(aln_ref[i]=='A'||aln_ref[i]=='G'||aln_ref[i]=='C')
		{
			sites[current_site++]=indel; indel=0;
		}
	}
	/***
	for(map<int,int>::iterator it=sites.begin(); it!=sites.end(); it++)
	{
		if(it->second>0) cout << "Del" << it->second << ", ";
		else if(it->second<0) cout << "Ins" << -1*it->second << ", ";
		else cout << "0, ";
	}
	***/
	
}
/***
	Given map<int,int> of classified sites,
	read and path aligned to the ref, returns true iff read do not 
	contradict the path
***/
bool CompareReadWithPath(string& read, string& path, string& ref, map<int, int>& sites, ofstream& ofile)
{
	int pos=0; int site=0; bool result=true; int tolerance=2;
	while(read[pos]=='-'&&pos<read.size()&&pos<path.size())
	{
		if(ref[pos]=='A'||ref[pos]=='G'||ref[pos]=='C') site++;
		pos++;
	}
	ostringstream offile;
	offile << "AB<" << pos << ">  ";
	int indel=0;
	while(pos<read.size()&&pos<path.size())
	{
		if(ref[pos]=='A'||ref[pos]=='G'||ref[pos]=='C') {
			if(result&&indel*sites[site]<0)
			{
			    //offile << "*Controversy in site " << site << endl;
			    result = false;
			}
			if(result&&sites[site]==0&&indel!=0)
			{
			    if(tolerance<2) tolerance++;
			    else {
				//offile << "*Path matched ref but read not at " << site << endl;
				result = false;
			    }
			}
			if(result&&abs(indel)>abs(sites[site]))
			{
			    if(tolerance<2) tolerance++;
			    else {
				//offile << "*Read has longer indel at " << site << endl;
				result = false;
			    }
			}
			if(indel>0) offile << "I" << indel << ", ";
			else if(indel<0) offile << "D" << -1*indel <<  ", ";
			else offile << "0, ";
			site++;indel=0;
		}
		if(read[pos]=='-'&&ref[pos]=='T') indel--;
		if(read[pos]=='T'&&ref[pos]=='-') indel++;
		pos++;
	}
	if(result) ofile << "\t&main&\t";
	ofile << offile.str() << endl;
	return result;
}

void CompareReadWithMainPath(string& read, string& path, string& ref, map<int, int>& sites, ofstream& ofile)
{
	int pos=0; int site=0;
	while(read[pos]=='-'&&pos<read.size()&&pos<path.size())
	{
		if(ref[pos]=='A'||ref[pos]=='G'||ref[pos]=='C') site++;
		pos++;
	}
	int indel=0;
	while(pos<read.size()&&pos<path.size())
	{
		if(ref[pos]=='A'||ref[pos]=='G'||ref[pos]=='C') {
			if(indel==sites[site]) ofile << "0, ";
			else if(indel*sites[site]<0) ofile << "!, ";
			else if(abs(indel)>abs(sites[site]))
			{
			    if(indel>0) ofile << "+"; else ofile << "-";
			    ofile << abs(indel)-abs(sites[site]) << ", ";
			}
			else if(sites[site]==0&&indel!=0) ofile << "|, ";
			else ofile << "0, ";
			site++;indel=0;
		}
		if(read[pos]=='-'&&ref[pos]=='T') indel--;
		if(read[pos]=='T'&&ref[pos]=='-') indel++;
		pos++;
	}
	ofile << endl;
}


void AlignReadsWithRef( int alnpos,  set<int>& mapped,  map<int,int>& fw,  map<int,int>& rc,  map<int,string>& rd,  string ref, map<int,string>& aln)
{
	map<int,int> gps;
	map<int,string> _alnseq;
	map<int,deque<char> > _rdsvec;
	for(set<int>::iterator it = mapped.begin(); it!=mapped.end(); it++)
	{
		string _rd = rd[*it];
		if(rc.count(*it)>0) _rd=RCDNA(_rd);
		while(_rd[0]=='T') _rd=_rd.substr(1);
		_rd=_rd.substr(0,_rd.size()-1);
		int _ntc=0;
		for(int i=0;i<_rd.size();i++)
		{
			_rdsvec[*it].push_back(_rd[i]);
			if(_rd[i]!='T') _ntc++; 
		}
		int _refoff=0;
		int _j=rc.count(*it)>0?rc[*it]:fw[*it];
		while(_ntc>0)
		{
			if(ref[_j+(seed_len_g-1)-_refoff]!='T') _ntc--; //29
			_refoff++;
		}
		gps[*it]=_j+(seed_len_g)-_refoff; //30
	}
	for(map<int, deque<char> >::iterator it=_rdsvec.begin(); it!=_rdsvec.end(); it++)
		for(int i=0;i<gps[it->first];i++)
			if(ref[i]!='T') it->second.push_front('-');
	
	for(int i=0; i<ref.size(); i++)
		_rdsvec[-2].push_back(ref[i]);
	
	int _rfsz=0;
	for(int i=0; i<ref.size(); i++)
		if(ref[i]!='T') _rfsz++;
	int _rftp=0;
	for(int i=0; i<_rfsz; i++)
	{
		bool indel=false;
		int _is=0;
		map<int,int> _Tmap;
		for(map<int, deque<char> >::iterator it=_rdsvec.begin(); it!=_rdsvec.end(); it++)
		{
			if(it->second.size()>0&&it->second.front()=='T')
			{
				int _t=1;
				it->second.pop_front();
				while(it->second.front()=='T') {_t++;it->second.pop_front();}
				indel=true;
				_Tmap[it->first]=_t;
				if(_t>_is) _is=_t;
			}
			else _Tmap[it->first]=0;
		}
		if(indel)
		{
			i--;
			for(int k=0;k<_is;k++)
				for(map<int, deque<char> >::iterator it=_rdsvec.begin(); it!=_rdsvec.end(); it++)
				{
					if(k<_Tmap[it->first]) _alnseq[it->first]+="T";
					else _alnseq[it->first]+="-";
				}
		}
		else
		{
			for(map<int, deque<char> >::iterator it=_rdsvec.begin(); it!=_rdsvec.end(); it++)
			{
				if(it->second.size()>0)
				{
					_alnseq[it->first]+=it->second.front();
					it->second.pop_front();
				}
				else _alnseq[it->first]+="-";
			}
		}
	}
	for(map<int,string>::iterator it=_alnseq.begin(); it!=_alnseq.end(); it++)
	{
		aln[it->first]=it->second; // 0 - reference;
	}
}

void TruncateRead(int alnpos, string& _rd, string& _ref)
{
	if(_rd.size()>alnpos+10)
	{
	cout << "Truncating read!\n";
		_rd=_rd.substr(_rd.size()-alnpos-5);
	}
}