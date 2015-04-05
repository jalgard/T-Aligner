#include "ht_map.h"
#include "sw_t.h"
#include <iostream>
#include <set>
#include <fstream>
#include <sstream>

using namespace std;

/***
	Takes exactly 3 arguments:
	ref-file name [fasta];
	reads file [fastq];
	alignment-start coordinate [int];
	FASTA 
	>ref
	AGGATGACTGAGTAAAA-A---GTTTGGG
	>hapX
	               AATATTTGTTTGGG
	int = 16 
	
	ref_hybrid:
	AGGATGACTGAGTAAAATATTTGTTTGGG
***/

int seed_step_g=15;
int seed_len_g=25;
int mismatch_max_g=2;
int cut_val_g=100000;

int main(int argc, char** argv)
{
	//ofstream debug("debug.txt");
	
	if(argc < 6) 
	{
		cout << "In this version you also need to set 3 more <integer> options:\nseed length, seed step, max mismatch\n";
		cout << "seed length [10-40]\nseed step [10-40]\nmax mismatch [0-5]\nOptionally set alignment length cut...\n";
		return 1;
	}
	seed_len_g=atoi(argv[3]);
	seed_step_g=atoi(argv[4]);
	mismatch_max_g=atoi(argv[5]);
	if(argc > 6)
	{
		cut_val_g=atoi(argv[6]);
	}
	map<int,string> reads;
	//map<string, vector<HRD> > hash_fw;
	//map<string, vector<HRD> > hash_rc;
			
	string ref;
	string hap0;
	string mainp="";
	int ap;
	
	ReadWindowFile(argv[2], ref, hap0, ap, mainp);
	map<int,int> fw_map;
	map<int,int> rc_map;
	
	int** H = new int*[8000];	
	for(int i=0;i<8000;i++)
		H[i]=new int[8000];
	cout << "\nTAligner v2.5. update 09/07/14 /(c) Jalgard/\n";
	//cout << "\nChecking 'hap0' vs ref alignment...\n" << endl;
	string _ref=ref;
	string _hap0=hap0;
	//hap0="-"+hap0;
	/***
		Now we directly create hap0 hybrid without checking...
	***/	
	
	if(AlignWaterReadnRefRegionT(H,hap0,ref,hap0.size(),ref.size()))
	{
		cout << "Alignment of hap0 on reference:\n";
		cout << "length of alignment: " << hap0.size() << "/" << ref.size() << endl;
		cout << "you set alignment start pos " << ap << ", and I got: ";
		int _a=0;
		while(1) { if(hap0[_a]!='-') break; _a++; }
		cout << _a+1 << endl;
		if(_a+1==ap) cout << "\nYour start position verified! Trusted 'hap0'. Accepting.\n\n";
	}
	else	{
		cout << "\nError: can't find alignment of hap0 on ref!\n";
		cout << "Input file must look like this:\n\
		>ref\nAAAAAGGCCCCC\n\
		>hap0\nGTTGTTCCCCC\n>start\n6\n\n";
		return 1;
	}
	
	string hybref=_ref.substr(0,ap-1)+_hap0;
	//cout << "Hashing reads from fastq... [this can take a long time]" << endl;
	LoadReadsFromFASTQ(argv[1], reads);
	//CreateHashMapReadsFromFASTQ(argv[1], reads, hash_fw, hash_rc);
	//cout << "\nComplete!\n";
	if(mainp!="")
	{
		reads[-1]=mainp+_hap0;
	}
	
	cout << "\nSeeding reads on whole reference hybrid [ref|hap0]\n";
	
	MapReadsOnReferenceMultiSeedUltra(hybref,ap, reads, fw_map, rc_map);
	
	cout << "\nComplete!\nProceeding with Water-T [ideal T-alignment]..." << endl;
	
	if(ref.size()>3000) {
		cout << "\nError: bad allocation, not enough space for so long reference sequence!\nRecompile with large A constant!\n";
		return 1;
	}
	
	map<int,string> fw_aligned_rd;
	map<int,string> fw_aligned_rf;
	map<int,string> rc_aligned_rd;
	map<int,string> rc_aligned_rf;
	
	int _mapped_before_ap=0;
	
	for(map<int,int>::iterator it=fw_map.begin(); it!= fw_map.end(); it++)
	{
		if(it->second/*>==*/>=ap)
		{
			string _rd=reads[it->first].substr(0,reads[it->first].size()-1-seed_len_g); //-31
			string _rf=hybref.substr(0,it->second);
			if(it->first!=-1) TruncateRead(it->second, _rd, _ref);
			reads[it->first]=_rd+reads[it->first].substr(reads[it->first].size()-1-seed_len_g); //-31
			if(AlignWaterReadnRefRegionT2(H,_rd,_rf,_rd.size(),_rf.size()))
			{
				fw_aligned_rd[it->first]=_rd;
				fw_aligned_rf[it->first]=_rf;
			}
		}
		else
		{
			_mapped_before_ap++;
		}
	}
	for(map<int,int>::iterator it=rc_map.begin(); it!= rc_map.end(); it++)
	{
		if(it->second/*>==*/>=ap)
		{
			string _rd=RCDNA(reads[it->first]).substr(0,reads[it->first].size()-1-seed_len_g); //-31
			string _rf=hybref.substr(0,it->second);
			TruncateRead(it->second, _rd, _ref);
			reads[it->first]=reads[it->first].substr(0,1+seed_len_g+_rd.size()); //+31
			if(AlignWaterReadnRefRegionT2(H,_rd,_rf,_rd.size(),_rf.size()))
			{
				
				rc_aligned_rd[it->first]=_rd;
				rc_aligned_rf[it->first]=_rf;
			}
		}
		else
		{
			_mapped_before_ap++;
		}
	}
	if(mainp!=""&&fw_aligned_rd.count(-1)==0)
	{
		cout << "Can't align main" << endl; return 0;
	}
	cout << "Reads aligned in fw " << fw_aligned_rd.size() << ", in rc " << rc_aligned_rd.size() << endl;
	//cout << "Seeded on 'ref' before 'hap0' starts: " << _mapped_before_ap << endl;
	
	cout << "\nMerging read classes...\n";
	set<int> _unique_rm; set<int> _unique_rm_1;
	map<int,bool> _all_rm;
	for(map<int,string>::iterator it=fw_aligned_rd.begin(); it!=fw_aligned_rd.end(); it++)
		_all_rm[it->first]=false;
	for(map<int,string>::iterator it=rc_aligned_rd.begin(); it!=rc_aligned_rd.end(); it++)
		_all_rm[it->first]=true;
	
	//cout << "\nPrimary haplotypes discovered: " << _unique_rm.size() << "\n";
	//cout << "\nBuilding T-profiles for mapped haplotypes:\n";
	//cout << "Reads matching ref: " << reads_matching_ref << endl;	
	
	for(map<int,bool>::iterator it=_all_rm.begin(); it!=_all_rm.end(); it++)
		_unique_rm_1.insert(it->first);
	
	if(mainp!="") {
		_unique_rm.insert(-1);
		_unique_rm_1.insert(-1);
		cout << "Main path has been specified in run file" << endl;
	}
	map<int,string> aln;
	AlignReadsWithRef(ap,_unique_rm_1,fw_map, rc_map, reads, hybref, aln);
	
	// File name generation
	ostringstream out_tmp_str;
	out_tmp_str << seed_len_g << "_" << seed_step_g << "_" << mismatch_max_g;
	if(cut_val_g<4000) out_tmp_str << "_" << cut_val_g; 
	vector<string> fa_fname; Split(string(argv[2]),'/',fa_fname);
	vector<string> fq_fname; Split(string(argv[1]),'/',fq_fname);	
	string output_fname = fa_fname[fa_fname.size()-1] + "_vs_" + fq_fname[fq_fname.size()-1] + "_" + out_tmp_str.str();
	string output_fname_reads2ref = output_fname+".classes.ref";
	string output_fname_reads2main = output_fname+".classes.main";
	output_fname+=".aln.fasta";
	ofstream o_aln(output_fname.c_str());
	ofstream o_clsref(output_fname_reads2ref.c_str());
	ofstream o_clsmain(output_fname_reads2main.c_str());
	// Output files ready!
	
	map<int,int> sites; set<int> exclusions;
	ClassifyEditingSitesOnGivenPath(hybref, aln[-2], aln[-1], sites);
	for(set<int>::iterator it=_unique_rm_1.begin(); it!=_unique_rm_1.end();it++)
	{
		o_clsref << "Alignment " << *it << " ";
		o_clsmain << "Alignment " << *it << " ";
		bool in_main=false;
		if((*it)>0&&CompareReadWithPath(aln[*it], aln[-1], aln[-2],sites, o_clsref))
		{
			exclusions.insert(*it); in_main=true;
		}
		if(in_main) o_clsmain << "&main&\t";
		CompareReadWithMainPath(aln[*it], aln[-1], aln[-2],sites,o_clsmain);
	}
	
	cout << "Exclusion set size " << exclusions.size() << endl;
	
	int reads_matching_ref=0;
	map<int,float> read_support_values;
	for(map<int,bool>::iterator it=_all_rm.begin(); it!=_all_rm.end(); it++)
		read_support_values[it->first]=1.0;
	
	set<int> matching_ref;
	for(map<int,bool>::iterator it=_all_rm.begin(); it!=_all_rm.end(); it++)
	{
		bool submit=true;
		string _q = reads[it->first];
		if(it->second) _q=RCDNA(_q);
		vector<int> support_vec;
		if(hybref.find(_q)!=string::npos) {
			submit=false;
			reads_matching_ref++;
			matching_ref.insert(it->first);
			if(exclusions.count(it->first))
				exclusions.erase(it->first);
		}
		
		for(map<int,bool>::iterator jt=_all_rm.begin(); jt!=_all_rm.end(); jt++)
		{
			if(exclusions.count(it->first)==0) {
				string _t = reads[jt->first];
				if(jt->second) _t=RCDNA(_t);
				if(it!=jt&&_t.find(_q)!=string::npos)	{
					support_vec.push_back(jt->first);
					submit=false;
				}
			}
		}
		if(exclusions.count(it->first)==0)
		for(int k=0; k<support_vec.size(); k++)
			read_support_values[support_vec[k]]+=1.0/support_vec.size();
		if(submit) _unique_rm.insert(it->first);
	}
	for(map<int,bool>::iterator it=_all_rm.begin(); it!=_all_rm.end(); it++)
	{
		if(exclusions.count(it->first)>0)
		{
			read_support_values[it->first]=exclusions.size();
		}
		if(matching_ref.count(it->first)>0)
		{
			read_support_values[it->first]=matching_ref.size();
		}
	}
	
	read_support_values[-2]=matching_ref.size();
	read_support_values[-1]=exclusions.size();
	
	int i=-1;
	int reads_dd_by_aln_len=0;
	for(map<int,string>::iterator it=aln.begin();it!=aln.end();it++)
	{
		//debug << "READ: " << it->first << " = " << i << endl;
		string __rf = _all_rm[it->first]?rc_aligned_rf[it->first]:fw_aligned_rf[it->first];
		string __rd = _all_rm[it->first]?rc_aligned_rd[it->first]:fw_aligned_rd[it->first];
		//debug <<  __rf << endl << __rd << endl;
		//debug << reads[it->first] << endl << endl;
		int _cp=0;
		while(_cp<__rd.size()&&__rd[_cp]=='-') _cp++;
		
		if(_cp>cut_val_g)
		{
			reads_dd_by_aln_len++;
		}
		else if(matching_ref.count(it->first) ||
			(_unique_rm.count(it->first)==0&&exclusions.count(it->first)==0&&
			it->first>0)
			)
		{
			// Do nothing;
		}
		else
		{
			if(i==0) o_aln << ">main_path";
			else if(i==-1) o_aln << ">REF";
			else o_aln << ">seq" << " [" << it->first << "]";
			if(exclusions.count(it->first)) o_aln << "(*)";
			if(matching_ref.count(it->first)) o_aln << "[ref]";
			o_aln << "_support_val=" << read_support_values[it->first]/*<<",support="<<_5prime_mapping[it->first]*/<<endl<<it->second<<endl;
		}
		i++;
	}
	cout << "Reads discarded by alignment length filter: " << reads_dd_by_aln_len << endl;	
}