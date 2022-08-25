// by J@:)
#ifndef ORFFINDER_H
#define ORFFINDER_H

#include <tuple>
#include <math.h>
#include "taldriver.h"
#include "geneticcode.h"
#include "taoptions.h"

using namespace std;

using overlapgraph = map<int, map<int, int> >;

// assumes filtration by reference id before function call
static vector<TAlignment> FilterMappingDuplicatesReads(vector<TAlignment>& inputAlignments)
{
	vector<TAlignment> filteredAlignments;
	for(int A = 0; A < inputAlignments.size()-1; A++)
	{
		bool isDuplicate = false;
		for(int B = A+1; B < inputAlignments.size(); B++)
		{
			if(inputAlignments[B].ref_S == inputAlignments[A].ref_S &&
			   inputAlignments[B].ref_E == inputAlignments[A].ref_E)
			{
				bool equal = true;
				for(int i = 0; i < inputAlignments[B].dR.size(); i++)
				{
					if(inputAlignments[A].dR[i] != inputAlignments[B].dR[i])
					{
						equal = false;
						break;
					}
				}
				if(equal)
				{
					isDuplicate = true;
					break;
				}
			}
		}
		if(!isDuplicate)
		{
			filteredAlignments.push_back(inputAlignments[A]);
		}
	}
	return filteredAlignments;
}

/*
    It is not obvious, but subreads
    (e.g. reads that are substrings of longer reads)
    should not be filtered out
*/

static inline tuple<int, TAlignment&, TAlignment&, bool>
    OverlapAlignmentPair(TAlignment& left, TAlignment& right)
{
    /*
    Let alignments be:
        +(a)------(b)+   L
            +(c)--------(d)+   R
    then lets denote vector of pairs (a,L), (b,L), (c,R), (d,R).
    sort vector by 1 element, look at the 2 element
    */

    vector<vector<int> > projection = { {left.ref_S, 1},
                                        {left.ref_E, 1},
                                        {right.ref_S, 2},
                                        {right.ref_E, 2} };

    sort(projection.begin(), projection.end(), [](vector<int>& L, vector<int>& R)
        -> bool { return L[0] < R[0]; } );

    /*
        after sorting projection vector can reflect one of three cases:
        1. 1122      +(a)----(b)+
           2211                  +(c)-----(d)+
        2. 1212      +(a)---------(b)+
           2121            +(c)---------(d)+
        3. 1221      +(a)---------(b)+
           2112         +(c)--(d)+
    */
    if(projection[1][1] == projection[2][1])
    {
        if(projection[1][1] == 1)
        {
            return make_tuple(left.ref_E - left.ref_S, std::ref(right), std::ref(left), true);
        }
        else
        {
            return make_tuple(right.ref_E - right.ref_S, std::ref(left), std::ref(right), false);
        }
    }
    else if(projection[0][1] == projection[1][1])
    {
        if(projection[0][1] == 1)
        {
            return make_tuple(0, std::ref(left), std::ref(right), false);
        }
        else
        {
            return make_tuple(0, std::ref(right), std::ref(left), true);
        }
    }
    else
    {
        if(projection[0][1] == 1)
        {
            return make_tuple(left.ref_E - right.ref_S, std::ref(left), std::ref(right), false);
        }
        else
        {
            return make_tuple(right.ref_E - left.ref_S, std::ref(right), std::ref(left), true);
        }
    }
}

static inline bool Joinable(TAlignment& nodeL, TAlignment& nodeR, int overlap)
{
    int deltaX = nodeR.ref_S - nodeL.ref_S;
    bool joinable = true;
    for(int i = 0; i < overlap; i++)
    {
        if(nodeR.dR[i] != nodeL.dR[deltaX + i])
        {
            joinable = false;
            break;
        }
    }
    return joinable;
}

static inline TAlignment MergeReadsOverlap(TAlignment& left, TAlignment& right)
{

    TAlignment mergedAlignment;

    auto overlap_layout = OverlapAlignmentPair(left, right);
	int overlap = get<0>(overlap_layout);
    auto nodeL  = get<1>(overlap_layout);
    auto nodeR  = get<2>(overlap_layout);

    // no overlap
	if(overlap == 0)
	{
		return mergedAlignment;
	}

    // test case '3' from OverlapAlignmentPair()
    /*
        +----------------+ L
           +------+        R
    */
	if(nodeL.ref_E >= nodeR.ref_E)
	{
		// check for overlap and return left
		if(Joinable(nodeL, nodeR, overlap))
		{
			mergedAlignment = nodeL;
		}
	}
    else
    {
    	if(Joinable(nodeL, nodeR, overlap))
    	{
    		mergedAlignment.ref_S = nodeL.ref_S;
    		mergedAlignment.ref_Idx = nodeL.ref_Idx;
    		mergedAlignment.ref_E = nodeR.ref_E;
    		mergedAlignment.flags = nodeL.flags;
            mergedAlignment.setVaild();
            for(int i = 0; i < nodeL.dR.size(); i++)
    		{
    			mergedAlignment.dR.push_back(nodeL.dR[i]);
    		}
    		for(int i = overlap; i < nodeR.dR.size(); i++)
    		{
    			mergedAlignment.dR.push_back(nodeR.dR[i]);
    		}
    	}
    }


	//     AGA A G  C C  G CG  A G C  GACC
	//     012 3 4  5 6  7 89  0 1 2  3456
	//     AGATATG--CTC--GTCG--ATGTC--GACC
	//      GA-A-GTTCTCTTG-CGTTA
	//                  tG-CGTTATG-CTTGA
	//
	//   left:   ref_S = 1; ref_E = 10; dR = [ 0, -1, -1, 2,  0, 2, -1, 0, 2]
	//   right:  ref_S = 7; ref_E = 14; dR = [-1,  0,  2, 0, -1, 2,  0]
	//   joined: ref_S = 1; ref_E = 14; |dR| = 13! ?
	//
	//   deltaX = rightS - leftS = 7-1 = 6
	//   >joinable?
	//   overlap = left.ref_E - right.ref_S = 10 - 7 = 3
	//   do for i = 0...(overlap-1) = 0, 1, 2:
	//    -------------------> (shift for deltaX)
	//   [ 0, -1, -1, 2,  0, 2, -1, 0, 2]
	//                         [-1, 0, 2, 0, -1, 2,  0]
	//                           +  +  +
	//   >merged read
	//   do for i = 0...(deltaX-1) = 0, 1, 2, 3, 4, 5:
	//   [ 0, -1, -1, 2,  0, 2, -1, 0, 2]
	//     +   +   +  +   +  +
	//   do for i = overlap...(right.dR.size()-1) = 3, 4, 5, 6:
	//   [-1,  0,  2, 0, -1, 2,  0]
	//                +   +  +   +
	//

    return mergedAlignment;
}

// returns true/false of whether the ORF is complete
// CDS start/end relative to mRNA sequence,
// mRNA sequence, peptide sequence

static inline tuple<bool, int, int, string, string> FindLongestOrfComplete(
    const string& alignedSeq,
    vector<TlessDNA>& refTless,
    map<string, string>& genetic_code)
{

    map<string, int> start_codons;
    map<string, int> stop_codons;
    for(auto& codon : genetic_code)
    {
        if(codon.second == "start")
        {
            start_codons[codon.first.substr(2)] = 1;
        }
        if(codon.second == "stop")
        {
            stop_codons[codon.first.substr(2)] = 1;
        }
    }

    vector<string> frames(3, "");
    frames[0] = Translate(alignedSeq, genetic_code, 0);
	frames[1] = Translate(alignedSeq, genetic_code, 1);
	frames[2] = Translate(alignedSeq, genetic_code, 2);



    int l_frame = 0;
	int l_start = 0;
	int l_end   = 0;

	int c_start = 0;
	int c_end   = 0;
	int c_frame;

	auto updateLongestOrf = [&]()
    {
		if(c_end - c_start > l_end - l_start)
		{
			l_start = c_start;
			l_end   = c_end;
			l_frame = c_frame;
		}
	};

	for(c_frame = 0; c_frame < 3; c_frame++)
	{
		for(size_t i = 0; i < frames[c_frame].size(); i++)
		{
			if(frames[c_frame][i] == '*')
			{

                c_end     = i;
                c_start   = i-1;

                int p_start   = i-1;
                int p = c_frame + (i * 3) - 3;
                bool start_found = false;
                while(p > 0)
                {
                    if(start_codons.count(alignedSeq.substr(p, 3)))
                    {
                        start_found = true;
                        c_start     = p_start;
                    }
                    if(stop_codons.count(alignedSeq.substr(p, 3)))
                    {
                        break;
                    }
                    p -= 3;
                    p_start--;
                }
                if(start_found)
                {
                    updateLongestOrf();
                }
			}

		}
	}

    bool complete = false;
    // internal cut-off on peptide length (30 aa)
    // function will NOT report any ORFs shorter then 90 nt
    if(l_end - l_start > 30)
    {
        complete = true;
    }

    int cds_begin = l_frame + l_start * 3;
    int cds_end   = 3 * (l_end - l_start + 1);

    return make_tuple(
        complete,                                                               //0
        cds_begin,                                                              //1
        cds_end,                                                                //2
        // CDS is reported from START to STOP inclusive                         //3
        alignedSeq.substr(cds_begin, cds_end),                                  //4
        // peptide is reported as is ('*' truncated)                            //5
        frames[l_frame].substr(l_start, l_end - l_start)
    );
}

map<int, int> CalculateAlignmentSupportIndex(vector<TAlignment>& alignments, bool edited = true)
{
    map<int, int> alignment_supp_index;
    map<int, int> most_supp_ES;
    map<int, map<int, int> > supp_matrix_ES;

    for(auto& align : alignments)
    {
        for(int p = 0 ; p < align.dR.size(); p++)
        {
            supp_matrix_ES[align.ref_S + p][align.dR[p]]++;
        }
    }

    for(auto& dX : supp_matrix_ES)
    {
        int max_val = 0;
        int max_dR  = 0;
        for(auto& dY : dX.second)
        {
            if(edited && dY.first == 0)
            {
                continue;
            }
            if(dY.second > max_val)
            {
                max_val = dY.second;
                max_dR  = dY.first;
            }
        }
        most_supp_ES[dX.first] = max_dR;
    }

    auto alignments_dedup = FilterMappingDuplicatesReads(alignments);
    alignments = alignments_dedup; alignments_dedup.clear();

    for(int i = 0; i < alignments.size(); i++)
    {
        int supp_match = 0;
        for(int p = 0 ; p < alignments[i].dR.size(); p++)
        {
            if(edited && alignments[i].dR[p] == 0)
            {
                continue;
            }
            if(alignments[i].dR[p] == most_supp_ES[alignments[i].ref_S + p])
            {
                supp_match++;
            }
        }
        alignment_supp_index[i] = supp_match;//int(100 * float(supp_match) / float(alignments[i].dR.size()));
    }

    return alignment_supp_index;
}

overlapgraph BuildOverlapGraph(vector<TAlignment>& alignments)
{
	overlapgraph overlapsGraph;

    int filter_ESD =  TAlignerOptions::Options().getOption("ORFinder|OverlapGraph|FilterEsd");
    int min_overlap = TAlignerOptions::Options().getOption("ORFinder|OverlapMin");
    int min_extension = TAlignerOptions::Options().getOption("ORFinder|ExtensionMin");
    string construction_mode =
                      TAlignerOptions::Options().getOption("ORFinder|PreferredMode");
    int edges_count = 0;

    bool prefer_extension = false;
    bool prefer_coverage = false;
    bool prefer_editing  = false;
    map<int, int> alignments_ES_supp;
    if(construction_mode == "extension")
    {
        prefer_extension = true;
    }
    if(construction_mode == "coverage")
    {
        prefer_coverage = true;
        alignments_ES_supp = CalculateAlignmentSupportIndex(alignments);
    }
    if(construction_mode == "coverage_unedited")
    {
        prefer_coverage = true;
        alignments_ES_supp = CalculateAlignmentSupportIndex(alignments, false);
    }
    if(construction_mode == "editing")
    {
        prefer_editing = true;
    }

    // deduplication is not done yet, so lets do it now
    if(construction_mode == "extension" || construction_mode == "editing")
    {
        auto alignments_dedup = FilterMappingDuplicatesReads(alignments);
        alignments = alignments_dedup; alignments_dedup.clear();
    }


	for(int i = 0; i < alignments.size()-1; i++)
	{

        if(alignments[i].ref_ESD < filter_ESD) continue;

		for(int j = i+1; j < alignments.size(); j++)
		{
            if(alignments[j].ref_ESD < filter_ESD) continue;

            auto overlap_layout = OverlapAlignmentPair(alignments[i], alignments[j]);

            int overlap  = get<0>(overlap_layout);
            auto nodeL   = get<1>(overlap_layout);
            auto nodeR   = get<2>(overlap_layout);
            bool swapped = get<3>(overlap_layout);

            int extend = nodeR.dR.size() - overlap;
            int value  = overlap;
            if(prefer_extension)
            {
                value = extend;
            }
            else if(prefer_coverage)
            {
                value = alignments_ES_supp[i];
            }

            if(overlap > min_overlap && extend > min_extension && Joinable(nodeL, nodeR, overlap))
            {
                if(swapped)
                {
                    if(prefer_editing)
                    {
                        value = alignments[i].ref_ESD;
                    }
                    overlapsGraph[j][i] = value;   //{overlap, extend, edits};
                }
                else
                {
                    if(prefer_editing)
                    {
                        value = alignments[j].ref_ESD;
                    }
                    overlapsGraph[i][j] = value;   //{overlap, extend, edits};
                }
                edges_count++;
            }
		}
	}
    overlapsGraph[-1][-1] = edges_count;
    overlapsGraph[-1][-2] = alignments.size();
	return overlapsGraph;
}


void Trace(int I, int N, vector<int> P, overlapgraph& OG,
    vector<vector<int> >& R, const string& mode, const int& sd)
{
	P[N] = I;
	N++;
    bool traced = false;
	if(OG[I].size() == 0)
	{
		vector<int> Orf(N, 0);
		for(int i = 0; i < N; i++)
		{
			Orf[i] = P[i];
		}
		R.push_back(Orf);
	}
	else
	{
		vector<vector<int> > candidates;
		for(auto& neighbour : OG[I])
		{
            candidates.push_back({neighbour.second, neighbour.first});
		}

        if(candidates.size() > 0)
        {
            sort(candidates.begin(), candidates.end(),
            [](const vector<int>& a, const vector<int>& b) {  return a[0] > b[0]; } );

            int depth = sd;
            if(sd < 0)
                //depth = 1 + int(0.22 * N);
                depth = (-1 * sd) + int(0.03 * N) + int(0.01 * pow(N, 2) ) + int(0.002 * pow(N, 3) );

            for(int tr_att = 0; tr_att < depth && tr_att < candidates.size(); tr_att++)
            {
                traced=true;
                Trace(candidates[tr_att][1], N, P, OG, R, mode, sd);
            }
        }
        else if(!traced && N > 1)
        {
            vector<int> Orf(N, 0);
    		for(int i = 0; i < N; i++)
    		{
    			Orf[i] = P[i];
    		}
    		R.push_back(Orf);
        }
	}
}

/*
    Bad decision: this function uses singleton to get the
    parameters 'overlap min' and 'preferred mode' and pass
    them to Trace() engine function
*/
void TraceORF(int I, overlapgraph& OG, vector<vector<int> >& R)
{
	vector<int> P(75, 0);
	int N = 0;
	static string preferredMode = TAlignerOptions::Options().getOption("ORFinder|PreferredMode");
	static int searchDepth      = TAlignerOptions::Options().getOption("ORFinder|Tracing|SearchDepth");
	Trace(I, N, P, OG, R, preferredMode, searchDepth);
}

// returns a tuple (mRNA sequence, peptide sequence, TAlignment object)
// uses FindLongestOrf wrapper with specified genetic code
tuple<string, string, TAlignment, bool> DecodeORF(
    vector<int>& orf, vector<TAlignment>& alignments,
    vector<TlessDNA>& refTless, map<string, string>& genetic_code)
{

    TAlignment ORF = alignments[orf[0]];
    // assemble orf alignment
    for(auto& alignment_id : orf)
    {
        ORF = MergeReadsOverlap(ORF, alignments[alignment_id]);
    }

    string mRNA   = Alignment2Seq(ORF, refTless);
    auto cds_data = FindLongestOrfComplete(mRNA, refTless, genetic_code);
    return make_tuple(get<3>(cds_data), get<4>(cds_data), ORF, get<0>(cds_data));
}

// Auxilary function called to trace a series of ORFs within a thread
// Used by 'StartTracingMT' for multithreaded ORF tracing algo

vector<tuple<string, string, TAlignment> > TraceVector(vector<int>& rids, overlapgraph& og,
    map<string, string>& genetic_code, vector<TAlignment>& Alignments, vector<TlessDNA>& refTless)
{
    vector<vector<int> > R;
    for(int i = 0; i < rids.size(); i++)
    {
        TraceORF(rids[i], og, R);
    }

    vector<tuple<string, string, TAlignment> > orfs_vector;

    for(int i = 0; i < R.size(); i++)
    {
        auto result = DecodeORF(R[i], Alignments, refTless, genetic_code);

        auto cds_seq   = get<0>(result);
        auto pep_seq   = get<1>(result);
        auto alignment = get<2>(result);
        bool complete  = get<3>(result);
        if(complete)
        {
            orfs_vector.push_back({cds_seq, pep_seq, alignment});
        }
    }

    return orfs_vector;
}

void StartTracingMT(vector<TAlignment>& Alignments, vector<TlessDNA>& refTless,
    overlapgraph OG, vector<tuple<string, string, TAlignment> >& ORFs)
{
    int nThreads = thread::hardware_concurrency() - 1;
	int userTheads = TAlignerOptions::Options().getOption("T-Aligner|Threads");
	if(userTheads > 0) nThreads = userTheads;

    string genetic_code_name = TAlignerOptions::Options().getOption("ORFinder|GeneticCode|Table");
    auto genetic_code = setGeneticCode(genetic_code_name);

    vector<vector<int> > rids_to_trace(nThreads);
	vector<future<vector<tuple<string, string, TAlignment> > > >  tracingFutures;

    int jobcounter = 0;
    for(int i = 0; i < Alignments.size(); i++)
    {
        if(Alignments[i].ref_S < refTless[0].T.size() * 0.4)
        {
            rids_to_trace[jobcounter % nThreads].push_back(i);
            jobcounter++;
        }
    }

	for(int i = 0; i < rids_to_trace.size(); i++)
	{
        if(rids_to_trace[i].size() > 0)
					tracingFutures.push_back(
                        async(launch::async, TraceVector,
                        std::ref(rids_to_trace[i]), std::ref(OG),
                        std::ref(genetic_code), std::ref(Alignments), std::ref(refTless) ));
	}
	for(auto &tracingTask : tracingFutures)
	{
				auto data = tracingTask.get();
        if(data.size() > 0) ORFs.insert(ORFs.end(), data.begin(), data.end());
	}
}



#ifdef DOUNITTEST


void JoinOverlapped(string ref_seq = "AAAGGTTTGGATTTCGCCGATGC",
	string read1_seq="AGTTGGTTGTTACTTGCTTTCG",
	string read2_seq="TGTTGTTACTTGCTTTCGAGC")
{
		// this function tests if the read merging functions work
		// properly

		// it will create and re-write files pseudoref.fa and pseudoread.fq,
		// which are the test files

    ofstream ofasta("pseudoref.fa");
    ofasta << ">Psref\n" << ref_seq << "\n";
    ofasta.close();

    ofstream ofastq("pseudoread.fq");
    ofastq << "@rd1\n" << read1_seq << "\n+\nread1qualsdummy\n";
    ofastq << "@rd2\n" << read2_seq << "\n+\nread2qualsdummy\n";
    ofastq.close();

    const char* referenceFile = "pseudoref.fa";
    const char* fastqFile     = "pseudoread.fq";

    // set options for seeding small sequence

    auto& TOP = TAlignerOptions::Options();
    TOP.setOption("Aligner|Seed", "4");

    vector<vector<string> > refs;
    vector<TlessDNA> refTless;
    unordered_map<string, vector<Rindex> > refIndex;

		// build reference index
    ProcessReferenceFile(referenceFile, refs, refTless);
    Build_Fasta_Index(refTless, refIndex);

    int rfid = 0;

		// print out the reference Tless profile
    cout << "Dumping reference:\n" << DumpTless(refTless[rfid]) << "\n" <<
    "Printing reference:\n" << PrintTless(refTless[rfid]) << "\nDone testing reference!\n\n";


    std::ifstream ifile(fastqFile); std::string line1, line2, line3, line4;

    // aligning read 1
    std::getline(ifile, line1, '\n');
	std::getline(ifile, line2, '\n');
	std::getline(ifile, line3, '\n');
	std::getline(ifile, line4, '\n');

		ReadAlignmentTask RAT;
    RAT.refIndex = &refIndex;
    RAT.refTless = &refTless;
    RAT.read     = MakeTless(line2);
    RAT.aligned_segment_min = 8;

    auto Alignment1 = TAlignRead(RAT);

    // aligning read 2
    std::getline(ifile, line1, '\n');
	std::getline(ifile, line2, '\n');
	std::getline(ifile, line3, '\n');
	std::getline(ifile, line4, '\n');

		// same RAT is used, only read sequence is needed to be changed
    RAT.read     = MakeTless(line2);
    auto Alignment2 = TAlignRead(RAT);


		// printing read alignments
    cout << "Printing read1:\n" << Alignment2Seq(Alignment1, refTless) << "\nDone!\n\n";
    cout << "Printing read2:\n" << Alignment2Seq(Alignment2, refTless) << "\nDone!\n\n";

    cout << "Printing the alignment 1 itself:\n" <<
    PrintAlignedReadFasta(Alignment1, refTless) << "\nDone!\n\n";
    cout << "Printing the alignment 2 itself:\n" <<
    PrintAlignedReadFasta(Alignment2, refTless) << "\nDone!\n\n";

    cout << "Printing TAF v2 lines:\n" << PrintReadTaf(Alignment1, refTless) << "\n" << PrintReadTaf(Alignment2, refTless) << "\nDone!\n\n";


	// testing auxilary functions 'OverlapAlignmentPair' and 'Joinable'
	// both functions are actually called inside main function 'MergeReadsOverlap'

	auto OverlapLR = OverlapAlignmentPair(Alignment1, Alignment2);
	auto OverlapRL = OverlapAlignmentPair(Alignment2, Alignment1);

	cout << "Overlap L-R: " << get<0>(OverlapLR) << ", suggest pair reverse is " << get<3>(OverlapLR) << "\n";
	cout << "Overlap R-L: " << get<0>(OverlapRL) << ", suggest pair reverse is " << get<3>(OverlapRL) << "\n";

	cout << "\n";

	cout << "L-R is joinable " << Joinable(Alignment1, Alignment2, get<0>(OverlapLR)) << ", R-L is joinable " << Joinable(Alignment2, Alignment1, get<0>(OverlapRL)) << "\n\n";

	cout << "\n";

	// testing 'MergeReadsOverlap' main function

	TAlignment overlap = MergeReadsOverlap(Alignment1, Alignment2);
	cout << "Printing merged read sequence:\n" << Alignment2Seq(overlap, refTless) << "\nDone!\n\n";
	cout << "Printing TAF v2 line for merged:\n" << PrintReadTaf(overlap, refTless) << "\nDone!\n\n";

    cout << "Printing merged alignment data:\nref_S=" << overlap.ref_S << ", ref_E=" << overlap.ref_E <<
        ", size of dR=" << overlap.dR.size() << "\ndR vector: ";
    for(int i =0; i < overlap.dR.size(); i++) cout << overlap.dR[i] << "; ";
    cout << "\nDone!\n\n";

	cout << "Printing the alignment merged itself:\n" <<
    PrintAlignedReadFasta(overlap, refTless) << "\nDone!\n\n";


    cout << "All done!\n";


}

#endif

#endif
