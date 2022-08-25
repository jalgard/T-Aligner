
#include "orffinder.h"
#include <fstream>


int FindOrfs()
{

    vector<vector<string> > refsBank;
    vector<TlessDNA> refTless;
    unordered_map<string, vector<Rindex> > refIndex;

    cerr << "Starting protocol: ORF finding\n";

    string inFastaFile = TAlignerOptions::Options().getOption("Input|Reference|Filename");
    ProcessReferenceFile(inFastaFile.c_str(), refsBank, refTless);
    Build_Fasta_Index(refTless, refIndex);

    ReadAlignmentTask RAT;
    RAT.refIndex = &refIndex;
    RAT.refTless = &refTless;

    int nThreads = thread::hardware_concurrency() - 1;
	int userThreads = TAlignerOptions::Options().getOption("T-Aligner|Threads");
    if(userThreads > 0) nThreads = userThreads;

    cerr << "Alignment stage: using " << nThreads << " threads\n";
    string fastqInputFile = TAlignerOptions::Options().getOption("Input|Library1|Filename");
    auto AlignmentsBank = SingleEndAlignFastqMT(fastqInputFile.c_str(), RAT);

    cerr << "Alignment stage: counting reads matching reference...\n";

    int matchingRef = 0;
    for(int i = 0; i < AlignmentsBank.size(); i++)
    {
        if(AlignmentsBank[i].doesMatchRef()) matchingRef++;
    }



    cerr << "Alignment stage:\tAligned:" << AlignmentsBank.size() << "\tMatching ref:" << matchingRef << "\n";
    cerr << "@TAL:\tAligned:\t" << AlignmentsBank.size() << "\tMatching ref:\t" << matchingRef << "\n";

    // Deduplication step was moved inside 'BuildOverlapGraph' function

    cerr << "Assembly stage: building overlap graph... ";
    auto OGr = BuildOverlapGraph(AlignmentsBank);    // previous version accepted 'dedupAlignBank' after read duplications
    cerr << "final graph has " << OGr[-1][-1] << " edges.\n";
    cerr << "Alignments after read merging: " << OGr[-1][-2] << ".\n";

    string prefMode = TAlignerOptions::Options().getOption("ORFinder|PreferredMode");
    cerr << "Assembly stage: searching ORFs with " << prefMode << " algorithm.\n";

    // Multi-threaded ORF tracer

    vector<tuple<string, string, TAlignment> > ORFs, unique_ORFs;
    StartTracingMT(AlignmentsBank, refTless, OGr, ORFs);

    cerr << "Total ORFs traced: " << ORFs.size() <<  "\nRemoving exact duplicates...\n";
    map<string, int> orf_seen;

    for(auto& orf : ORFs)
    {
        string pep_seq = get<1>(orf);
        if(orf_seen[pep_seq] == 0)
        {
            orf_seen[pep_seq] = 1;
            unique_ORFs.push_back(orf);
        }
    }

    cerr << "Unique ORFs left: " << unique_ORFs.size() << "\nSorting ORFs...\n";

    sort(unique_ORFs.begin(), unique_ORFs.end(),
        [](const tuple<string, string, TAlignment>& a,
            const tuple<string, string, TAlignment>& b)
        {  return get<1>(a).size() > get<1>(b).size(); });

    cerr << "Writing output files...\n";
    string output_prefix = TAlignerOptions::Options().getOption("Output|Prefix");
    string output_cds_fasta = TAlignerOptions::Options().getOption("Output|Files|OrfFinderCdsFasta");
    string output_pep_fasta = TAlignerOptions::Options().getOption("Output|Files|OrfFinderPepFasta");
    string output_mrna_fasta = TAlignerOptions::Options().getOption("Output|Files|OrfFinderMrnaFasta");
    string output_mrna_taf = TAlignerOptions::Options().getOption("Output|Files|OrfFinderMrnaTaf");

    ofstream ofile_cds_fasta((output_prefix + output_cds_fasta).c_str());
    ofstream ofile_pep_fasta((output_prefix + output_pep_fasta).c_str());
    ofstream ofile_mrna_fasta((output_prefix + output_mrna_fasta).c_str());
    ofstream ofile_mrna_taf((output_prefix + output_mrna_taf).c_str());

    int orf_name_suffix = 1;
    string reference_name_prefix = refsBank[0][0];

    string print_each_orf_alignment =
        TAlignerOptions::Options().getOption("ORFinder|OutputEachOrfAlignment");

    string orf_mode = TAlignerOptions::Options().getOption("ORFinder|PreferredMode");
    int orf_length_min_aa = TAlignerOptions::Options().getOption("ORFinder|MinOrfLength|aa");

    for(auto&& orf : unique_ORFs)
    {
        if(get<1>(orf).size() < orf_length_min_aa) continue;

        ofile_cds_fasta << ">" << reference_name_prefix << "_" << orf_mode << "_ORF_" << orf_name_suffix << "\n";
        ofile_pep_fasta << ">" << reference_name_prefix << "_" << orf_mode << "_ORF_" << orf_name_suffix << " Length=" << get<1>(orf).length() << "\n";
        ofile_mrna_fasta << ">" << reference_name_prefix << "_" << orf_mode << "_mRNA_" << orf_name_suffix << "\n";

        ofile_cds_fasta << get<0>(orf) << "\n";
        ofile_pep_fasta << get<1>(orf) << "\n";
        ofile_mrna_fasta << Alignment2Seq(get<2>(orf), refTless) << "\n";

        ofile_mrna_taf << PrintReadTaf(get<2>(orf), refTless) << "\n";

        if(print_each_orf_alignment == "true")
        {

            ofstream ofile_alignment((output_prefix + "orf_" +
                std::to_string(orf_name_suffix) + "_alignment.fasta").c_str());
            string orf_name = "ORF" + std::to_string(orf_name_suffix);
            ofile_alignment << PrintAlignedReadFasta(get<2>(orf), refTless, "ref", orf_name);
            ofile_alignment.close();
        }

        orf_name_suffix++;
    }

    ofile_cds_fasta.close();
    ofile_pep_fasta.close();
    ofile_mrna_fasta.close();
    ofile_mrna_taf.close();

    return 0;

}


int main(int argc, char** argv)
{
    cerr << ParseFromCommandLine(argc, argv) << "\n";
    int result = FindOrfs();
    return result;
}
