/*
    Updated kinetoplastid gRNA finding
    algorithm based on the seed extension

    by default seed of length 25 is extended
    in both directions

    Written by Evgeny Gerasimov, 2022
*/

#include <algorithm>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <map>
#include <fstream>
#include <unordered_map>

using namespace std;

void ReadCommandLine(int argc, char** argv, map<string, string>& keys)
{
    for(int argId = 0; argId < argc; argId++)
    {
        string key = string(argv[argId]);
        if(key == "--mm")
        {
            keys["MaxMismatch"] = string(argv[++argId]);
        }
        else if(key == "--gu")
        {
            keys["MaxG:U"] = string(argv[++argId]);
        }
        else if(key == "--seed_score")
        {
            keys["SeedScoreMin"] = string(argv[++argId]);
        }
        else if(key == "--seed_length")
        {
            keys["SeedLength"] = string(argv[++argId]);
        }
        else if(key == "--grna")
        {
            keys["GrnaFile"] = string(argv[++argId]);
        }
        else if(key == "--mrna")
        {
            keys["MrnaFile"] = string(argv[++argId]);
        }
        else if(key == "--score")
        {
            keys["AlignScoreMin"] = string(argv[++argId]);
        }
        else if(key == "--length")
        {
            keys["AlignLengthMin"] = string(argv[++argId]);
        }
        else if(key == "--flanks")
        {
            keys["FlankLength"] = string(argv[++argId]);
        }
        else if(key == "--anchor")
        {
            keys["MinAnchor"] = string(argv[++argId]);
        }
        else if(key == "--anchor_gu")
        {
            keys["AllowG:uInAnchor"] = "true";
        }
    }

}


inline vector<string> split(string str, char dlm)
{
    stringstream ss(str);
    string item;
    vector<string> splitted;
    while (std::getline(ss, item, dlm))
    {
       splitted.push_back(item);
    }
    return splitted;
}

inline string UpperDNA(const string& DNA)
{
    string upperDNA = DNA;
    transform(upperDNA.begin(), upperDNA.end(), upperDNA.begin(), ::toupper);
    return upperDNA;
}


static inline vector<vector<string> > FastaReader(const char* filename, bool split_space = true)
{
    vector<vector<string> > entries;

    string line;
    ifstream ifile(filename);
    while(getline(ifile, line, '\n'))
    {
        if(line[line.size()-1] == '\r')
        {
            line.pop_back();
        }
        if(line[0] == '>')
        {
            string name = line.substr(1);
            if(split_space)
            {
                name = split(line.substr(1), ' ')[0];
            }
            entries.push_back({name, ""});
        }
        else
        {
            entries.back()[1] += UpperDNA(line);
        }
    }
    return entries;
}

vector<int> CheckSegmentMatch(const string& gRnaMer, const string& mRnaMer)
{
    int merScore    = 0;
    int guPairCount = 0;
    int mmPairCount = 0;

    for(int i = 0; i < gRnaMer.size(); i++)
    {
        if(gRnaMer[i] == mRnaMer[i])
        {
            merScore += 2;
        }
        else if( (gRnaMer[i] == 'A' && mRnaMer[i] == 'G') ||
                 (gRnaMer[i] == 'C' && mRnaMer[i] == 'T') )
        {
            merScore += 1;
            guPairCount++;
        }
        else
        {
            merScore -= 2;
            mmPairCount++;
        }
    }

    vector<int> merResult = { merScore, guPairCount, mmPairCount };

    return merResult;
}


// Trim anchor region to the right

void TrimToAnchor(const string& gRnaText, int& gl, int& gr,
                  const string& mRnaText, int& ml, int& mr,
                  int& gu, int& mm, int minAnchor = 5, bool allow_gu_anchor = false)
{
    int cAnchor = 0;


    // Scan right-to-left and count only exatch matches in 5+ anchor
    while(gr > gl && mr > ml)
    {
        if(gRnaText[gr] == mRnaText[mr])
        {
            cAnchor++;
        }
        else if(allow_gu_anchor && ((gRnaText[gr] == 'A' && mRnaText[mr] == 'G') ||
                (gRnaText[gr] == 'C' && mRnaText[mr] == 'T')) )
        {
            cAnchor++;
        }
        else
        {
            cAnchor = 0;
        }
        if(cAnchor == minAnchor)
        {
            break;
        }
        gr--; mr--;
    }
    gr += minAnchor;
    mr += minAnchor;
}

inline bool ValidateHit(const string& alnStr, int gu, int mm, int minScore = 32, int minLength = 24)
{
    int score = (-2 * mm) + gu + (2.0 * (alnStr.length() - gu - mm));

    // normal mode: alnStr.length() should be > 25
    if(score > minScore && alnStr.length() > minLength)  return true;

    return false;
}

void PrintMatch(const string& gRnaName, const string& mRnaName,
                const string& gRnaText, int gl, int gr,
                const string& mRnaText, int ml, int mr, int gu, int mm,
                int flankLength = 50, int minScore = 32, int minLength = 24, int minAnchor = 5,
                bool allow_gu_anchor = false)
{

    TrimToAnchor(gRnaText, gl, gr, mRnaText, ml, mr, gu, mm, minAnchor, allow_gu_anchor);

    string alnStr = "";

    for(int i = 0; i < gr-gl; i++)
    {
        if(gRnaText[i+gl] == mRnaText[i+ml])
        {
            alnStr += "|";
        }
        else if( (gRnaText[i+gl] == 'A' && mRnaText[i+ml] == 'G') ||
                 (gRnaText[i+gl] == 'C' && mRnaText[i+ml] == 'T') )
        {
            alnStr += ":";
        }
        else
        {
            alnStr += "#";
        }
    }

    if(ValidateHit(alnStr, gu, mm, minScore, minLength))  // normal mode: deault min score = 32
    {
        cout << gRnaName << "\t" << gl << "\t" << gr << "\t" << gRnaText.substr(gl, gr-gl) << "\t" <<
                mRnaName << "\t" << ml << "\t" << mr << "\t" << mRnaText.substr(ml, mr-ml) << "\t" <<
                gu << "\t" << mm << "\t" << gr-gl <<  "\n";

        cout << "## Printing alignment:\n##\t" << gRnaText.substr(gl, gr-gl) << "\n##\t" <<
                alnStr << "\n##\t" << mRnaText.substr(ml, mr-ml) << "\n\n";

        // report upstream sequence
        if(gl - flankLength > 0)
        {
            cout << "@@ Upstream " << flankLength <<  "bp: " << gRnaText.substr(gl - flankLength, flankLength) << "\n";
        }
        else
        {
            cout << "@@ Upstream: 5' end" << "\n";
        }
        if(gr + flankLength < gRnaText.size())
        {
            cout << "@@ Downstream " << flankLength << "bp: " << gRnaText.substr(gr, flankLength) << "\n\n";
        }
        else
        {
            cout << "@@ Downstream: 3' end" << "\n";
        }

    }


}

struct gRnaHit
{
    int gLeft;
    int gRight;
    int mLeft;
    int mRight;
    int gu;
    int mm;
};

vector<gRnaHit> SearchForGuiding(
    const string& gRnaText, const string& mRnaText,
    const int kSeed = 25, const int kStep=4, const int sCrit = 32,
    // the critical score for 25 nt match is 32, which is roughly
    // 13 matches + 10 G:U pairs + 2 mismatches
    const int guMax = 15, const int mmMax = 4)

{

    vector<gRnaHit> gRnaHits;

    for(int iGrna = 0; iGrna < gRnaText.size() - kSeed; iGrna++)
    {
        for(int iMrna = 0; iMrna < mRnaText.size() - kSeed; iMrna++)
        {
            // match pair of segments
            vector<int> match = CheckSegmentMatch(gRnaText.substr(iGrna, kSeed),
                                            mRnaText.substr(iMrna, kSeed));
            // try to extend match
            if(match[0] > sCrit)
            {
                int segLeftG  = iGrna;
                int segRightG = iGrna + kSeed;
                int segLeftM  = iMrna;
                int segRightM = iMrna + kSeed;

                int gu = match[1];
                int mm = match[2];
                while(gu < guMax && mm < mmMax)
                {
                    vector<int> leftStep(3,0), rightStep(3,0);

                    if(segLeftG - kStep > 0 && segLeftM - kStep > 0)
                    {
                        leftStep = CheckSegmentMatch(gRnaText.substr(segLeftG - kStep, kStep),
                                               mRnaText.substr(segLeftM - kStep, kStep));
                    }

                    if(segRightG + kStep < gRnaText.size() &&
                        segRightM + kStep < mRnaText.size())
                    {
                        rightStep = CheckSegmentMatch(gRnaText.substr(segRightG, kStep),
                                               mRnaText.substr(segRightM, kStep));
                    }


                    if(rightStep[0] <= 0 && leftStep[0] <= 0)
                    {
                        break;
                    }
                    else if(leftStep[0] > rightStep[0])
                    {
                        gu += leftStep[1];
                        mm += leftStep[2];
                        segLeftG -= kStep;
                        segLeftM -= kStep;
                    }
                    else
                    {
                        gu += rightStep[1];
                        mm += rightStep[2];
                        segRightG += kStep;
                        segRightM += kStep;
                    }
                }

                gRnaHit gRna;
                gRna.gLeft  = segLeftG;
                gRna.gRight = segRightG;
                gRna.mLeft  = segLeftM;
                gRna.mRight = segRightM;
                gRna.gu=gu;
                gRna.mm=mm;

                gRnaHits.push_back(gRna);
            }
        }
    }

    return gRnaHits;

}

void ClusterAndPrintHits(
const vector<gRnaHit>& gRnaHits,
const string& gRnaName, const string& mRnaName,
const string& gRnaText, const string& mRnaText,
int flankLength = 50, int minScore = 32, int minLength = 24, int minAnchor = 5,
bool allow_gu_anchor = false)
{
    // vector of gRnaHits is assumed to be sorted
    // by start coordinate in 'minicircle' scaffold
    if(gRnaHits.size() > 1)
    {
        // create overlapping clusters
        vector<vector<gRnaHit> > hitClusters;

        vector<gRnaHit> clr = { gRnaHits[0] };
        int lastRight = gRnaHits[0].gRight;
        for(int i = 1; i < gRnaHits.size(); i++)
        {
            if(gRnaHits[i].gLeft < lastRight)
            {
                clr.push_back(gRnaHits[i]);
            }
            else
            {
                hitClusters.push_back(clr);
                clr.clear();
                clr.push_back(gRnaHits[i]);
            }
            lastRight = gRnaHits[i].gRight;
        }

        for(int c = 0; c < hitClusters.size(); c++)
        {
            int bestL = 0;
            int bestI = 0;
            for(int i = 0; i < hitClusters[c].size(); i++)
            {
                gRnaHit& grna = hitClusters[c][i];
                int l = grna.gRight - grna.gLeft;
                if(l > bestL)
                {
                    bestI = i;
                    bestL = l;
                }
            }
            gRnaHit& grna = hitClusters[c][bestI];
            PrintMatch(gRnaName, mRnaName,
                       gRnaText, grna.gLeft, grna.gRight,
                       mRnaText, grna.mLeft, grna.mRight,
                       grna.gu, grna.mm,
                       flankLength, minScore, minLength, minAnchor, allow_gu_anchor);
        }
    }
    else
    {
        const gRnaHit& grna = gRnaHits[0];
        PrintMatch(gRnaName, mRnaName,
                   gRnaText, grna.gLeft, grna.gRight,
                   mRnaText, grna.mLeft, grna.mRight,
                   grna.gu, grna.mm,
                   flankLength, minScore, minLength, minAnchor, allow_gu_anchor);
    }
}

int main(int argc, char** argv)
{
    map<string, string> commandLine;
    /*
        Defaults
    */
    commandLine["SeedLength"]="25";
    commandLine["SeedScoreMin"]="32";
    commandLine["MaxG:U"]="14";
    commandLine["MaxMismatch"]="4";
    commandLine["FlankLength"]="50";
    commandLine["AlignScoreMin"]="32";
    commandLine["AlignLengthMin"]="26";
    commandLine["MinAnchor"]="5";
    commandLine["AllowG:uInAnchor"]="false";

    ReadCommandLine(argc, argv, commandLine);


    cerr << "Starting G-Finder (T-Aligner suite) with options: ";
    for(auto& key : commandLine)
    {
        cerr << " [";
        cerr << key.first << " -> " << key.second << "]; ";
    }
    cerr << "\n\n";

    auto gRnaFasta = FastaReader(commandLine["GrnaFile"].c_str(), false);
    auto mRnaFasta = FastaReader(commandLine["MrnaFile"].c_str(), false);

    cerr << "Reading gRNA targets from " << commandLine["GrnaFile"] << "\n";
    cerr << "Reading mRNAs from " << commandLine["MrnaFile"] << "\n";


    bool allow_gu_anchor = false;
    if(commandLine["AllowG:uInAnchor"] == "true")
    {
        allow_gu_anchor = true;
    }

    for(auto& mRna : mRnaFasta)
    {
        cerr << "Searching guide RNAs for mRNA " << mRna[0] << " sequence of length " << mRna[1].size() << " bp.\n";
        for(auto& gRna : gRnaFasta)
        {
            cerr << "...   target molecule >" << gRna[0] << "< of length " << gRna[1].size() << " bp.\n";
            vector<gRnaHit> rawHits = SearchForGuiding(gRna[1], mRna[1],
                    atoi(commandLine["SeedLength"].c_str()), 4,
                    atoi(commandLine["SeedScoreMin"].c_str()),
                    atoi(commandLine["MaxG:U"].c_str()),
                    atoi(commandLine["MaxMismatch"].c_str()));


            if(rawHits.size() > 0)
            {
                cerr <<  "......   clustering and printing hits. \n\n";
                ClusterAndPrintHits(rawHits, gRna[0], mRna[0], gRna[1], mRna[1],
                    atoi(commandLine["FlankLength"].c_str()),
                    atoi(commandLine["AlignScoreMin"].c_str()),
                    atoi(commandLine["AlignLengthMin"].c_str()),
                    atoi(commandLine["MinAnchor"].c_str()),
                    allow_gu_anchor);
            }
            else
            {
                cerr << "...... no hits found. \n\n";
            }
        }
    }
    cerr << "All done!\n";
}
