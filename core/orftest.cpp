#include "aligner.h"
#include "orffinder.h"



int main(int argc, char** argv)
{

    // default test

    /***  COMMENTED OUT TEST 1 ***
    JoinOverlapped();
    */

    // test 2.
    //
    //   REF            AG--G---AG---ATTA--CGAGTG---GTTTA-A-G---CTTTGAG---C---CACAAG
    //   READ1          AGTTGTTTAGTTTA--ATTCGAG-GTTTG-TTATATGTTTCT--G
    //   READ2                                G-GTTTG-TTATATGTTTCT--GAGTTTCTTTCA

    /***  COMMENTED OUT TEST 2 ***
    JoinOverlapped("AGGAGATTACGAGTGGTTTAAGCTTTGAGCCACAAG",
    "AGTTGTTTAGTTTAATTCGAGGTTTGTTATATGTTTCTG",
    "GGTTTGTTATATGTTTCTGAGTTTCTTTCA");
    */

    // test 3.
    //
    //   REF    TTTC---GA---C--A--GTTG--AC--G--GATTTCC-A-A---CGTGG---AC---GCAGCAGTTTTT
    //   READ1    TCTTTGATTTCTTATTG--GTTACTTGTT
    //   READ2                 ATTG--GTTACTTGTTGA---CCTATATTTCG-GGTTTACTTT

    JoinOverlapped("TTTCGACAGTTGACGGATTTCCAACGTGGACGCAGCAGTTTTT",
    "TCTTTGATTTCTTATTGGTTACTTGTT",
    "ATTGGTTACTTGTTGACCTATATTTCGGGTTTACTTT");
}
