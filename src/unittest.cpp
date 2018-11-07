#include "unittest.h"
#include "bamutil.h"
#include "pair.h"
#include <time.h>

UnitTest::UnitTest(){

}

void UnitTest::run(){
    bool passed = true;
    passed &= BamUtil::test();
    passed &= diffWithOneGapTest();
    passed &= Pair::test();
    printf("\n==========================\n");
    printf("%s\n\n", passed?"PASSED":"FAILED");
}