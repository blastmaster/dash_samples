#include <iostream>
#include <libdash.h>

#include "dashutils.h"

/**
 * Simple circular shift
 * its necessary to use extra stackvariable if we use global subscript operator
 */
int main(int argc, char** argv)
{
    dash::init(&argc, &argv);
    int myid = dash::myid();
    int size = dash::size();
    dash::Array<int> vals { size };

    vals[myid] = myid;
    if (myid==0) std::cout << "BEFORE\n" << dump(vals);

    vals.barrier();
    int value = vals[(myid == 0) ? size - 1 : myid - 1];
    vals.barrier();
    vals[myid] = value;

    if (myid==0) std::cout << "AFTER\n" << dump(vals);
    dash::finalize();
    return 0;
}
