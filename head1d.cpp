#include <iostream>
#include <memory>
#include <stdint.h>
#include <libdash.h>

#include <mpi.h>

#include "dashutils.h"

/**
 * 1D heat equation using dash library
 * tested with MPI-3.0 dart runtime
 * usage: $ mpirun -n #procs ./heat [problemsize]
 * problemsize defaults to 16
 */

class Heatmap
{

public:
    dash::Pattern pat;
    std::unique_ptr<dash::Array<double>> u;
    std::unique_ptr<dash::Array<double>> v;
    dash::Array<double> a;

    Heatmap(long long n, dash::DistSpec ds = dash::BLOCKED) : pat(n, ds),
        u(new dash::Array<double>(n, ds)),
        v(new dash::Array<double>(n, ds)),
        a(n, ds)
    { }

    ~Heatmap() {
        dash::finalize();
    }

    void initialize(const double alpha=22.0,
                    const double u_min=10.0,
                    const double u_max=100.0)
    {

        int myid = dash::myid();
        int size = dash::size();

        if ( myid < (int)(size / 2) ) {
            std::fill(u->lbegin(), u->lend(), u_max);
            std::fill(v->lbegin(), v->lend(), u_max);
        }
        else {
            std::fill(u->lbegin(), u->lend(), u_min);
            std::fill(v->lbegin(), v->lend(), u_min);
        }

        std::fill(a.lbegin(), a.lend(), alpha);
        dash::TeamAll.barrier();
    }
};

void printstate(const double residual, const uint64_t step) {
    std::cout << "Step: " << step << " residual: " << residual << std::endl;
}

int main(int argc, char* argv[])
{
    dash::init(&argc, &argv);
    int myid = dash::myid();
    int size = dash::size();

    size_t ps = 16;

    if (argc >= 2) {
        ps = atoi(argv[1]);
        std::cout << "Using problem size: " << ps << std::endl;
    }
    else {
        std::cout << "Using default problem size: " << ps << std::endl;
    }

    double alpha     = 22.0;
    double l         = 2.0;
    double h         = l / (double) ps;
    double hi        = 1.0 / h;
    double hi2       = hi * hi;
    double dt        = h * h / 4 / alpha;
    double eps       = 1.0e-4;
    double residual  = eps + 1.0;

    Heatmap ht_map(ps, dash::BLOCKED);
    ht_map.initialize();

    if ( myid == 0 ) {
        std::cout << "After init\n";
        std::cout << "U\n" << dump(*ht_map.u);
        std::cout << "V\n" << dump(*ht_map.v);
    }

    uint64_t step = 0;
    dash::Array<double> borders(size * 2);

    dash::TeamAll.barrier();
    borders_from_values(borders, *ht_map.u);
    borders_iexchange(borders);

    dash::TeamAll.barrier();
    int blksz = ht_map.u->size() / size;
    if(myid==0)
        std::cout << "blksz: " << blksz << "\n";

    while ( residual > eps ) {
        residual = 0.0;
        auto it_u = ht_map.u->lbegin();
        auto it_v = ht_map.v->lbegin();
        auto it_a = ht_map.a.lbegin();
        double after, before = 0;
        for (int i = 0; i < blksz; ++i, ++it_u, ++it_v, ++it_a)
        {
            if ( it_u == --ht_map.u->lend() )
                after = *(--borders.lend());
            else
                after = *(it_u + 1);

            if ( it_u == ht_map.u->lbegin() )
                before = *borders.lbegin();
            else
                before = *(it_u - 1);

            if ( it_u != ht_map.u->begin() && it_u != --ht_map.u->end() ) {
                double du = ( before * *it_a + after * *(it_a + 1) -
                        *it_u * (*it_a + *(it_a + 1))) * dt * hi2;
                *it_v = du + *it_u;
                du = std::max(du, -du);
                residual = std::max(residual, du);
            }
        }

        std::swap(ht_map.u, ht_map.v);
        step++;

        dash::TeamAll.barrier();
        borders_from_values(borders, *ht_map.u);
        borders_iexchange(borders);

        MPI_Allreduce(&residual, &residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if ((step % 1000) == 0 && myid == 0)
            printstate(residual, step);
    }

    if (myid == 0)
        std::cout << "Done ... steps: " << step << std::endl;

    if ( myid == 0 ) {
        std::cout << "Results:\n";
        std::cout << "U\n" << dump(*ht_map.u);
        std::cout << "V\n" << dump(*ht_map.v);
    }

    return 0;
}
