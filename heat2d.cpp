#include <iostream>
#include <memory>
#include <stdint.h>

#include "libdash.h"
#include "dashutils.h"

#include <mpi.h>

using dash::size_type;

template<class PATTERN>
class Heatmap
{
    private:
        std::unique_ptr<dash::Matrix<double, PATTERN>> m_u;
        std::unique_ptr<dash::Matrix<double, PATTERN>> m_v;

    public:
        Heatmap(const size_type rows=0, const size_type columns=0) :
            m_u(new dash::Matrix<double, PATTERN>(rows, columns)),
            m_v(new dash::Matrix<double, PATTERN>(rows, columns))
        {}

        ~Heatmap() {}

        void
        initialize_data(const double hi_v=100.0,
                        const double lo_v=10.0,
                        const double alpha=0.0)
        {
            //int myid = dash::myid();
            //int n_units = dash::size();

            //auto size = m_u->size();
            //if (myid < n_units / 2) {
                //std::fill(m_u->lbegin(), m_u->lend(), hi_v);
                //std::fill(m_v->lbegin(), m_v->lend(), hi_v);
            //}
            //else {
                //std::fill(m_u->lbegin(), m_u->lend(), lo_v);
                //std::fill(m_v->lbegin(), m_v->lend(), lo_v);
            //}

            auto rows = m_u->rows();
            for (int r = 0; r < rows; ++r) {
                char is_hot = (r < (rows / 3));
                for (int c = 0; c < m_u->cols(); ++c) {
                    if (is_hot) {
                        m_u->operator()(r, c) = hi_v;
                        m_v->operator()(r, c) = hi_v;
                    }
                    else {
                        m_u->operator()(r, c) = lo_v;
                        m_v->operator()(r, c) = lo_v;
                    }
                }
            }

            m_u->barrier();
            m_v->barrier();

        }

        template<typename STENCIL>
        void compute(STENCIL &s, double *residual, uint64_t& steps=0)
        {
            double l        = 2.0;
            double alpha    = 22.0;
            double h        = l / (double) m_u->size();
            double hi       = 1.0 / h;
            double hi2      = hi * hi;
            double dt       = h * h / 4 / alpha;

            *residual = 0.0;

            auto viter = m_v->st_begin(s);
            for (auto iter = m_u->st_begin(s);
                        iter != m_u->st_end(s);
                        ++iter, ++viter)
            {
                double du = ( *iter.west() + *iter.east() +
                        *iter.north() + *iter.south() -
                        4.0 * *iter.self() ) * alpha * dt * hi2;

                *viter.self() = *iter.self() + du;
                du = std::max(du, -du);
                *residual = std::max(*residual, du);
            }
            std::swap(m_u, m_v);
            steps++;
            if (dash::myid() == 0 && steps < 4)
                std::cout << pretty_print(*m_u) << "\n";
        }
};


int main(int argc, char** argv)
{
    dash::init(&argc, &argv);

    double eps      = 1.0e-4;
    double residual = 0.0;
    static uint64_t step = 0;
    size_type rows, columns = 8;

    if (argc >= 3) {
        rows = atoi(argv[1]);
        columns = atoi(argv[2]);
    }

    Heatmap<dash::Pattern2DRowLayout> hm { rows, columns };
    hm.initialize_data();
    dash::TeamAll.barrier();

    dash::Stencil_5P stencil;

    do
    {
        hm.compute(stencil, &residual, step);
        MPI_Allreduce(&residual, &residual, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    } while (residual > eps);


    using std::cout;
    if (dash::myid() == 0) {
        cout << "done ...\n";
        cout << "steps " << step << "\n";
        cout << "residual " << residual <<  " < " << eps << "\n";
    }

    dash::finalize();
    return 0;
}
