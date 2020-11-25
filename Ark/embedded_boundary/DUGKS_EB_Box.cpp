#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <DUGKSAmr.H>

using namespace amrex;

/****************************************************************************
 * Function to create a simple rectangular box with EB walls.               *
 *                                                                          *
 ****************************************************************************/
void DUGKSAmr::make_eb_box()
{
    // Get box information from inputs file
    ParmParse pp("box");

    if(geom[0].isAllPeriodic())
    {
        make_eb_regular();
    }
    else
    {
        std::cout << "IN MAKE_EB_BOX NOT PERIODIC " << std::endl;
        exit(0);
        /************************************************************************
         *                                                                      *
         * Define Box geometry:                                                 *
         *        -> box.{Lo,Hi} vector storing box lo/hi                       *
         *        -> box.offset  vector storing box offset                      *
         * NOTE: walls are placed _outside_ domain for periodic directions.     *
         *                                                                      *
         ************************************************************************/

        Vector<Real> boxLo(3), boxHi(3);
        Real offset = 1.0e-15;

        for(int i = 0; i < 3; i++)
        {
            boxLo[i] = geom[0].ProbLo(i);
            boxHi[i] = geom[0].ProbHi(i);
        }

        pp.queryarr("Lo", boxLo, 0, 3);
        pp.queryarr("Hi", boxHi, 0, 3);

        pp.query("offset", offset);

        Real xlo = boxLo[0] + offset;
        Real xhi = boxHi[0] - offset;

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if(geom[0].isPeriodic(0))
        {
            xlo = 2.0 * geom[0].ProbLo(0) - geom[0].ProbHi(0);
            xhi = 2.0 * geom[0].ProbHi(0) - geom[0].ProbLo(0);
        }

        Real ylo = boxLo[1] + offset;
        Real yhi = boxHi[1] - offset;

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if(geom[0].isPeriodic(1))
        {
            ylo = 2.0 * geom[0].ProbLo(1) - geom[0].ProbHi(1);
            yhi = 2.0 * geom[0].ProbHi(1) - geom[0].ProbLo(1);
        }

        Real zlo = boxLo[2] + offset;
        Real zhi = boxHi[2] - offset;

        // This ensures that the walls won't even touch the ghost cells. By
        // putting them one domain width away
        if(geom[0].isPeriodic(2))
        {
            zlo = 2.0 * geom[0].ProbLo(2) - geom[0].ProbHi(2);
            zhi = 2.0 * geom[0].ProbHi(2) - geom[0].ProbLo(2);
        }

        Array<Real, AMREX_SPACEDIM> point_lox{AMREX_D_DECL(xlo, 0.0, 0.0)};
        Array<Real, AMREX_SPACEDIM> normal_lox{AMREX_D_DECL(-1.0, 0.0, 0.0)};
        Array<Real, AMREX_SPACEDIM> point_hix{AMREX_D_DECL(xhi, 0.0, 0.0)};
        Array<Real, AMREX_SPACEDIM> normal_hix{AMREX_D_DECL(1.0, 0.0, 0.0)};

        Array<Real, AMREX_SPACEDIM> point_loy{AMREX_D_DECL(0.0, ylo, 0.0)};
        Array<Real, AMREX_SPACEDIM> normal_loy{AMREX_D_DECL(0.0, -1.0, 0.0)};
        Array<Real, AMREX_SPACEDIM> point_hiy{AMREX_D_DECL(0.0, yhi, 0.0)};
        Array<Real, AMREX_SPACEDIM> normal_hiy{AMREX_D_DECL(0.0, 1.0, 0.0)};

        Array<Real, AMREX_SPACEDIM> point_loz{AMREX_D_DECL(0.0, 0.0, zlo)};
        Array<Real, AMREX_SPACEDIM> normal_loz{AMREX_D_DECL(0.0, 0.0, -1.0)};
        Array<Real, AMREX_SPACEDIM> point_hiz{AMREX_D_DECL(0.0, 0.0, zhi)};
        Array<Real, AMREX_SPACEDIM> normal_hiz{AMREX_D_DECL(0.0, 0.0, 1.0)};

        EB2::PlaneIF plane_lox(point_lox, normal_lox);
        EB2::PlaneIF plane_hix(point_hix, normal_hix);

        EB2::PlaneIF plane_loy(point_loy, normal_loy);
        EB2::PlaneIF plane_hiy(point_hiy, normal_hiy);

        EB2::PlaneIF plane_loz(point_loz, normal_loz);
        EB2::PlaneIF plane_hiz(point_hiz, normal_hiz);

        // Generate GeometryShop
        auto gshop = EB2::makeShop(EB2::makeUnion(plane_lox, plane_hix,
                                                  plane_loy, plane_hiy,
                                                  plane_loz, plane_hiz));

        // Build index space
        int max_level_here = 0;
        int max_coarsening_level = 100;
        EB2::Build(gshop, geom.back(), max_level_here, max_level_here + max_coarsening_level);
    }
}
