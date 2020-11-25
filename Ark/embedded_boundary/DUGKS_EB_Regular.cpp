#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <DUGKS_EB_IF.H> 
#include <DUGKSAmr.H>

using namespace amrex;

void DUGKSAmr::make_eb_regular()
{
    EB2::AllRegularIF my_regular;
    auto gshop = EB2::makeShop(my_regular);
    EB2::Build(gshop, geom.back(), 0, 100);
}
