#include <DUGKSAmr.H>
#include <DUGKSField.H>

using amrex::FabType;

#ifdef AMREX_USE_EB
#include <AMReX_EBCellFlag.H>
using amrex::EBCellFlagFab;
#endif

void DUGKSAmr::ComputeDivU(Real time_in)
{
    //to be done
}
void DUGKSAmr::ComputeDrag()
{
    //to be done
}
//
Real DUGKSAmr::ComputeKineticEnergy () const
{
    //to be done
    return -1;
}

//The following codes are directly copied from incflo
void DUGKSAmr::ComputeStrainrateAtLevel (int lev,
                                         MultiFab* strainrate,
                                         MultiFab* vel,
                                         Geometry& lev_geom,
                                         Real time, int nghost)
{

#ifdef AMREX_USE_EB
    auto const& fact = EBFactory(lev);
    auto const& flags = fact.getMultiEBCellFlagFab();
#endif

    AMREX_D_TERM(Real idx = 1.0 / lev_geom.CellSize(0);,
                 Real idy = 1.0 / lev_geom.CellSize(1);,
                 Real idz = 1.0 / lev_geom.CellSize(2););

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*strainrate,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.growntilebox(nghost);
        Array4<Real> const& sr_arr = strainrate->array(mfi);
        Array4<Real const> const& vel_arr = vel->const_array(mfi);
#ifdef AMREX_USE_EB
        auto const& flag_fab = flags[mfi];
        auto typ = flag_fab.getType(bx);
        if (typ == FabType::covered)
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                sr_arr(i,j,k) = 0.0;
            });
        }
        else if (typ == FabType::singlevalued)
        {
            auto const& flag_arr = flag_fab.const_array();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                sr_arr(i,j,k) = incflo_strainrate_eb(i,j,k,AMREX_D_DECL(idx,idy,idz),
                                                     vel_arr,flag_arr(i,j,k));
            });
        }
        else
#endif
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                sr_arr(i,j,k) = incflo_strainrate(i,j,k,AMREX_D_DECL(idx,idy,idz),vel_arr);
            });
        }
    }
}

//
#if (AMREX_SPACEDIM == 2)
void DUGKSAmr::ComputeVorticity (int lev, Real t, MultiFab& vort, MultiFab const& vel)
{
    BL_PROFILE("DUGKSAmr::ComputeVorticity");
    const Real idx = Geom(lev).InvCellSize(0);
    const Real idy = Geom(lev).InvCellSize(1);

#ifdef AMREX_USE_EB
    const auto& fact = EBFactory(lev);
    const auto& flags_mf = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for(MFIter mfi(vel, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();
        Array4<Real const> const& ccvel_fab = vel.const_array(mfi);
        Array4<Real> const& vort_fab = vort.array(mfi);

#ifdef AMREX_USE_EB
        const EBCellFlagFab& flags = flags_mf[mfi];
        auto typ = flags.getType(bx);
        if (typ == FabType::covered)
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                vort_fab(i,j,k) = 0.0;
            });
        }
        else if (typ == FabType::singlevalued)
        {
            const auto& flag_fab = flags.const_array();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                constexpr Real c0 = -1.5;
                constexpr Real c1 = 2.0;
                constexpr Real c2 = -0.5;

                if (flag_fab(i,j,k).isCovered())
                {
                    vort_fab(i,j,k) = 0.0;
                }
                else
                {
                    Real vx, uy;
                    // Need to check if there are covered cells in neighbours --
                    // -- if so, use one-sided difference computation (but still quadratic)
                    if (!flag_fab(i,j,k).isConnected( 1,0,0))
                    {
                        // Covered cell to the right, go fish left
                        vx = - (c0 * ccvel_fab(i  ,j,k,1)
                              + c1 * ccvel_fab(i-1,j,k,1)
                              + c2 * ccvel_fab(i-2,j,k,1)) * idx;
                    }
                    else if (!flag_fab(i,j,k).isConnected(-1,0,0))
                    {
                        // Covered cell to the left, go fish right
                        vx = (c0 * ccvel_fab(i  ,j,k,1)
                            + c1 * ccvel_fab(i+1,j,k,1)
                            + c2 * ccvel_fab(i+2,j,k,1)) * idx;
                    }
                    else
                    {
                        // No covered cells right or left, use standard stencil
                        vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                    }
                    // Do the same in y-direction
                    if (!flag_fab(i,j,k).isConnected(0, 1,0))
                    {
                        uy = - (c0 * ccvel_fab(i,j  ,k,0)
                              + c1 * ccvel_fab(i,j-1,k,0)
                              + c2 * ccvel_fab(i,j-2,k,0)) * idy;
                    }
                    else if (!flag_fab(i,j,k).isConnected(0,-1,0))
                    {
                        uy = (c0 * ccvel_fab(i,j  ,k,0)
                            + c1 * ccvel_fab(i,j+1,k,0)
                            + c2 * ccvel_fab(i,j+2,k,0)) * idy;
                    }
                    else
                    {
                        uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                    }
                    vort_fab(i,j,k) = vx-uy;
                }
            });
        }
        else
#endif
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                Real uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                vort_fab(i,j,k) = vx-uy;
            });
        }
    }
}
    
#elif (AMREX_SPACEDIM == 3)
void DUGKSAmr::ComputeVorticity (int lev, Real t, MultiFab& vort, MultiFab const& vel)
{
    BL_PROFILE("incflo::ComputeVorticity");
    const Real idx = Geom(lev).InvCellSize(0);
    const Real idy = Geom(lev).InvCellSize(1);
    const Real idz = Geom(lev).InvCellSize(2);

#ifdef AMREX_USE_EB
    const auto& fact = EBFactory(lev);
    const auto& flags_mf = fact.getMultiEBCellFlagFab();
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for(MFIter mfi(vel, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.tilebox();
        Array4<Real const> const& ccvel_fab = vel.const_array(mfi);
        Array4<Real> const& vort_fab = vort.array(mfi);

#ifdef AMREX_USE_EB
        const EBCellFlagFab& flags = flags_mf[mfi];
        auto typ = flags.getType(bx);
        if (typ == FabType::covered)
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                vort_fab(i,j,k) = 0.0;
            });
        }
        else if (typ == FabType::singlevalued)
        {
            const auto& flag_fab = flags.const_array();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                constexpr Real c0 = -1.5;
                constexpr Real c1 = 2.0;
                constexpr Real c2 = -0.5;

                if (flag_fab(i,j,k).isCovered())
                {
                    vort_fab(i,j,k) = 0.0;
                }
                else
                {
                    Real vx, wx, uy, wy, uz, vz;
                    // Need to check if there are covered cells in neighbours --
                    // -- if so, use one-sided difference computation (but still quadratic)
                    if (!flag_fab(i,j,k).isConnected( 1,0,0))
                    {
                        // Covered cell to the right, go fish left
                        vx = - (c0 * ccvel_fab(i  ,j,k,1)
                              + c1 * ccvel_fab(i-1,j,k,1)
                              + c2 * ccvel_fab(i-2,j,k,1)) * idx;
                        wx = - (c0 * ccvel_fab(i  ,j,k,2)
                              + c1 * ccvel_fab(i-1,j,k,2)
                              + c2 * ccvel_fab(i-2,j,k,2)) * idx;
                    }
                    else if (!flag_fab(i,j,k).isConnected(-1,0,0))
                    {
                        // Covered cell to the left, go fish right
                        vx = (c0 * ccvel_fab(i  ,j,k,1)
                            + c1 * ccvel_fab(i+1,j,k,1)
                            + c2 * ccvel_fab(i+2,j,k,1)) * idx;
                        wx = (c0 * ccvel_fab(i  ,j,k,2)
                            + c1 * ccvel_fab(i+1,j,k,2)
                            + c2 * ccvel_fab(i+2,j,k,2)) * idx;
                    }
                    else
                    {
                        // No covered cells right or left, use standard stencil
                        vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                        wx = 0.5 * (ccvel_fab(i+1,j,k,2) - ccvel_fab(i-1,j,k,2)) * idx;
                    }
                    // Do the same in y-direction
                    if (!flag_fab(i,j,k).isConnected(0, 1,0))
                    {
                        uy = - (c0 * ccvel_fab(i,j  ,k,0)
                              + c1 * ccvel_fab(i,j-1,k,0)
                              + c2 * ccvel_fab(i,j-2,k,0)) * idy;
                        wy = - (c0 * ccvel_fab(i,j  ,k,2)
                              + c1 * ccvel_fab(i,j-1,k,2)
                              + c2 * ccvel_fab(i,j-2,k,2)) * idy;
                    }
                    else if (!flag_fab(i,j,k).isConnected(0,-1,0))
                    {
                        uy = (c0 * ccvel_fab(i,j  ,k,0)
                            + c1 * ccvel_fab(i,j+1,k,0)
                            + c2 * ccvel_fab(i,j+2,k,0)) * idy;
                        wy = (c0 * ccvel_fab(i,j  ,k,2)
                            + c1 * ccvel_fab(i,j+1,k,2)
                            + c2 * ccvel_fab(i,j+2,k,2)) * idy;
                    }
                    else
                    {
                        uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                        wy = 0.5 * (ccvel_fab(i,j+1,k,2) - ccvel_fab(i,j-1,k,2)) * idy;
                    }
                    // Do the same in z-direction
                    if (!flag_fab(i,j,k).isConnected(0,0, 1))
                    {
                        uz = - (c0 * ccvel_fab(i,j,k  ,0)
                              + c1 * ccvel_fab(i,j,k-1,0)
                              + c2 * ccvel_fab(i,j,k-2,0)) * idz;
                        vz = - (c0 * ccvel_fab(i,j,k  ,1)
                              + c1 * ccvel_fab(i,j,k-1,1)
                              + c2 * ccvel_fab(i,j,k-2,1)) * idz;
                    }
                    else if (!flag_fab(i,j,k).isConnected(0,0,-1))
                    {
                        uz = (c0 * ccvel_fab(i,j,k  ,0)
                            + c1 * ccvel_fab(i,j,k+1,0)
                            + c2 * ccvel_fab(i,j,k+2,0)) * idz;
                        vz = (c0 * ccvel_fab(i,j,k  ,1)
                            + c1 * ccvel_fab(i,j,k+1,1)
                            + c2 * ccvel_fab(i,j,k+2,1)) * idz;
                    }
                    else
                    {
                        uz = 0.5 * (ccvel_fab(i,j,k+1,0) - ccvel_fab(i,j,k-1,0)) * idz;
                        vz = 0.5 * (ccvel_fab(i,j,k+1,1) - ccvel_fab(i,j,k-1,1)) * idz;
                    }
                    vort_fab(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
                }
            });
        }
        else
#endif
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real vx = 0.5 * (ccvel_fab(i+1,j,k,1) - ccvel_fab(i-1,j,k,1)) * idx;
                Real wx = 0.5 * (ccvel_fab(i+1,j,k,2) - ccvel_fab(i-1,j,k,2)) * idx;

                Real uy = 0.5 * (ccvel_fab(i,j+1,k,0) - ccvel_fab(i,j-1,k,0)) * idy;
                Real wy = 0.5 * (ccvel_fab(i,j+1,k,2) - ccvel_fab(i,j-1,k,2)) * idy;

                Real uz = 0.5 * (ccvel_fab(i,j,k+1,0) - ccvel_fab(i,j,k-1,0)) * idz;
                Real vz = 0.5 * (ccvel_fab(i,j,k+1,1) - ccvel_fab(i,j,k-1,1)) * idz;

                vort_fab(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
            });
        }
    }
}
#endif

void DUGKSAmr::ReverseVelocityField (MultiFab &VelFab, const Real *dx, 
                                     const RealBox &rb, int const token)
{
    double di = 0;
    double dj = 0;

    if(token == -1)
    {
        di = 0.5,dj = 0.5;
    }
    else if(token == Ark::X)
    {
        di = 0,dj = 0.5;
    }
    else if(token == Ark::Y)
    {
        di = 0.5,dj = 0;
    }
    else {
        amrex::Abort("meaningless face id");
    }

    const double real_lo_x = rb.lo(Ark::X);
    const double real_hi_x = rb.hi(Ark::X);

    const double real_lo_y = rb.lo(Ark::Y);
    const double real_hi_y = rb.hi(Ark::Y);

    for(amrex::MFIter mfi(VelFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        const Real dlx = dx[0];
        const Real dly = dx[1];
        //
        Array4<Real> const &vel = VelFab.array(mfi);

        amrex::For(bx,[=](int i,int j,int k)
        {
            const double x = (i + di)*dlx;
            const double y = (j + dj)*dly;
            
            vel(i,j,k,Ark::X) = 
            -U0()*Ark::pi*sin(Ark::pi*x/ChLength())*cos(Ark::pi*y/ChLength());
            vel(i,j,k,Ark::Y) =  
            U0()*Ark::PI*cos(Ark::PI*x/ChLength())*sin(Ark::PI*y/ChLength());
    
            if(token == Ark::X){
            if(x == real_lo_x or x == real_hi_x)
            {
                vel(i,j,k,Ark::X) = 0;
                vel(i,j,k,Ark::Y) = 0; 
            }}
            //
            if(token == Ark::Y){
            if(y == real_lo_y or y == real_hi_y)
            {
                vel(i,j,k,Ark::X) = 0;
                vel(i,j,k,Ark::Y) = 0; 
            }}
        });
    }
}
void DUGKSAmr::SmoothedVelocityField (MultiFab &VelFab, const Real *dx, 
                                      const Real Time, int const token)
{
    const Real dlx = dx[0];
    const Real dly = dx[1];

    Real di = 0,dj = 0;

    if(token == -1)
    {
        di = 0.5,dj = 0.5;
    }
    else if(token == Ark::X)
    {
        di = 0,dj = 0.5;
    }
    else if(token == Ark::Y)
    {
        di = 0.5,dj = 0;
    }
    else {
        amrex::Abort("meaningless token");
    }

    for(amrex::MFIter mfi(VelFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &vel = VelFab.array(mfi);

        amrex::For(bx,[=] (int i,int j,int k)
        {
            const double x = (i + di)*dlx;
            const double y = (j + dj)*dly;
            
            vel(i,j,k,Ark::X) = 
            -U0()*sin(4.0*Ark::PI*x/ChLength())*sin(4.0*Ark::PI*y/ChLength())
            *cos(Ark::PI*Time*U0()/ChLength());

            vel(i,j,k,Ark::Y) = 
            -U0()*cos(4.0*Ark::PI*x/ChLength())*cos(4.0*Ark::PI*y/ChLength())
            *cos(Ark::PI*Time*U0()/ChLength());
        });
    }
}