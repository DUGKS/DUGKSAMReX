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
//

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

double DUGKSAmr::PhaseProfile (double const d, LV token)
{
    return 0.5*(PhiL()+PhiV() + static_cast<int>(token)*(PhiL()-PhiV())*tanh(d*2.0/WidthI()));
}
void DUGKSAmr::Nullify_Source (int const lev)
{
    auto &ld = Data(lev);

    // ld.msv[Ark::GradPHI].setVal(0);
    ld.msv[Ark::GradRHO].setVal(0);
    ld.msv[Ark::FORCE].setVal(0);
}
void DUGKSAmr::BoundPhase(MultiFab &PhaseFab)
{
    const Real bound = 5e-6;
    // const Real bound = 0;

    for(amrex::MFIter mfi(PhaseFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        //
        Array4<Real> const &phi  = PhaseFab.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if(phi(i,j,k) < bound) {
                phi(i,j,k) = 0;
            }
            else if(phi(i,j,k) > 1 - bound) {
                phi(i,j,k) = 1;
            }
        });
    }
}

void DUGKSAmr::BoundGradient(MultiFab &GradPhaseFab)
{
    for(amrex::MFIter mfi(GradPhaseFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        //
        Array4<Real> const &GradPhi  = GradPhaseFab.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if(std::fabs(GradPhi(i,j,k,Ark::X)) < 1e-12) {
                GradPhi(i,j,k,Ark::X) = 0;
            }
            if(std::fabs(GradPhi(i,j,k,Ark::Y)) < 1e-12) {
                GradPhi(i,j,k,Ark::Y) = 0;
            }
        });
    }
}
//Sudden Reverse 2D
void DUGKSAmr::ElongationVelocityField (MultiFab &VelFab, const Real *dx, const RealBox &rb, int const token)
{
    const double multiplier = (nstep <=0 ? 1.0 : -1.0);

    double di = 0;
    double dj = 0;

    if(token == -1) {
        di = 0.5,dj = 0.5;
    }
    else if(token == Ark::X) {
        di = 0,dj = 0.5;
    }
    else if(token == Ark::Y) {
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

        amrex::ParallelFor(bx,[=](int i,int j,int k)
        {
            const double x = (i + di)*dlx + real_lo_x;
            const double y = (j + dj)*dly + real_lo_y;

            const double xT = x/ChLength()*Ark::PI;
            const double yT = y/ChLength()*Ark::PI;
            
            vel(i,j,k,Ark::X) =  U0()*Ark::PI*sin(xT)*cos(yT)*multiplier;
            vel(i,j,k,Ark::Y) = -U0()*Ark::PI*cos(xT)*sin(yT)*multiplier;

            if(token == Ark::X) {
                if(x == real_lo_x or x == real_hi_x)
                {
                    vel(i,j,k,Ark::X) = 0;
                    vel(i,j,k,Ark::Y) = 0; 
                }
            }
            //
            if(token == Ark::Y) {
                if(y == real_lo_y or y == real_hi_y)
                {
                    vel(i,j,k,Ark::X) = 0;
                    vel(i,j,k,Ark::Y) = 0; 
                }
            }
        });
    }
}

//Sudden Reverse 3D
void DUGKSAmr::ElongationVelocityField3D  (MultiFab &VelFab, const Real *dx, const RealBox &rb, int const token)
{
#if(AMREX_SPACEDIM == 3)

    // const double multiplier = (nstep <=0 ? 1.0 : -1.0);
    const double tT = cur_time*U0()/(2*ChLength());

    double di = 0;
    double dj = 0;
    double dk = 0;

    if(token == -1) {
        di = 0.5, dj = 0.5, dk = 0.5;
    }
    else if(token == Ark::X) {
        di = 0, dj = 0.5, dk = 0.5;
    }
    else if(token == Ark::Y) {
        di = 0.5, dj = 0, dk = 0.5;
    }
    else if(token == Ark::Z) {
        di = 0.5, dj = 0.5, dk = 0.0;
    }
    else {
        amrex::Abort("meaningless face id");
    }

    const Real dlx = dx[Ark::X];
    const Real dly = dx[Ark::Y];
    const Real dlz = dx[Ark::Z];

    const double real_lo_x = rb.lo(Ark::X);
    const double real_hi_x = rb.hi(Ark::X);

    const double real_lo_y = rb.lo(Ark::Y);
    const double real_hi_y = rb.hi(Ark::Y);

    const double real_lo_z = rb.lo(Ark::Z);
    const double real_hi_z = rb.hi(Ark::Z);

    for(amrex::MFIter mfi(VelFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &vel = VelFab.array(mfi);

        amrex::ParallelFor(bx,[=](int i,int j,int k)
        {
            const double x = (i + di)*dlx + real_lo_x;
            const double y = (j + dj)*dly + real_lo_y;
            const double z = (k + dk)*dlz + real_lo_z;

            const double xT = (x/ChLength()-0.5)*Ark::PI;
            const double yT = (y/ChLength()-0.5)*Ark::PI;
            const double zT = (z/ChLength()-0.5)*Ark::PI;
            
            vel(i,j,k,Ark::X) = Ark::PI*U0()*cos(xT)*(sin(yT)-sin(zT))*cos(Ark::PI*tT);
            vel(i,j,k,Ark::Y) = Ark::PI*U0()*cos(yT)*(sin(zT)-sin(xT))*cos(Ark::PI*tT);
            vel(i,j,k,Ark::Z) = Ark::PI*U0()*cos(zT)*(sin(xT)-sin(yT))*cos(Ark::PI*tT);

            if(token == Ark::X) {
                if(x == real_lo_x or x == real_hi_x) {
                    vel(i,j,k,Ark::X) = 0;
                    vel(i,j,k,Ark::Y) = 0; 
                    vel(i,j,k,Ark::Z) = 0; 
                }
            }
            //
            if(token == Ark::Y) {
                if(y == real_lo_y or y == real_hi_y)
                {
                    vel(i,j,k,Ark::X) = 0;
                    vel(i,j,k,Ark::Y) = 0; 
                    vel(i,j,k,Ark::Z) = 0; 
                }
            }

            if(token == Ark::Z) {
                if(z == real_lo_z or z == real_hi_z)
                {
                    vel(i,j,k,Ark::X) = 0;
                    vel(i,j,k,Ark::Y) = 0; 
                    vel(i,j,k,Ark::Z) = 0; 
                }
            }
        });
    }
#endif
}

//smoothed
void DUGKSAmr::DeformationVelocityField (MultiFab &VelFab, const Real *dx, 
                                         const Real time, int const token)
{
    const Real dlx = dx[0];
    const Real dly = dx[1];

    double di = 0;
    double dj = 0;

    if(token == -1) {
        di = 0.5,dj = 0.5;
    }
    else if(token == Ark::X) {
        di = 0,dj = 0.5;
    }
    else if(token == Ark::Y) {
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

        amrex::ParallelFor(bx,[=] (int i,int j,int k)
        {
            const double x = (i + di)*dlx/ChLength();
            const double y = (j + dj)*dly/ChLength();
            
            vel(i,j,k,Ark::X) = 
            -U0()*sin(4.0*Ark::PI*x)*sin(4.0*Ark::PI*y)*cos(Ark::PI*time*U0()/ChLength());

            vel(i,j,k,Ark::Y) = 
            -U0()*cos(4.0*Ark::PI*x)*cos(4.0*Ark::PI*y)*cos(Ark::PI*time*U0()/ChLength());
        });
    }
}

// Trailing tail
void DUGKSAmr::ReversedVortexVelocityField (MultiFab &VelFab, const Real *dx, 
                                            const Real time, int const token)
{
    using Ark::PI;

    const Real dlx = dx[0];
    const Real dly = dx[1];

    double di = 0;
    double dj = 0;

    if(token == -1) {
        di = 0.5, dj = 0.5;
    }
    else if(token == Ark::X) {
        di = 0, dj = 0.5;
    }
    else if(token == Ark::Y) {
        di = 0.5, dj = 0;
    }
    else {
        amrex::Abort("meaningless token");
    }

    for(amrex::MFIter mfi(VelFab);mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        //
        Array4<Real> const &vel = VelFab.array(mfi);

        amrex::ParallelFor(bx,[=] (int i,int j,int k)
        {
            const double x = (i + di)*dlx/ChLength();
            const double y = (j + dj)*dly/ChLength();
            
            vel(i,j,k,Ark::X) = 
            U0()*sin(PI*x)*sin(PI*x)*sin(2.0*PI*y)*cos(PI*time*U0()/ChLength()/4.0);

            vel(i,j,k,Ark::Y) = 
            -U0()*sin(PI*y)*sin(PI*y)*sin(2.0*PI*x)*cos(PI*time*U0()/ChLength()/4.0);
        });
    }
}
void DUGKSAmr::FillLevelUnsteadyTaylorGreenVortex(int const lev, MultiFab &VelFab, MultiFab &PresFab)
{
    using Ark::PI;
    const Real *dl = Geom(lev).CellSize();

    for(amrex::MFIter mfi(VelFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();

        Array4<Real> const &vel  = VelFab.array(mfi);
        Array4<Real> const &pres = PresFab.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double x = dl[Ark::X]*(i+0.5)/ChLength();
            const double y = dl[Ark::Y]*(j+0.5)/ChLength();

            vel(i,j,k,Ark::X) = 
            -U0()*cos(2.0*PI*x)*sin(2.0*PI*y)*exp(-8.0*PI*PI*Nu0()*cur_time)/(2.0*PI);

            vel(i,j,k,Ark::Y) = 
            U0()*sin(2.0*PI*x)*cos(2.0*PI*y)*exp(-8.0*PI*PI*Nu0()*cur_time)/(2.0*PI);

            pres(i,j,k) = -0.25*U0()*U0()*(cos(4.0*PI*x) + cos(4.0*PI*y))
                        *exp(-16.0*PI*PI*Nu0()*cur_time)/(4.0*PI*PI);
        });
    }
}
void DUGKSAmr::FillLevelForceTaylorGreenVortex (int const lev, MultiFab &VelFab, MultiFab &PresFab, 
                                                MultiFab &ForceFab, int const token)
{
    using Ark::PI;
    const Real *dl = Geom(lev).CellSize();
    const Real multiplier = 8.0*PI*PI*Nu0()*U0()/ChLength()/ChLength();

    double di = 0;
    double dj = 0;

    if(token == -1) {
        di = 0.5, dj = 0.5;
    }
    else if(token == Ark::X) {
        di = 0, dj = 0.5;
    }
    else if(token == Ark::Y) {
        di = 0.5, dj = 0;
    }
    else {
        amrex::Abort("meaningless token in function [FillLevelForceTaylorGreenVortex]");
    }

    for(amrex::MFIter mfi(VelFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();

        Array4<Real> const &vel  = VelFab.array(mfi);
        Array4<Real> const &pres = PresFab.array(mfi);
        Array4<Real> const &force = ForceFab.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double x = dl[Ark::X]*(i + di)/ChLength();
            const double y = dl[Ark::Y]*(j + dj)/ChLength();

            vel(i,j,k,Ark::X) = U0()*sin(2.0*PI*x)*sin(2.0*PI*y);
            vel(i,j,k,Ark::Y) = U0()*cos(2.0*PI*x)*cos(2.0*PI*y);

            force(i,j,k,Ark::X) = multiplier*sin(2.*PI*x)*sin(2.*PI*y);
            force(i,j,k,Ark::Y) = multiplier*cos(2.*PI*x)*cos(2.*PI*y);

            pres(i,j,k) = 0.25*U0()*U0()*(cos(4.0*PI*x)-cos(4.0*PI*y));            
        });
    }
}
void DUGKSAmr::FillLevelPhaseSphere(MultiFab &phase, const Real *dl, LV token)
{
    for(amrex::MFIter mfi(phase, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        //
        Array4<Real> const &phi  = phase.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double x = dl[Ark::X]*(i+0.5);
            const double y = dl[Ark::Y]*(j+0.5);
#if(AMREX_SPACEDIM == 3)
            const double z = dl[Ark::Z]*(k+0.5);
            const double distance = 
            Ark::Hypot(x-CenterPF(Ark::X), y-CenterPF(Ark::Y), z-CenterPF(Ark::Z)) - RadiusPF();
#else
            const double distance = std::hypot(x-CenterPF(Ark::X), y-CenterPF(Ark::Y)) - RadiusPF();           
#endif
            phi(i,j,k) = PhaseProfile(distance, token);
        });
    }
}
void DUGKSAmr::FillLevelPhaseZalesakDisk(MultiFab &phase, const Real *dl)
{
    for(amrex::MFIter mfi(phase, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        //
        Array4<Real> const &phi  = phase.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double x = dl[Ark::X]*(i+0.5);
            const double y = dl[Ark::Y]*(j+0.5);

            Real slot_xL = CenterPF(Ark::X) - 0.1*RadiusPF();
            Real slot_xR = CenterPF(Ark::X) + 0.1*RadiusPF();
#if(AMREX_SPACEDIM == 3)
            const double z = dl[Ark::Z]*(k+0.5);
            Real distance = std::sqrt(
            Ark::Square(x-CenterPF(Ark::X), y-CenterPF(Ark::Y), z-CenterPF(Ark::Z))
            );
#else
            Real distance = std::hypot(x-CenterPF(Ark::X), y-CenterPF(Ark::Y));
#endif
            if (distance > RadiusPF()) {
                phi(i,j,k) = PhiV();
            }
            else if (x > slot_xL and x < slot_xR and y < CenterPF(Ark::Y)) {
                phi(i,j,k) = PhiV();
            }
            else {
                phi(i,j,k) = PhiL();
            }
        });
    }
}
void DUGKSAmr::ReinitializePhase(MultiFab &phase)
{
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for(amrex::MFIter mfi(phase, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        //
        Array4<Real> const &phi  = phase.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phi(i,j,k) = phi(i,j,k) > 0.5*(PhiL() + PhiV()) ? PhiL() : PhiV();
        }); 
    }   
}