#include <DUGKSAmr.H>

double DUGKSAmr::L2NormSubtract(const MultiFab& xFabA, const MultiFab& xFabB, int const comp)
{
    AMREX_ASSERT(xFabA.boxArray() == xFabB.boxArray());
    AMREX_ASSERT(xFabA.DistributionMap() == xFabB.DistributionMap());

    std::unique_ptr<MultiFab> diff_xFab;

    diff_xFab.reset(new MultiFab(xFabA.boxArray(), xFabA.DistributionMap(), 1, 0, MFInfo()));

    MultiFab::LinComb(*diff_xFab,
                       1.0, xFabA, comp,
                      -1.0, xFabB, comp,
                      0, 1, 0);

    return diff_xFab->norm2();
}
double DUGKSAmr::L2NormScalar(MultiFab &scalarFabOld, const MultiFab &scalarFabNew, const bool updateOld)
{
    const Real norm_diff   = L2NormSubtract(scalarFabOld, scalarFabNew);
    const Real norm_scalar = scalarFabNew.norm2();

    if(updateOld) {
        MultiFab::Copy(scalarFabOld, scalarFabNew, 0, 0, 1, 0)
    }
    return norm_diff / (norm_scalar + Ark::tiny);
}
double DUGKSAmr::ResidualPhase(int lev)
{
	auto &ld = Data(lev);

	MultiFab &phiFabOld = ld.phase_old;
	const MultiFab &phiFabNew = ld.msv[Ark::PHI];

    return L2NormScalar(phiFabOld, phiFabNew, true);
}

double DUGKSAmr::ResidualDensity(int lev)
{
	auto &ld = Data(lev);

	MultiFab &rhoFabOld = ld.density_old;
	const MultiFab &rhoFabNew = ld.msv[Ark::RHO];

	return L2NormScalar(rhoFabOld, rhoFabNew, true);
}
double DUGKAmr::L2NormVector(MultiFab &vectorFabOld, const MultiFab &vectorFabNew, bool updateOld)
{
    Real norm_diff[AMREX_SPACEDIM];
	Real norm_vector[ AMREX_SPACEDIM];

	norm_diff[Ark::X] = L2NormSubtract(vectorFabOld, vectorFabNew, Ark::X);
	norm_diff[Ark::Y] = L2NormSubtract(vectorFabOld, vectorFabNew, Ark::Y);

	norm_vector[Ark::X]  = vectorFabNew.norm2(Ark::X);
	norm_vector[ Ark::Y] = vectorFabNew.norm2(Ark::Y);

#if(AMREX_SPACEDIM == 3)
	norm_diff[Ark::Z]   = L2NormSubtract(vectorFabOld, vectorFabNew, Ark::Z);
	norm_vector[Ark::Z] = vectorFabNew.norm2(Ark::Z);

    if(updateOld) {
        MultiFab::Copy(vectorFabOld, vectorFabNew, 0, 0, AMREX_SPACEDIM, 0);
    }

    WriteL2NormVelocity (lev,
        std::sqrt(Ark::Square(norm_diff[Ark::X], norm_diff[Ark::Y], norm_diff[Ark::Z])),
        std::sqrt(Ark::Square(norm_vel[Ark::X], norm_vel[Ark::Y], norm_vel[Ark::Z])));

	return 	std::sqrt( Ark::Square(norm_diff[Ark::X], norm_diff[Ark::Y], norm_diff[Ark::Z])
                      /Ark::Square(norm_vector[Ark::X], norm_vector[Ark::Y], norm_vector[Ark::Z]));
#else
    if(updateOld) {
        MultiFab::Copy(vectorFabOld, vectorFabNew, 0, 0, AMREX_SPACEDIM, 0);
    }

    return std::sqrt( Ark::Square(norm_diff[Ark::X], norm_diff[Ark::Y])
			         /Ark::Square(norm_vector[Ark::X], norm_vector[Ark::Y]));
#endif
}
double DUGKSAmr::ResidualVelocity(int lev)
{
	auto &ld = Data(lev);

	const MultiFab &velFabNew = ld.msv[Ark::VEL];
	MultiFab &velFabOld = ld.velocity_old;

    return L2NormVector(velFabOld, velFabNew, true);
    // WriteL2NormVelocity (lev,
    //     std::sqrt(Ark::Square(norm_diff[Ark::X], norm_diff[Ark::Y], norm_diff[Ark::Z])),
    //     std::sqrt(Ark::Square(norm_vel[Ark::X], norm_vel[Ark::Y], norm_vel[Ark::Z])));
}

Vector<double> DUGKSAmr::MultiLevelSumPhase()
{
    Vector<double> cur_sum;

    for(int lev = finest_level; lev >= 0; --lev)
    {
        cur_sum.push_back(Data(lev).msv[Ark::PHI].sum());
    }
    return cur_sum;
}
Vector<double> DUGKSAmr::MultiLevelSumDensity()
{
    Vector<double> cur_sum;

    for(int lev = finest_level; lev >= 0; --lev)
    {
        cur_sum.push_back(Data(lev).msv[Ark::RHO].sum());
    }
    return cur_sum;
}
Vector<double> DUGKSAmr::MultiLevelResidual()
{
    Vector<double> cur_residual;

    for(int lev = finest_level; lev >=0; --lev)
    {
#ifdef ARK_USE_MOMENTUM
        cur_residual.push_back(ResidualVelocity(lev));
        // cur_residual.push_back(ResidualDensity(lev));
#elif defined ARK_USE_PHASE
        cur_residual.push_back(ResidualPhase(lev));
#endif
    }
    return cur_residual;
}
Vector<double> DUGKSAmr::MultiLevelL2Norm()
{
    Vector<double> cur_L2_norm;

    //* * * * * * * * * * * * * * * * Drop * * * * * * * * * * * * * * * * * * * *
    // do nothing

    //* * * * * * * * * * * * * * * * Zalesak Disk * * * * * * * * * * * * * * * *
    // for(int lev = finest_level; lev >=0; --lev)
    // {
    //     MultiFab::Copy(Data(lev).phase_old, Data(lev).msv[Ark::PHI],0,0,1,zero_ghost());
    //     StaircasePhase(Data(lev).phase_old);
    // }

    // * * * * * * * * * * * * * * * * Unsteady TGV * * * * * * * * * * * * * * * *
    for(int lev = finest_level; lev >= 0; --lev) {
        FillLevelUnsteadyTaylorGreenVortex(lev, Data(lev).velocity_prev, Data(lev).phase_prev);
    }

    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    // //* * * * * * * * * * * * Scalar L2-norm Phase Field  * * * * * * * * * * * * *
    // for(int lev = finest_level; lev >=0; --lev)
    // {
    //     // - General
    //     Real const norm_diff = L2NormSubtract(Data(lev).density_old, Data(lev).msv[Ark::PHI]);

    //     // - Zalesak Disk
    //     // Real const norm_diff = L2NormSubtract(Data(lev).density_old, Data(lev).phase_old);

    //     Real const norm_phase = Data(lev).density_old.norm2();
    //     cur_L2_norm.push_back(norm_diff/norm_phase);
    // }

    //* * * * * * * * * * * * Scalar L2-norm Density Field  * * * * * * * * * * * * *
    // for(int lev = finest_level; lev >=0; --lev)
    // {
    //     // - General
    //     Real const norm_diff = L2NormSubtract(Data(lev).phase_prev, Data(lev).msv[Ark::P]);

    //     Real const norm_scalar = Data(lev).density_old.norm2();
    //     cur_L2_norm.push_back(norm_diff/norm_scalar);
    // }

    //* * * * * * * * * * * * * * * * Vector L2-norm Velocity Field * * * * * * * * * * * * * * *
    for(int lev = finest_level; lev >=0; --lev) {
        Real const norm_diff_x =
        L2NormSubtract(Data(lev).velocity_prev, Data(lev).msv[Ark::VEL], Ark::X);

        Real const norm_diff_y =
        L2NormSubtract(Data(lev).velocity_prev, Data(lev).msv[Ark::VEL], Ark::Y);

        // Real const norm_vel_x = Data(lev).msv[Ark::VEL].norm2(Ark::X);
        // Real const norm_vel_y = Data(lev).msv[Ark::VEL].norm2(Ark::Y);
        Real const norm_vel_x = Data(lev).velocity_prev.norm2(Ark::X);
        Real const norm_vel_y = Data(lev).velocity_prev.norm2(Ark::Y);

        cur_L2_norm.push_back(
        std::sqrt(Ark::Square(norm_diff_x, norm_diff_y)/Ark::Square(norm_vel_x, norm_vel_y))
        );
    }

    return cur_L2_norm;
}
Vector<double> DUGKSAmr::L2NormUnsteadyTaylorGreenVortex()
{
    Vector<double> cur_L2_norm;

    // analytical velocity : velocity_prev    analytical pressure : phase_prev
    for(int lev = finest_level; lev >= 0; --lev) {
        FillLevelUnsteadyTaylorGreenVortex(lev, Data(lev).velocity_prev, Data(lev).phase_prev);
    }

    for(int lev = finest_level; lev >=0; --lev) {
        const MultiFab &velFabAnalytic = Data(lev).velocity_prev;
        MultiFab &velFabNumeric = Data(lev).msv[Ark::VEL];

        cur_L2_norm.push_back(L2NormVector(velFabNumeric, velFabAnalytic));
    }
    return cur_L2_norm;
}

Vector<double> DUGKSAmr::L2NormZalesakDisk()
{
    Vector<double> cur_L2_norm;

    // analytical solution -> density_old    computational solution -> phase_old
    for(int lev = finest_level; lev >= 0; --lev)
    {
        MultiFab::Copy(Data(lev).phase_old, Data(lev).msv[Ark::PHI],0,0,1,zero_ghost());
        StaircasePhase(Data(lev).phase_old);
    }

    for(int lev = finest_level; lev >= 0; --lev)
    {
        const MultiFab &phaseFabAnalytic = Data(lev).density_old;
        MultiFab &phaseFabNumeric  = Data(lev).phase_old;

        cur_L2_norm.push_back(L2NormScalar(phaseFabNumeric, phaseFabAnalytic));
    }

    return cur_L2_norm;
}
Vector<double> DUGKSAmr::L2NormPhase()
{
    Vector<double> cur_L2_norm;
    // analytical solution -> density_old
    for(int lev = finest_level; lev >=0; --lev)
    {
        const MultiFab &phaseFabAnalytic = Data(lev).msv[Ark::PHI];
        MultiFab &phaseFabNumeric = Data(lev).density_old;
        cur_L2_norm.push_back(L2NormScalar(phaseFabNumeric, phaseFabAnalytic));
    }
    return cur_L2_norm;
}

Vector<double> DUGKSAmr::L2NormPressure()
{
    Vector<double> cur_L2_norm;

    // analytical solution ->
    // for(int lev = finest_level; lev >=0; --lev)
    // {
    //     // - General
    //     Real const norm_diff = L2NormSubtract(Data(lev).phase_prev, Data(lev).msv[Ark::P]);

    //     Real const norm_scalar = Data(lev).density_old.norm2();
    //     cur_L2_norm.push_back(norm_diff/norm_scalar);
    // }

    //* * * * * * * * * * * * * * * * Vector L2-norm Velocity Field * * * * * * * * * * * * * * *
    for(int lev = finest_level; lev >=0; --lev) {
        Real const norm_diff_x =
        L2NormSubtract(Data(lev).velocity_prev, Data(lev).msv[Ark::VEL], Ark::X);

        Real const norm_diff_y =
        L2NormSubtract(Data(lev).velocity_prev, Data(lev).msv[Ark::VEL], Ark::Y);

        // Real const norm_vel_x = Data(lev).msv[Ark::VEL].norm2(Ark::X);
        // Real const norm_vel_y = Data(lev).msv[Ark::VEL].norm2(Ark::Y);
        Real const norm_vel_x = Data(lev).velocity_prev.norm2(Ark::X);
        Real const norm_vel_y = Data(lev).velocity_prev.norm2(Ark::Y);

        cur_L2_norm.push_back(
        std::sqrt(Ark::Square(norm_diff_x, norm_diff_y)/Ark::Square(norm_vel_x, norm_vel_y))
        );
    }

    return cur_L2_norm;
}
double DUGKSAmr::InterfacePosition (const MultiFab& PhaseFab, const double probDomain_lo_y,
                                    const Real* dx,
                                    int const indexI, int const indexJ, int const indexK)
{
    const double phiMid = (PhiL() + PhiV())/2;
    double position = -1;
    double deltaPhi = PhiL() - PhiV();
    int min_proc = 0;
//
    for(amrex::MFIter mfi(PhaseFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        //
        Array4<Real const> const &phi  = PhaseFab.const_array(mfi);

        amrex::ParallelFor(bx,
        [=, &position, &deltaPhi] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if(j > indexJ and i == indexI and k == indexK)
            {
                if (std::fabs(phi(i,j,k) - phiMid) < deltaPhi)
                {
                    deltaPhi = std::fabs(phi(i,j,k) - phiMid);
                    position = ((j + 0.5)*dx[Ark::Y] + probDomain_lo_y)/ChLength();
                }
            }
        });
    }

    std::vector<double> dphiv = amrex::ParallelDescriptor::Gather(deltaPhi,
                                amrex::ParallelDescriptor::IOProcessorNumber());
    if (amrex::ParallelDescriptor::IOProcessor()) {
        min_proc = std::min_element(dphiv.begin(), dphiv.end()) - dphiv.begin();
    }
//
    amrex::ParallelDescriptor::Bcast(&min_proc,1,amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::Bcast(&position,1,min_proc);

    return position;
}

double DUGKSAmr::InterfacePositionHighest (const MultiFab& PhaseFab, const double probDomain_lo_y,
                                           const Real* dx, int const index)
{
#ifdef ARK_USE_PHASE
    const double phiMid = (PhiL() + PhiV())/2;
    const double deltaPhi = (PhiL() - PhiV())/8.0;
#else
    const double phiMid = (RhoL() + RhoV())/2;
    const double deltaPhi = (RhoL() - RhoV())/8.0;
#endif
    double position = -1;
//
    for(amrex::MFIter mfi(PhaseFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        double temp_position = -1;
        //
        Array4<Real const> const &phi  = PhaseFab.const_array(mfi);

        amrex::ParallelFor(bx,
        [=, &position, &temp_position] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i == index && std::fabs(phi(i,j,k) - phiMid) < deltaPhi) {

               temp_position = ((j + 0.5)*dx[Ark::Y] + probDomain_lo_y)/ChLength();

               if(temp_position > position) {
                    position = temp_position;
               }
            }
        });
    }

    amrex::ParallelDescriptor::ReduceRealMax(position);

    return position;
}

double DUGKSAmr::InterfacePositionLowest (const MultiFab& PhaseFab, const double probDomain_lo_y,
                                          const Real* dx, int const index)
{
#ifdef ARK_USE_PHASE
    const double phiMid = (PhiL() + PhiV())/2;
    const double deltaPhi = (PhiL() - PhiV())/8.0;
#else
    const double phiMid = (RhoL() + RhoV())/2;
    const double deltaPhi = (RhoL() - RhoV())/8.0;
#endif
    double position = std::numeric_limits<double>::max();
//
    for(amrex::MFIter mfi(PhaseFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        double temp_position = std::numeric_limits<double>::max();
        //
        Array4<Real const> const &phi  = PhaseFab.const_array(mfi);

        amrex::ParallelFor(bx,
        [=, &position, &temp_position] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i == index and std::fabs(phi(i,j,k) - phiMid) < deltaPhi)
            {
                temp_position = ((j + 0.5)*dx[Ark::Y] + probDomain_lo_y)/ChLength();

                if(temp_position < position) {
                    position = temp_position;
               }
            }
        });
    }

    amrex::ParallelDescriptor::ReduceRealMin(position);

    return position;
}

// double DUGKSAmr::InterfacePositionLowest (const MultiFab& PhaseFab, const double probDomain_lo_y,
//                                           const Real* dx)
// {
//     const double phiMid = (PhiL() + PhiV())/2;
//     const double deltaPhi = (PhiL() - PhiV())/10.0;
//     double position = std::numeric_limits<double>::max();
// //
//     for(amrex::MFIter mfi(PhaseFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
//     {
//         const Box& bx = mfi.tilebox();
//         double temp_position = std::numeric_limits<double>::max();
//         //
//         Array4<Real const> const &phi = PhaseFab.const_array(mfi);

//         amrex::ParallelFor(bx,
//         [=, &position, &temp_position] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//         {
//             if (std::fabs(phi(i,j,k) - phiMid) < deltaPhi)
//             {
//                 temp_position = ((j + 0.5)*dx[Ark::Y] + probDomain_lo_y)/ChLength();

//                 if(temp_position < position) {
//                     position = temp_position;
//                }
//             }
//         });
//     }

//     amrex::ParallelDescriptor::ReduceRealMin(position);

//     return position;
// }

Vector<double> DUGKSAmr::InterfacePositionRayleighTaylor(int lev)
{
    const Real *dx = Geom(lev).CellSize();
    const Box &domain = Geom(lev).Domain();
    const RealBox &rb = Geom(lev).ProbDomain();
    const auto &ld = Data(lev);

    const double prob_lo_y = rb.lo(Ark::Y);
    const int midX = domain.length(Ark::X)/2;
    const int boundY = 0;

    const MultiFab &PhaseFab = ld.msv[Ark::PHI];

#if(AMREX_SPACEDIM == 3)
    const int midZ = domain.length(Ark::Z)/2;
    const double bubble = InterfacePosition(PhaseFab, prob_lo_y, dx, 0, boundY, 0);
    const double saddle = InterfacePosition(PhaseFab, prob_lo_y, dx, 0, boundY, midZ);
    const double spike  = InterfacePosition(PhaseFab, prob_lo_y, dx, midX, boundY, midZ);

    return {bubble,saddle,spike};
#else
    const double bubble = InterfacePosition(PhaseFab, prob_lo_y, dx, 0, boundY, 0);
    const double spike  = InterfacePosition(PhaseFab, prob_lo_y, dx, midX, boundY, 0);

    return {bubble,spike};
#endif
}

Vector<double> DUGKSAmr::InterfacePositionTaylorBubble(int lev)
{
    const Real *dx = Geom(lev).CellSize();
    const Box &domain = Geom(lev).Domain();
    const RealBox &rb = Geom(lev).ProbDomain();
    const auto &ld = Data(lev);

    const double prob_lo_y = rb.lo(Ark::Y);
    const int midX = domain.length(Ark::X)/2;
    const int boundY = static_cast<int>((CenterPF(Ark::Y)/ChLength() + (nstep*2.0)/end_step)*ChNum());

    const MultiFab &PhaseFab = ld.msv[Ark::PHI];

#if(AMREX_SPACEDIM == 3)
    const int midZ = domain.length(Ark::Z)/2;
    const double spike  = InterfacePosition(PhaseFab, prob_lo_y, dx, midX, boundY, midZ);
#else
    const double spike  = InterfacePosition(PhaseFab, prob_lo_y, dx, midX, boundY, 0);
#endif

    return {spike};
}

Vector<double> DUGKSAmr::InterfacePositionRisingBubble(int lev)
{
    const Real *dx = Geom(lev).CellSize();
    const Box &domain = Geom(lev).Domain();
    const RealBox &rb = Geom(lev).ProbDomain();
    const auto &ld = Data(lev);

    const int midX = domain.length(Ark::X)/2;
    const double prob_lo_y = rb.lo(Ark::Y);

    const MultiFab &PhaseFab = ld.msv[Ark::PHI];
    const MultiFab &DensityFab = ld.msv[Ark::RHO];

    const double highPosition = InterfacePositionHighest(DensityFab, prob_lo_y, dx, midX);
    const double lowPosition  = InterfacePositionLowest (DensityFab, prob_lo_y, dx, midX);

    return {highPosition, lowPosition};
}

Vector<double> DUGKSAmr::InterfacePositionLiquidBridge(int lev)
{
    const Real *dx = Geom(lev).CellSize();
    const Box &domain = Geom(lev).Domain();
    const RealBox &rb = Geom(lev).ProbDomain();
    const auto &ld = Data(lev);

    const int midX = domain.length(Ark::X)/2;
    const double prob_lo_y = rb.lo(Ark::Y);

    const MultiFab &PhaseFab = ld.msv[Ark::PHI];
    const MultiFab &DensityFab = ld.msv[Ark::RHO];

    const double highPosition = InterfacePositionHighest(DensityFab, prob_lo_y, dx, midX);
    const double lowPosition  = InterfacePositionLowest (DensityFab, prob_lo_y, dx, midX);

    const double bridge = (highPosition - lowPosition)*ChLength()/RadiusPF()/2.0;

    return {bridge};
}