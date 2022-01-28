#include <DUGKSAmr.H>

double DUGKSAmr::L2NormSubtract(const MultiFab& xFabA, const MultiFab& xFabB, 
							  int const comp)
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
double DUGKSAmr::L2NormSubtractZalesakDisk (const MultiFab& xFabA, const MultiFab& xFabB)
{
    AMREX_ASSERT(xFabA.boxArray() == xFabB.boxArray());
    AMREX_ASSERT(xFabA.DistributionMap() == xFabB.DistributionMap());

    std::unique_ptr<MultiFab> xFabAmend;

    xFabAmend.reset(new MultiFab(xFabA.boxArray(), xFabA.DistributionMap(), 1, 0, MFInfo()));

    MultiFab::Copy(*xFabAmend, xFabB, 0, 0, 1, zero_ghost());

    for(amrex::MFIter mfi(*xFabAmend,true);mfi.isValid();++mfi)
    {
        const Box &bx = mfi.tilebox();

        Array4<Real> const &phi = (*xFabAmend).array(mfi);

        amrex::ParallelFor(bx,[=](int i,int j,int k)
        {
            phi(i,j,k) = phi(i,j,k) > 0.5*(PhiL() + PhiV()) ? PhiL() : PhiV();
        });
    }

    std::unique_ptr<MultiFab> diff_xFab;

    diff_xFab.reset(new MultiFab(xFabA.boxArray(), xFabA.DistributionMap(), 1, 0, MFInfo()));



    MultiFab::LinComb(*diff_xFab, 
                      1.0, xFabA, 0, 
                      -1.0, *xFabAmend, 0, 
                      0, 1, 0);

    return diff_xFab->norm2();
}

double DUGKSAmr::L2NormVelocity(int lev)
{
	auto &ld = Data(lev);

	Real norm_diff[AMREX_SPACEDIM];
	Real norm_vel[ AMREX_SPACEDIM];

	const MultiFab &velFabNew = ld.msv[Ark::VEL];
	MultiFab &velFabOld = ld.velocity_old;

	norm_diff[Ark::X] = L2NormSubtract(velFabNew, velFabOld, Ark::X);
	norm_diff[Ark::Y] = L2NormSubtract(velFabNew, velFabOld, Ark::Y);

	norm_vel[ Ark::X] = velFabNew.norm2(Ark::X);
	norm_vel[ Ark::Y] = velFabNew.norm2(Ark::Y);

#if(AMREX_SPACEDIM == 3)	
	norm_diff[Ark::Z] = L2NormSubtract(velFabNew, velFabOld, Ark::Z);
	norm_vel[ Ark::Z] = velFabNew.norm2(Ark::Z);

    MultiFab::Copy(velFabOld, velFabNew, 0, 0, AMREX_SPACEDIM, 0);

    WriteL2NormVelocity (lev, 
        std::sqrt(Ark::Square(norm_diff[Ark::X], norm_diff[Ark::Y], norm_diff[Ark::Z])),
        std::sqrt(Ark::Square(norm_vel[Ark::X], norm_vel[Ark::Y], norm_vel[Ark::Z])));

	return 	std::sqrt(
			Ark::Square(norm_diff[Ark::X], norm_diff[Ark::Y], norm_diff[Ark::Z]) 
			/Ark::Square(norm_vel[Ark::X], norm_vel[Ark::Y], norm_vel[Ark::Z]));
#else

    MultiFab::Copy(velFabOld, velFabNew, 0, 0, AMREX_SPACEDIM, 0);
//
    WriteL2NormVelocity (lev, std::sqrt(Ark::Square(norm_diff[Ark::X], norm_diff[Ark::Y])),
                         std::sqrt(Ark::Square(norm_vel[Ark::X], norm_vel[Ark::Y])));
    return std::sqrt(
			Ark::Square(norm_diff[Ark::X], norm_diff[Ark::Y]) 
			/Ark::Square(norm_vel[Ark::X], norm_vel[Ark::Y]));

#endif
}

double DUGKSAmr::L2NormPhase(int lev)
{
	auto &ld = Data(lev);

	Real norm_diff   = 0;
	Real norm_phase  = 0;

	const MultiFab &phiFabNew = ld.msv[Ark::PHI];
	MultiFab &phiFabOld = ld.phase_old;

	norm_diff = L2NormSubtract(phiFabNew, phiFabOld, 0);

	norm_phase = phiFabNew.norm2();

    MultiFab::Copy(phiFabOld, phiFabNew, 0, 0, 1, 0);

	return norm_diff/(norm_phase);
}

double DUGKSAmr::L2NormDensity(int lev)
{
	auto &ld = Data(lev);

	Real norm_diff     = 0;
	Real norm_density  = 0;

	const MultiFab &rhoFabNew = ld.msv[Ark::RHO];
	MultiFab &rhoFabOld = ld.density_old;

	norm_diff = L2NormSubtract(rhoFabNew, rhoFabOld, 0);

	norm_density = rhoFabNew.norm2();

    MultiFab::Copy(rhoFabOld, rhoFabNew, 0, 0, 1, 0);

	return norm_diff/(norm_density + Ark::tiny);
}
double DUGKSAmr::SumPhase(int lev)
{
    return Data(lev).msv[Ark::PHI].sum();
}

Vector<double> DUGKSAmr::MultiLevelSum()
{
    Vector<double> cur_sum;

    for(int lev = finest_level; lev >=0; --lev)
    {
        cur_sum.push_back(SumPhase(lev));
    }
    return cur_sum;
}
Vector<double> DUGKSAmr::MultiLevelResidual()
{
    Vector<double> cur_residual;

    for(int lev = finest_level; lev >=0; --lev)
    {
#ifdef ARK_USE_MOMENTUM
        cur_residual.push_back(L2NormVelocity(lev));
#elif defined ARK_USE_PHASE
        cur_residual.push_back(L2NormPhase(lev));
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
    //     ReinitializePhase(Data(lev).phase_old);
    // }

    //* * * * * * * * * * * * * * * * Unsteady TGV * * * * * * * * * * * * * * * * 
    // for(int lev = finest_level; lev >=0; --lev) {
    //     FillLevelUnsteadyTaylorGreenVortex(lev, Data(lev).velocity_prev, Data(lev).phase_prev);
    // }

    //* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
    //* * * * * * * * * * * * * * * * Scalar L2-norm * * * * * * * * * * * * * * *
    for(int lev = finest_level; lev >=0; --lev)
    {
        // - General
        Real const norm_diff = L2NormSubtract(Data(lev).density_old, Data(lev).msv[Ark::PHI]);
        
        // - Zalesak Disk
        // Real const norm_diff = L2NormSubtract(Data(lev).density_old, Data(lev).phase_old);
        
        Real const norm_phase = Data(lev).density_old.norm2();        
        cur_L2_norm.push_back(norm_diff/norm_phase);
    }
    //* * * * * * * * * * * * * * * * Vector L2-norm * * * * * * * * * * * * * * *
    // for(int lev = finest_level; lev >=0; --lev) {   
    //     Real const norm_diff_x = 
    //     L2NormSubtract(Data(lev).velocity_prev, Data(lev).msv[Ark::VEL], Ark::X);

    //     Real const norm_diff_y = 
    //     L2NormSubtract(Data(lev).velocity_prev, Data(lev).msv[Ark::VEL], Ark::Y);

    //     Real const norm_vel_x = Data(lev).msv[Ark::VEL].norm2(Ark::X);
    //     Real const norm_vel_y = Data(lev).msv[Ark::VEL].norm2(Ark::Y);

    //     cur_L2_norm.push_back(
    //     std::sqrt(Ark::Square(norm_diff_x, norm_diff_y)/Ark::Square(norm_vel_x, norm_vel_y))
    //     );          
    // }

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
    const double phiMid = (PhiL() + PhiV())/2;
    const double deltaPhi = (PhiL() - PhiV())/8.0;
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
    const double phiMid = (PhiL() + PhiV())/2;
    const double deltaPhi = (PhiL() - PhiV())/8.0;
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
    const int boundY = (CenterPF(Ark::Y)/ChLength() + (nstep*2.0)/end_step)*ChNum();

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

    const double highPosition = InterfacePositionHighest(PhaseFab, prob_lo_y, dx, midX);
    const double lowPosition  = InterfacePositionLowest (PhaseFab, prob_lo_y, dx, midX);

    return {highPosition, lowPosition};
}