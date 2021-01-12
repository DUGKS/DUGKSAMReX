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

		phi(i,j,k) = 
		0.5*(PhiL() + PhiV()) + 0.5*static_cast<int>(token)*(PhiL() - PhiV())*
        tanh((
        	std::sqrt(
        	Ark::Square(x-CenterPF(Ark::X), y-CenterPF(Ark::Y), z-CenterPF(Ark::Z))
        	) - RadiusPF())*2.0/WidthI()
        );

#else
        phi(i,j,k) = 
        0.5*(PhiL() + PhiV()) + 0.5*static_cast<int>(token)*(PhiL() - PhiV())*
        tanh((std::hypot(x-CenterPF(Ark::X), y-CenterPF(Ark::Y)) - RadiusPF()
        	 )*2.0/WidthI()
        	
        );			 
#endif
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
void DUGKSAmr::reinitialize(MultiFab &phase)
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
//unsteady taylor-green vortex
void DUGKSAmr::FillLevelVelocityUnsteadyTGV(MultiFab &Velocity, MultiFab &Pressure,
											const Real *dl,const Real time)
{
	using namespace Ark;

	for(amrex::MFIter mfi(Velocity, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        //
        Array4<Real> const &vel   = Velocity.array(mfi);
        Array4<Real> const &pres  = Pressure.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double x = dl[0]*(i+0.5);
            const double y = dl[1]*(j+0.5);

            vel(i,j,k,Ark::X) = 
            -U0()*cos(2.0*pi*x)*sin(2.0*pi*y)*exp(-8.0*pi*pi*Nu0()*time)/(2.0*pi);

            vel(i,j,k,Ark::Y) = 
            U0()*sin(2.0*pi*x)*cos(2.0*pi*y)*exp(-8.0*pi*pi*Nu0()*time)/(2.0*pi);

            pres(i,j,k) = -0.25*U0()*U0()*(cos(4.0*pi*x) + cos(4.0*pi*y))
                                  *exp(-16.0*pi*pi*Nu0()*time)/(4.0*pi*pi);
        });
    }
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
        cur_residual.push_back(L2NormVelocity(lev));
    }
    return cur_residual;
}
Vector<double> DUGKSAmr::MultiLevelL2Norm()
{
    Vector<double> cur_L2_norm;

    for(int lev = finest_level; lev >=0; --lev)
    {
        FillLevelPhaseSphere(Data(lev).density_old, Geom(lev).CellSize(), LV::drop);
    }
    // for(int lev = finest_level; lev >=0; --lev)
    // {
    //     MultiFab::Copy(Data(lev).density_old, Data(lev).msv[Ark::PHI],0,0,1,zero_ghost());
    //     reinitialize(Data(lev).density_old);
    // }
    for(int lev = finest_level; lev >=0; --lev)
    {
        Real const norm_diff = L2NormSubtract(Data(lev).density_old, Data(lev).msv[Ark::PHI]);
        // Real const norm_diff = L2NormSubtract(Data(lev).density_old, Data(lev).phase_old);
        
        Real const norm_phase = Data(lev).density_old.norm2();        
        cur_L2_norm.push_back(norm_diff/norm_phase);
    }
    return cur_L2_norm;
}
double DUGKSAmr::InterfacePosition (const MultiFab& PhaseFab, const double probDomain_lo_y,
                                    const Real* dx, int const indexI, int const indexK)
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
            if(i == indexI and k == indexK)
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
Vector<double> DUGKSAmr::InterfacePosition(int lev)
{
    const Real *dx = Geom(lev).CellSize();
    const Box &domain = Geom(lev).Domain();
    const RealBox &rb = Geom(lev).ProbDomain();
    const auto &ld = Data(lev);

    const MultiFab &PhaseFab = ld.msv[Ark::PHI];

#if(AMREX_SPACEDIM == 3)   
    const double bubble = InterfacePosition(PhaseFab, rb.lo(Ark::Y), dx, 0, 0);
    const double saddle = InterfacePosition(PhaseFab, rb.lo(Ark::Y), dx, 
                                            0, domain.length(Ark::Z)/2);
    const double spike  = InterfacePosition(PhaseFab, rb.lo(Ark::Y), dx, 
                                            domain.length(Ark::X)/2, domain.length(Ark::Z)/2);
    return {bubble,saddle,spike};
#else
    const double bubble = InterfacePosition(PhaseFab, rb.lo(Ark::Y), dx, 0, 0);
    const double spike  = InterfacePosition(PhaseFab, rb.lo(Ark::Y), dx,  
                                            domain.length(Ark::X)/2, 0);

    return {bubble,spike};
#endif 
}