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

double DUGKSAmr::L2NormVelocity(int lev)
{
	auto &ld = Data(lev);

	Real norm_diff[AMREX_SPACEDIM];
	Real norm_vel[ AMREX_SPACEDIM];

	const MultiFab &velFabNew = ld.msv[Ark::VEL];
	MultiFab &velFabOld = ld.velocity_old;

	norm_diff[Ark::dir_x] = L2NormSubtract(velFabNew, velFabOld, Ark::dir_x);
	norm_diff[Ark::dir_y] = L2NormSubtract(velFabNew, velFabOld, Ark::dir_y);

	norm_vel[ Ark::dir_x] = velFabNew.norm2(Ark::dir_x);
	norm_vel[ Ark::dir_y] = velFabNew.norm2(Ark::dir_y);

#if(AMREX_SPACEDIM == 3)	
	norm_diff[Ark::dir_z] = L2NormSubtract(velFabNew, velFabOld, Ark::dir_z);
	norm_vel[ Ark::dir_z] = velFabNew.norm2(Ark::dir_z);

    MultiFab::Copy(velFabOld, velFabNew, 0, 0, AMREX_SPACEDIM, 0);

	return 	amrex::Math::sqrt(
			Ark::Square(norm_diff[Ark::dir_x], norm_diff[Ark::dir_y], norm_diff[Ark::dir_z]) 
			/Ark::Square(norm_vel[Ark::dir_x], norm_vel[Ark::dir_y], norm_vel[Ark::dir_z]));
#else

    MultiFab::Copy(velFabOld, velFabNew, 0, 0, AMREX_SPACEDIM, 0);
    return amrex::Math::sqrt(
			Ark::Square(norm_diff[Ark::dir_x], norm_diff[Ark::dir_y], 0) 
			/Ark::Square(norm_vel[Ark::dir_x], norm_vel[Ark::dir_y], 0));

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
void DUGKSAmr::FillLevelPhaseAnalytical(MultiFab &phase, const Real *dl, LV token)
{
	for(amrex::MFIter mfi(phase, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const Box& bx = mfi.tilebox();
        //
        Array4<Real> const &phi  = phase.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
        	const double x = dl[Ark::dir_x]*(i+0.5);
            const double y = dl[Ark::dir_y]*(j+0.5);
#if(AMREX_SPACEDIM == 3)
            const double z = dl[Ark::dir_z]*(k+0.5);

		phi(i,j,k) = 
		0.5*(PhiL() + PhiV()) + 0.5*static_cast<int>(token)*(PhiL() - PhiV())*
        tanh((
        	amrex::Math::sqrt(
        	Ark::Square(x-CenterPF(Ark::dir_x), y-CenterPF(Ark::dir_y), z-CenterPF(Ark::dir_z))
        	) - RadiusPF())*2.0/WidthI()
        );

#else
        phi(i,j,k) = 
        0.5*(PhiL() + PhiV()) + 0.5*static_cast<int>(token)*(PhiL() - PhiV())*
        tanh((std::hypot(x-CenterPF(Ark::dir_x), y-CenterPF(Ark::dir_y)) - RadiusPF()
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
        	const double x = dl[Ark::dir_x]*(i+0.5);
            const double y = dl[Ark::dir_y]*(j+0.5);

			Real slot_xL = CenterPF(Ark::dir_x) - 0.1*RadiusPF();
			Real slot_xR = CenterPF(Ark::dir_x) + 0.1*RadiusPF();
#if(AMREX_SPACEDIM == 3)
			Real distance = amrex::Math::sqrt(
        	Ark::Square(x-CenterPF(Ark::dir_x), y-CenterPF(Ark::dir_y), z-CenterPF(Ark::dir_z))
    		);
#else
			Real distance = std::hypot(x-CenterPF(Ark::dir_x), y-CenterPF(Ark::dir_y));
#endif
			if (distance > RadiusPF()) {
				phi(i,j,k) = PhiV();
			}
			else if (x > slot_xL and x < slot_xR and y < CenterPF(Ark::dir_y)) {
				phi(i,j,k) = PhiV();
			}
			else {
				phi(i,j,k) = PhiL();
			}
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

            vel(i,j,k,Ark::dir_x) = 
            -U0()*cos(2.0*pi*x)*sin(2.0*pi*y)*exp(-8.0*pi*pi*Nu0()*time)/(2.0*pi);

            vel(i,j,k,Ark::dir_y) = 
            U0()*sin(2.0*pi*x)*cos(2.0*pi*y)*exp(-8.0*pi*pi*Nu0()*time)/(2.0*pi);

            pres(i,j,k) = -0.25*U0()*U0()*(cos(4.0*pi*x) + cos(4.0*pi*y))
                                  *exp(-16.0*pi*pi*Nu0()*time)/(4.0*pi*pi);
        });
    }
}
Vector<double> DUGKSAmr::MultiLevelResidual()
{
    Vector<double> cur_residual;

    for(int lev = finest_level; lev >=0; --lev)
    {
        cur_residual.push_back(L2NormPhase(lev));
    }
    return cur_residual;
}
Vector<double> DUGKSAmr::MultiLevelL2Norm()
{
    Vector<double> cur_L2_norm;

    for(int lev = finest_level; lev >=0; --lev)
    {
        FillLevelPhaseAnalytical(Data(lev).density_old, Geom(lev).CellSize(), LV::drop);
    }
    for(int lev = finest_level; lev >=0; --lev)
    {
        Real const norm_diff  = L2NormSubtract(Data(lev).density_old, Data(lev).msv[Ark::PHI]);
        // Real const norm_phase = Data(lev).msv[Ark::PHI].norm2();
        Real const norm_phase = Data(lev).density_old.norm2();        
        cur_L2_norm.push_back(norm_diff/norm_phase);
        // cur_L2_norm.push_back(norm_diff);
    }
    return cur_L2_norm;
}