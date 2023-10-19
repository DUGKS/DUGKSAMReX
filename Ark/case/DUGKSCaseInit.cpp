#include <DUGKSAmr.H>

void DUGKSAmr::CaseInitFluid(int const lev)
{
    case_init_macroscopic_variables(lev);

    switch(CaseIndex())
    {
        case 0:
        {
            UniformFlow(lev);
            break;
        }
        case 1:
        {
            UnsteadyTaylorGreenVortex(lev);
            break;
        }
        case 5:
        {
            ForceTaylorGreenVortex(lev);
            break;
        }
        case 219:
        {
            ReversedSingleVortex(lev);
            break;
        }
        case 220:
        {
            QuiescentDrop(lev);
            break;
        }
        case 221:
        {
            TransDrop(lev);
            break;
        }
        case 222:
        {
            ZalesakDisk(lev);
            break;
        }
        case 223:
        {
            ElongatedInterface(lev);
            break;
        }
        case 224:
        {
            SmoothDeformation(lev);
            break;
        }
        case 225:
        {
            LayeredPoiseuille(lev);
            break;
        }
        case 226:
        {
            SpinodalDecomposition(lev);
            break;
        }
        case 227:
        {
            RayleighTaylorInstability(lev);
            break;
        }
        case 228:
        {
            FlatInterface(lev);
            break;
        }
        case 229:
        {
            DropletImpact(lev);
            break;
        }
        case 230:
        {
            TaylorBubble(lev);
            break;
        }
        case 231:
        {
            RisingBubble(lev);
            break;
        }
        case 232:
        {
            FallingDrop(lev);
            break;
        }
        case 233:
        {
            DropletCoalescence(lev);
            break;
        }
        case -666:
        {
            BasicCase(lev);
            break;
        }
    	default :
    	{
    		std::ostringstream ss;
            ARK_PRINT_ERROR_MESSAGE(ss,__FILE__,__LINE__,__func__);
            ss  <<"Invalid Case ID\n";
            amrex::Abort(ss.str());
    	}
    }

    case_init_mesoscopic_equilibrium(lev);

    post_case_init(lev);
}
void DUGKSAmr::case_init_macroscopic_variables(int const lev)
{
	auto& ld = Data(lev);

    ld.phase_old.setVal(0);
    ld.density_old.setVal(0);
    ld.velocity_old.setVal(0);

    ld.phase_prev.setVal(0);
    ld.velocity_prev.setVal(0);

    for(auto &mf : ld.msv){
        mf.setVal(0.0);
    }

    // ld.msv[Ark::PHI].setVal(Rho0());
    // ld.msv[Ark::RHO].setVal(Rho0());
    // ld.msv[Ark::VEL].setVal(U0());
    // ld.msv[Ark::P].setVal(P0());
    // ld.msv[Ark::T].setVal(T0());

    for(auto &mf : ld.msv_xF) {
        mf.setVal(0.0);
    }

    for(auto &mf : ld.msv_yF) {
        mf.setVal(0.0);
    }
#if(AMREX_SPACEDIM == 3)
    for(auto &mf : ld.msv_zF) {
        mf.setVal(0.0);
    }
#endif

#ifdef ARK_USE_PHASE
    ld.hDVDF[Ark::TAU ].setVal(TauPhase());
	ld.hDVDFxF[Ark::TAUF].setVal(TauPhase());
	ld.hDVDFyF[Ark::TAUF].setVal(TauPhase());
#if(AMREX_SPACEDIM == 3)
    ld.hDVDFzF[Ark::TAUF].setVal(TauPhase());
#endif
    for(int n = Ark::BPx;n < Ark::DVDFNUM;++n) {
		ld.hDVDF[n].setVal(0);
	}
#endif

#ifdef ARK_USE_MOMENTUM
    ld.fDVDF[Ark::TAU ].setVal(TauMomentum());
	ld.fDVDFxF[Ark::TAUF].setVal(TauMomentum());
	ld.fDVDFyF[Ark::TAUF].setVal(TauMomentum());
#if(AMREX_SPACEDIM == 3)
    ld.fDVDFzF[Ark::TAUF].setVal(TauMomentum());
#endif
    for(int n = Ark::BPx;n < Ark::DVDFNUM;++n) {
		ld.fDVDF[n].setVal(0);
	}
#endif

#ifdef ARK_USE_ENERGY
    ld.gDVDF[Ark::TAU ].setVal(TauEnergy());
	ld.gDVDFxF[Ark::TAUF].setVal(TauEnergy());
	ld.gDVDFyF[Ark::TAUF].setVal(TauEnergy());
#if(AMREX_SPACEDIM == 3)
    ld.gDVDFzF[Ark::TAUF].setVal(TauEnergy());
#endif
    for(int n = Ark::BPx;n < Ark::DVDFNUM;++n) {
		ld.gDVDF[n].setVal(0);
	}
#endif
}

void DUGKSAmr::case_init_mesoscopic_equilibrium(int const lev)
{
	auto& ld = *m_leveldata[lev];

	// Force initialized to 0
    Update_DVDF_Eq(lev);

#ifdef ARK_USE_PHASE
	ld.hDVDF[Ark::WAVE].ParallelCopy(ld.hDVDF[Ark::EQ], 0, 0, Ark::DVM_PHASE,
		IntVect{0}, IntVect{0}, amrex::Periodicity::NonPeriodic());
#endif

#ifdef ARK_USE_MOMENTUM
	ld.fDVDF[Ark::WAVE].ParallelCopy(ld.fDVDF[Ark::EQ], 0, 0, Ark::DVM_MOMENTUM,
		IntVect{0}, IntVect{0}, amrex::Periodicity::NonPeriodic());
#endif
}
void DUGKSAmr::post_case_init(int const lev)
{
    auto& ld = Data(lev);

    //The following would result in a slight increase in L-2 norm error;
    MultiFab::Copy(ld.phase_prev, ld.msv[Ark::PHI], 0, 0, 1, unit_ghost());
    MultiFab::Copy(ld.velocity_prev, ld.msv[Ark::VEL], 0, 0, AMREX_SPACEDIM, unit_ghost());

#ifdef ARK_USE_PHASE
    FillPatchPhase(lev, cur_time, ld.msv[Ark::PHI]);

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for(amrex::MFIter mfi(ld.msv[Ark::PHI], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        Array4<Real const> const &phi      = ld.msv[Ark::PHI].const_array(mfi);
        Array4<Real      > const &grad_phi = ld.msv[Ark::GradPHI].array(mfi);

        GradFunc::gradient_operator_isotropy(bx, Geom(lev).CellSize(), phi, grad_phi);
    }
#else
    FillPatchDensity(lev, cur_time, ld.msv[Ark::RHO]);

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    for(amrex::MFIter mfi(ld.msv[Ark::RHO], amrex::TilingIfNotGPU());mfi.isValid();++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        Array4<Real const> const &rho      = ld.msv[Ark::RHO].const_array(mfi);
        Array4<Real      > const &grad_rho = ld.msv[Ark::GradRHO].array(mfi);

        GradFunc::gradient_operator_isotropy(bx, Geom(lev).CellSize(), rho, grad_rho);
    }
#endif


#ifdef ARK_USE_MOMENTUM
    if(CaseIndex() == 1)
    {
        using Ark::PI;
        const Real *dl = Geom(lev).CellSize();
        MultiFab const &VelFab = ld.msv[Ark::VEL];

        for(amrex::MFIter mfi(VelFab);mfi.isValid();++mfi)
        {
            const Box &bx = mfi.tilebox();
            Array4<Real const> const &vel  = VelFab.const_array(mfi);

            Array4<Real> const &fWave  = ld.fDVDF[Ark::WAVE].array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const double x = dl[Ark::X]*(i+0.5)/ChLength();
                const double y = dl[Ark::Y]*(j+0.5)/ChLength();

                const double U = vel(i,j,k,Ark::X);
                const double V = vel(i,j,k,Ark::Y);

                const double grad_ux = U0()*sin(2.0*PI*x)*sin(2.0*PI*y);
                const double grad_vx = U0()*cos(2.0*PI*x)*cos(2.0*PI*y);
//
                const double grad_uy = -grad_vx;
                const double grad_vy = -grad_ux;

                const double grad_ut = Nu0()*4.0*PI*U0()*cos(2.0*PI*x)*sin(2.0*PI*y);
                const double grad_vt = -Nu0()*4.0*PI*U0()*sin(2.0*PI*x)*cos(2.0*PI*y);

                LoopVS(Ark::DmVnMomentum::Q)
                {
                    const double xiU   = Ark::DmVnMomentum::xi_u[n];
                    const double xiV   = Ark::DmVnMomentum::xi_v[n];
                    const double omega = Ark::DmVnMomentum::omega[n];
                    const double RT    = Ark::DmVnMomentum::RT;

                    const double u1 = xiU*U + xiV*V;
                    const double A_u = xiU + u1*xiU/RT - U;
                    const double A_v = xiV + u1*xiV/RT - V;
                    const double A = omega*1.0/RT; //Rho0 = 1.0
                    const double fEq_t = A * (A_u*grad_ut + A_v*grad_vt);
                    const double fEq_x = A * (A_u*grad_ux + A_v*grad_vx)*xiU;
                    const double fEq_y = A * (A_u*grad_uy + A_v*grad_vy)*xiV;
                    fWave(i,j,k,n) -= TauMomentum()*(fEq_t + fEq_x + fEq_y);
                }
            });
        }
    }
#endif
}
