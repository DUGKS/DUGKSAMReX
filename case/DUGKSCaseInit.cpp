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
            break;
        }
        case 10:
        {
            QuiescentDrop(lev);
            break;
        }
        case 11:
        {
            TransDrop(lev);
            break;
        }
        case 12:
        {
            ZalesakDisk(lev);
            break;
        }
        case 13:
        {
            InterfaceElongation(lev);
            break;
        }
        case 14:
        {
            SmoothDeformation(lev);
            break;
        }
        case 15:
        {
            LayeredPoiseuille(lev);
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
#endif

#ifdef ARK_USE_MOMENTUM
    ld.fDVDF[Ark::TAU ].setVal(TauMomentum());
	ld.fDVDFxF[Ark::TAUF].setVal(TauMomentum());
	ld.fDVDFyF[Ark::TAUF].setVal(TauMomentum());
#if(AMREX_SPACEDIM == 3)
    ld.fDVDFzF[Ark::TAUF].setVal(TauMomentum());
#endif
#endif

#ifdef ARK_USE_ENERGY
    ld.gDVDF[Ark::TAU ].setVal(TauEnergy());
	ld.gDVDFxF[Ark::TAUF].setVal(TauEnergy());
	ld.gDVDFyF[Ark::TAUF].setVal(TauEnergy());
#if(AMREX_SPACEDIM == 3)
    ld.gDVDFzF[Ark::TAUF].setVal(TauEnergy());
#endif
#endif	
}

void DUGKSAmr::case_init_mesoscopic_equilibrium(int const lev)
{
	auto& ld = *m_leveldata[lev];

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
    auto& ld = *m_leveldata[lev];

    ld.phase_prev.setVal(0);
    ld.velocity_prev.setVal(0);
    
    //The following would result in a slight increase in L-2 norm error;
    // MultiFab::Copy(ld.phase_prev, ld.msv[Ark::PHI], 0, 0, 1, unit_ghost());
    // MultiFab::Copy(ld.velocity_prev, ld.msv[Ark::VEL], 0, 0, AMREX_SPACEDIM, unit_ghost());

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
}
