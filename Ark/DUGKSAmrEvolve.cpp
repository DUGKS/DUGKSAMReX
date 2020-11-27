#include <DUGKSAmr.H>

void DUGKSAmr::Update_DVDF_Wave (int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Wave (ld.hDVDF[Ark::WAVE], ld.hDVDF[Ark::BP],
											ld.get_hDVDFnF_flux_const(),
									  	    Geom(lev).CellSize(), hDt(lev), Dt());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Wave (ld.fDVDF[Ark::WAVE], ld.fDVDF[Ark::BP],
											   ld.get_fDVDFnF_flux_const(),
									  	       Geom(lev).CellSize(), hDt(lev), Dt());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Wave (ld.gDVDF[Ark::WAVE], ld.gDVDF[Ark::BP], 
											 ld.get_gDVDFnF_flux_const(),
									   		 Geom(lev).CellSize(), hDt(lev), Dt());
#endif
}

void DUGKSAmr::Update_DVDF_BarPlus(int const lev, double const hTime)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_BarPlus(ld.hDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	// Update_Tau_Momentum(ld.fDVDF[Ark::TAU], ld.msv[Ark::PHI], ld.msv[Ark::RHO]);
	DUGKSCore<Ark::DmVnMomentum>::Update_BarPlus(ld.fDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.fDVDF[Ark::BP], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_BarPlus(ld.gDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.gDVDF[Ark::BP], &LevelData::gDVDF);
#endif
}

void DUGKSAmr::Update_DVDF_Grad_BarPlus (int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Grad_BarPlus(ld.hDVDF, Geom(lev).CellSize());
	FillInteriorBoundaryGradBarPlus(ld.hDVDF, Geom(lev).periodicity());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Grad_BarPlus(ld.fDVDF, Geom(lev).CellSize());
	FillInteriorBoundaryGradBarPlus(ld.fDVDF, Geom(lev).periodicity());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Grad_BarPlus(ld.gDVDF, Geom(lev).CellSize());
	FillInteriorBoundaryGradBarPlus(ld.gDVDF, Geom(lev).periodicity());
#endif	
}

void DUGKSAmr::Update_DVDF_Barh_UW (int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW (ld.hDVDFxF[Ark::FLUX], ld.hDVDF,
						 				  Geom(lev).CellSize(), hDt(lev),
						 				  amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW (ld.hDVDFyF[Ark::FLUX], ld.hDVDF,
						 				  Geom(lev).CellSize(), hDt(lev),
						 				  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_UW (ld.hDVDFzF[Ark::FLUX], ld.hDVDF,
						 				  Geom(lev).CellSize(), hDt(lev),
						 				  amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW (ld.fDVDFxF[Ark::FLUX], ld.fDVDF,
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW (ld.fDVDFyF[Ark::FLUX], ld.fDVDF,
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_UW (ld.fDVDFzF[Ark::FLUX], ld.fDVDF,
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW (ld.gDVDFxF[Ark::FLUX], ld.gDVDF,
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW (ld.gDVDFyF[Ark::FLUX], ld.gDVDF,
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_UW (ld.gDVDFzF[Ark::FLUX], ld.gDVDF,
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::z);
#endif

#endif
}

void DUGKSAmr::Update_DVDF_Barh_CD (int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_CD (ld.hDVDFxF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  Geom(lev).CellSize(), hDt(lev),
						 				  amrex::Direction::x);

	DUGKSCore<Ark::DmVnPhase>::Update_Barh_CD (ld.hDVDFyF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  Geom(lev).CellSize(), hDt(lev),
						 				  amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_Barh_CD (ld.hDVDFzF[Ark::FLUX], ld.hDVDF[Ark::BP],
						 				  Geom(lev).CellSize(), hDt(lev),
						 				  amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_CD (ld.fDVDFxF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::x);

	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_CD (ld.fDVDFyF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_Barh_CD (ld.fDVDFzF[Ark::FLUX], ld.fDVDF[Ark::BP],
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::z);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_CD (ld.gDVDFxF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::x);

	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_CD (ld.gDVDFyF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::y);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_Barh_CD (ld.gDVDFzF[Ark::FLUX], ld.gDVDF[Ark::BP],
						 				  	 Geom(lev).CellSize(), hDt(lev),
						 				  	 amrex::Direction::z);
#endif

#endif
}
void DUGKSAmr::Update_DVDF_Flux (int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_xFluxX(ld.hDVDFxF, hDt(lev));
	DUGKSCore<Ark::DmVnPhase>::Update_xFluxX(ld.hDVDFyF, hDt(lev));
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_xFluxX(ld.hDVDFzF, hDt(lev));
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_xFluxX(ld.fDVDFxF, hDt(lev));
	DUGKSCore<Ark::DmVnMomentum>::Update_xFluxX(ld.fDVDFyF, hDt(lev));
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_xFluxX(ld.fDVDFzF, hDt(lev));
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_xFluxX(ld.gDVDFxF, hDt(lev));
	DUGKSCore<Ark::DmVnEnergy>::Update_xFluxX(ld.gDVDFyF, hDt(lev));
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_xFluxX(ld.gDVDFzF, hDt(lev));
#endif
#endif

	FillPhysicalBoundaryDVDF(lev);
}

void DUGKSAmr::Update_DVDF_MesoFlux (int const lev)
{
	(this->*Update_DVDF_Barh) (lev);
	Update_MacroVars_hDt (lev);
	Update_DVDF_Eq_hDt (lev);
	Update_DVDF_Flux (lev);

}

void DUGKSAmr::Advance()
{
	AverageDownMacroVars();

	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Eq (lev);
	}

	for (int lev = finest_level; lev >= 0; --lev) {

		int const crse_lev = (lev == 0 ? lev : lev-1);
		
		Update_DVDF_BarPlus (crse_lev, hDt(lev));
		Update_DVDF_BarPlus (lev, hDt(lev));
	}

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Grad_BarPlus (lev);
		Update_DVDF_MesoFlux (lev);
	}

	// CheckForPhysicalBoundary(0);

	AverageDownDVDFFlux();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Wave (lev);
	}

	AverageDownDVDFWave();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_MacroVars (lev);
	}
}
void DUGKSAmr::Evolve()
{
	BL_PROFILE("DUGKSAmr::Evolve()");

	Ark::PrintSplitLine();
	amrex::Print() 	<< " case : " << CaseName() <<"\n"
					<< " start on " << amrex::ParallelDescriptor::NProcs()
					<< " processes\n";
	Ark::PrintSplitLine('\n');
//
	while(nstep < end_step or steady_state)
	{
        if (m_regrid_int > 0 and nstep%m_regrid_int == 0)
        {
            if (m_verbose > 2) amrex::Print() << "Regridding...\n";
            regrid(0, cur_time);
            if (m_verbose > 2 and amrex::ParallelDescriptor::IOProcessor()) {
                printGridSummary(amrex::OutStream(), 0, finest_level);
            }
        }

		Advance();

		++nstep;
		cur_time += Dt();

		if(write_conv_interval > 0 and nstep%write_conv_interval == 0) {
			WriteL2NormFile(MultiLevelL2Norm());
			WriteSumPhaseFile(MultiLevelSum());
		}
		if(write_residual_interval > 0 and nstep%write_residual_interval == 0) {
			const auto &level_residuals = MultiLevelResidual();

			if(steady_state) {
				WriteResidualFile(level_residuals);
				SteadyStateReached(level_residuals);
			}
		}
		
		if(write_plot_interval > 0 and nstep%write_plot_interval == 0) {
			WritePlotFile();
		}

		if(write_check_interval > 0 and nstep%write_check_interval == 0) {
			WriteCheckPointFile();
		}

	}

	WritePlotFile();
//
	Ark::PrintSplitLine();
	amrex::Print() 	<< " case : " << CaseName() << " complete\n";
}
void DUGKSAmr::SteadyStateReached(const Vector<double> & residualv)
{
	if (m_verbose > 0) {
        amrex::Print().SetPrecision(8) 	
        			<< " step = "<<nstep <<"    "
        			<< " time = "<<cur_time<<"    "
        			<< " residual : ";            
        
        for(auto d : residualv) {	
        	amrex::Print() << d <<"    ";
        }
        amrex::Print() << "\n";
    }
	for(auto d : residualv)
	{
		
		if(d > residual_tolerance) 
			break;
		steady_state = false;
	}
	
}