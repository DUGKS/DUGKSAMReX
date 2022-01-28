#include <DUGKSAmr.H>

void DUGKSAmr::Evolve()
{
	BL_PROFILE_VAR("DUGKSAmr::Evolve()", Evolve);

	Ark::PrintSplitLine();
	amrex::Print() 	<< " case : " << CaseName() <<"\n"
					<< " start on " << amrex::ParallelDescriptor::NProcs()
					<< " processes\n";
	Ark::PrintSplitLine('\n');

	// WriteTecplotGrid(0);
	// WriteTecplotFlowField(0);

	Real begin = amrex::ParallelDescriptor::second();
	Real elapsed_time = 0;

	while(nstep < end_step or steady_state)
	{
        // if (m_regrid_int > 0 and nstep%m_regrid_int == 0)
        // {
        //     if (m_verbose > 2) amrex::Print() << "Regridding...\n";
        //     regrid(0, cur_time);
        //     if (m_verbose > 2 and amrex::ParallelDescriptor::IOProcessor()) {
        //         printGridSummary(amrex::OutStream(), 0, finest_level);
        //     }
        // }

		(this->*Advance) ();

		++nstep;
		cur_time += Dt();

		if(write_conv_interval > 0 and nstep%write_conv_interval == 0) {

			elapsed_time = amrex::ParallelDescriptor::second() - begin;
			begin = amrex::ParallelDescriptor::second();

			WriteComputingTime(elapsed_time);
			// WriteL2NormFile(MultiLevelL2Norm());
			WriteSumPhaseFile(MultiLevelSum());

			// WriteInterfacePosition(InterfacePositionTaylorBubble(finest_level));
			// WriteInterfacePosition(InterfacePositionRayleighTaylor(finest_level));
			// WriteInterfacePosition(InterfacePositionRisingBubble(finest_level));
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
			// WriteTecplotFlowField(0);

			if (m_verbose > 2 and amrex::ParallelDescriptor::IOProcessor()) {
            	printGridSummary(amrex::OutStream(), 0, finest_level);
        	}
		}

		if(nstep > 0 and write_check_interval > 0 and nstep%write_check_interval == 0) {
			WriteCheckPointFile();
		}

	}
	if(end_step%write_plot_interval != 0) {
		WritePlotFile();
	}
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
		if(d > residual_tolerance) {
			steady_state = true;
			break;
		} 
		steady_state = false;
	}

// amrex::ParallelDescriptor::Bcast(&steady_state,1,amrex::ParallelDescriptor::IOProcessorNumber());	
}
void DUGKSAmr::Update_DVDF_Wave (int const lev, const double hTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Wave()", Update_DVDF_Wave);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Wave (ld.hDVDF[Ark::WAVE], ld.hDVDF[Ark::BP],
											ld.get_hDVDFnF_flux_const(),
									  	    Geom(lev).CellSize(), hTime, Dt());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Wave (ld.fDVDF[Ark::WAVE], ld.fDVDF[Ark::BP],
											   ld.get_fDVDFnF_flux_const(),
									  	       Geom(lev).CellSize(), hTime, Dt());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Wave (ld.gDVDF[Ark::WAVE], ld.gDVDF[Ark::BP], 
											 ld.get_gDVDFnF_flux_const(),
									   		 Geom(lev).CellSize(), hTime, Dt());
#endif
}

void DUGKSAmr::Update_DVM_Wave (int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DVMCore<Ark::DmVnPhase>::Update_Wave(ld.hDVDF, Dt(), MinL());
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

	FillPhysicalBoundaryDVDF(lev);
}

// void DUGKSAmr::Update_DVDF_Wave_NonEq (int const lev)
// {
// 	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Wave_NonEq()", Update_DVDF_Wave_NonEq);

// 	auto &ld = Data(lev);

// #ifdef ARK_USE_PHASE
// 	DUGKSCore<Ark::DmVnPhase>::Update_Wave (ld.hDVDF, ld.get_hDVDFnF_flux_const(),
// 									  	    Geom(lev).CellSize(), Dt());
// #endif

// #ifdef ARK_USE_MOMENTUM
// 	DUGKSCore<Ark::DmVnMomentum>::Update_Wave (ld.fDVDF, ld.get_fDVDFnF_flux_const(),
// 									  	       Geom(lev).CellSize(), Dt());
// #endif

// #ifdef ARK_USE_ENERGY
// 	DUGKSCore<Ark::DmVnEnergy>::Update_Wave (ld.gDVDF, ld.get_gDVDFnF_flux_const(),
// 									   		 Geom(lev).CellSize(), Dt());
// #endif
// }

void DUGKSAmr::Update_DVDF_Prime (int const lev, double const upTime, double const downTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Prime()", Update_DVDF_Prime);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	// DUGKSCore<Ark::DmVnPhase>::Update_Prime(ld.hDVDF, upTime, downTime);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Prime(ld.fDVDF, upTime, downTime);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Prime(ld.gDVDF, upTime, downTime);
#endif
}

//hdTime = half of dTime
void DUGKSAmr::Update_DVDF_NonEq (int const lev, double const dTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_NonEq()", Update_DVDF_NonEq);

	const double hdTime = dTime/2.0;

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_NonEq(ld.hDVDF, hdTime);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_NonEq(ld.fDVDF, hdTime);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_NonEq(ld.gDVDF, hdTime);
#endif
}

void DUGKSAmr::Construct_DVDF_Wave (int const lev, const Real dTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Wave()", Update_DVDF_Wave);

	const double hdTime = dTime/2.0;
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Construct_Wave(ld.hDVDF, hdTime);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Construct_Wave(ld.fDVDF, hdTime);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Construct_Wave(ld.gDVDF, hdTime);
#endif
}

void DUGKSAmr::Update_DVDF_BarPlus (int const lev, const double hTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_BarPlus()", Update_DVDF_BarPlus);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_BarPlus(ld.hDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	// Update_Tau_Momentum (ld.fDVDF[Ark::TAU], ld.msv);
	DUGKSCore<Ark::DmVnMomentum>::Update_BarPlus(ld.fDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.fDVDF[Ark::BP], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_BarPlus(ld.gDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.gDVDF[Ark::BP], &LevelData::gDVDF);
#endif
}

void DUGKSAmr::Update_DVDF_BarPlus_NonEq (int const lev, const double hTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_BarPlus_NonEq()", Update_DVDF_BarPlus_NonEq);

	auto &ld = Data(lev);
	const double hhTime = hTime/2.0;

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_BarPlus(ld.hDVDF, hhTime);
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	// Update_Tau_Momentum (ld.fDVDF[Ark::TAU], ld.msv);
	DUGKSCore<Ark::DmVnMomentum>::Update_BarPlus(ld.fDVDF, hhTime);
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.fDVDF[Ark::BP], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_BarPlus(ld.gDVDF, hhTime);
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.gDVDF[Ark::BP], &LevelData::gDVDF);
#endif
}

void DUGKSAmr::Update_DVM_BarPlus(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_BarPlus(ld.hDVDF, Dt(), Dt());
	FillPatchDVDFBarPlus(lev, cur_time, ld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_BarPlus(ld.fDVDF, hDt(lev), Dt());
	FillPatchDVDFBarPlus(lev, cur_time, ld.fDVDF[Ark::BP], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_BarPlus(ld.gDVDF, hDt(lev), Dt());
	FillPatchDVDFBarPlus(lev, cur_time, ld.gDVDF[Ark::BP], &LevelData::gDVDF);
#endif
}

void DUGKSAmr::Update_FD_DVDF_BarPlus (int const lev, double const hTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_FD_DVDF_BarPlus()", Update_FD_DVDF_BarPlus);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_BarPlus(ld.hDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFWave(lev, cur_time, ld.hDVDF[Ark::WAVE], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	// Update_Tau_Momentum(ld.fDVDF[Ark::TAU], ld.msv[Ark::PHI], ld.msv[Ark::RHO]);
	DUGKSCore<Ark::DmVnMomentum>::Update_BarPlus(ld.fDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFWave(lev, cur_time, ld.fDVDF[Ark::WAVE], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_BarPlus(ld.gDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFWave(lev, cur_time, ld.gDVDF[Ark::WAVE], &LevelData::gDVDF);
#endif
}


void DUGKSAmr::Update_DVDF_Flux (int const lev, double const hTime)
{

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Flux()", Update_DVDF_Flux);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_xFluxX(ld.hDVDFxF, hTime);
	DUGKSCore<Ark::DmVnPhase>::Update_xFluxX(ld.hDVDFyF, hTime);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Update_xFluxX(ld.hDVDFzF, hTime);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_xFluxX(ld.fDVDFxF, hTime);
	DUGKSCore<Ark::DmVnMomentum>::Update_xFluxX(ld.fDVDFyF, hTime);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Update_xFluxX(ld.fDVDFzF, hTime);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_xFluxX(ld.gDVDFxF, hTime);
	DUGKSCore<Ark::DmVnEnergy>::Update_xFluxX(ld.gDVDFyF, hTime);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Update_xFluxX(ld.gDVDFzF, hTime);
#endif
#endif
BL_PROFILE_VAR_STOP(Update_DVDF_Flux);

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Flux_Boundary()", Update_DVDF_Flux_Boundary);

	FillPhysicalBoundaryDVDF(lev);

BL_PROFILE_VAR_STOP(Update_DVDF_Flux_Boundary);

}

void DUGKSAmr::Update_DVDF_MesoFlux (int const lev, double const hTime)
{
BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_MesoFlux()", Update_DVDF_MesoFlux);

	(this->*Update_DVDF_Barh) (lev);
	// Update_DVDF_Barh_UW_Limiter(lev);

	Update_MacroVars_hDt (lev, 0.5*hTime);
	Update_DVDFnF_Tau (lev);
	Update_DVDF_Eq_hDt (lev);
	Update_DVDF_Flux (lev, hTime);

BL_PROFILE_VAR_STOP(Update_DVDF_MesoFlux);
}

void DUGKSAmr::Update_DVDF_Wave_SS (int const lev, double const ssTime)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	// Aborted
	// DUGKSCore<Ark::DmVnPhase>::Update_Wave_SS(ld.hDVDF[Ark::WAVE],ld.hDVDF[Ark::SO], ssTime);
#endif
	
	Update_MacroSourceVars(ld.msv, ssTime);

#ifdef ARK_USE_MOMENTUM	
	DUGKSCore<Ark::DmVnMomentum>::Update_Eq_Momentum(ld.fDVDF[Ark::EQ], ld.fDVDF[Ark::SO], ld.msv);

	DUGKSCore<Ark::DmVnMomentum>::Update_Wave_SS(ld.fDVDF[Ark::WAVE],ld.fDVDF[Ark::SO], ssTime);
#endif

#ifdef ARK_USE_ENERGY	
	DUGKSCore<Ark::DmVnEnergy>::Update_Eq_Energy(ld.gDVDF[Ark::EQ], ld.gDVDF[Ark::SO], ld.msv);

	DUGKSCore<Ark::DmVnEnergy>::Update_Wave_SS(ld.gDVDF[Ark::WAVE],ld.gDVDF[Ark::SO], ssTime);
#endif
}

void DUGKSAmr::Update_DVDF_Tau (int const lev)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Tau()", Update_DVDF_Tau);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	Update_Tau_Phase(ld.hDVDF[Ark::TAU], ld.msv);
#endif

#ifdef ARK_USE_MOMENTUM
	Update_Tau_Momentum(ld.fDVDF[Ark::TAU], ld.msv);
#endif

#ifdef ARK_USE_Energy
	Update_Tau_Energy(ld.gDVDF[Ark::TAU], ld.msv);
#endif
}

void DUGKSAmr::Update_DVDFnF_Tau (int const lev)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDFnF_Tau()", Update_DVDFnF_Tau);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	Update_Tau_Phase(ld.hDVDFxF[Ark::TAUF], ld.msv_xF);
	Update_Tau_Phase(ld.hDVDFyF[Ark::TAUF], ld.msv_yF);
#if(AMREX_SPACEDIM == 3)
	Update_Tau_Phase(ld.hDVDFzF[Ark::TAUF], ld.msv_zF);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	Update_Tau_Momentum(ld.fDVDFxF[Ark::TAUF], ld.msv_xF);
	Update_Tau_Momentum(ld.fDVDFyF[Ark::TAUF], ld.msv_yF);
#if(AMREX_SPACEDIM == 3)
	Update_Tau_Momentum(ld.fDVDFzF[Ark::TAUF], ld.msv_zF);
#endif
#endif

#ifdef ARK_USE_Energy
	Update_Tau_Energy(ld.gDVDFxF[Ark::TAUF], ld.msv_xF);
	Update_Tau_Energy(ld.gDVDFyF[Ark::TAUF], ld.msv_yF);
#if(AMREX_SPACEDIM == 3)
	Update_Tau_Energy(ld.gDVDFzF[Ark::TAUF], ld.msv_zF);
#endif
#endif
}

void DUGKSAmr::Advance_CN()
{
BL_PROFILE_VAR("DUGKSAmr::Advance_CN_FillMacroVars()", Advance_CN_FillMacroVars);

	AverageDownMacroVars();

	FillPatchMacroVars();

BL_PROFILE_VAR_STOP(Advance_CN_FillMacroVars);

BL_PROFILE_VAR("DUGKSAmr::Advance_CN_DVDF()", Advance_CN_DVDF);

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Eq (lev);
	}

	for (int lev = finest_level; lev > 0; --lev) {

		int const crse_lev = lev-1;
		
		Update_DVDF_BarPlus (crse_lev, hDt(lev));
		Update_DVDF_BarPlus (lev, hDt(lev));
	}
		Update_DVDF_BarPlus (0, hDt(0));


	for (int lev = 0; lev <= finest_level; ++lev) {

		(this->*Update_DVDF_Grad_BarPlus) (lev);
		Update_DVDF_MesoFlux (lev, hDt(lev));
	}

	// CheckForPhysicalBoundary(0);

	AverageDownDVDFFlux();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Wave (lev, hDt(lev));
	}

	AverageDownDVDFWave();

	if (m_regrid_int > 0 and nstep%m_regrid_int == 0)
    {
        if (m_verbose > 100) {
        	amrex::Print() << "Regridding...\n";
      	}
      	
        regrid(0, cur_time);
    }

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);
	}

BL_PROFILE_VAR_STOP(Advance_CN_DVDF);

}
// void DUGKSAmr::Advance_ForwardEuler()
// {
// 	AverageDownMacroVars();

// 	FillPatchMacroVars();

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_DVDF_Eq (lev);

// 		Update_DVDF_BarPlus (lev, 2*hDt(lev));

// 		(this->*Update_DVDF_Grad_BarPlus) (lev);

// 		Update_DVDF_MesoFlux (lev, 0);
// 	}

// 	// CheckForPhysicalBoundary(0);

// 	AverageDownDVDFFlux();

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_DVDF_Wave (lev, 2*hDt(lev));
// 	}

// 	AverageDownDVDFWave();

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_MacroVars (lev);
// 	}
// }

void DUGKSAmr::Advance_BackwardEuler()
{
	AverageDownMacroVars();

	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Eq (lev);

		Update_DVDF_BarPlus (lev, 0);

		(this->*Update_DVDF_Grad_BarPlus) (lev);

		Update_DVDF_MesoFlux (lev, 2*hDt(lev));
	}

	// CheckForPhysicalBoundary(0);

	AverageDownDVDFFlux();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Wave (lev, 0);
	}

	AverageDownDVDFWave();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);
	}
}

void DUGKSAmr::Advance_IMEX() {}

void DUGKSAmr::Advance_SS()
{
	AverageDownMacroVars();

	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {

		Update_DVDF_Wave_SS (lev, 0.5*Dt());

		Nullify_Source (lev);
	}

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Eq (lev);

		Update_DVDF_Prime (lev, -Dt(), 0);
	}

	for (int lev = finest_level; lev > 0; --lev) {

		int const crse_lev = lev-1;
		
		Update_DVDF_BarPlus (crse_lev, hDt(lev));
		Update_DVDF_BarPlus (lev, hDt(lev));
	}

		Update_DVDF_BarPlus (0, hDt(0));

	for (int lev = 0; lev <= finest_level; ++lev) {

		(this->*Update_DVDF_Grad_BarPlus) (lev);

		Update_DVDF_MesoFlux (lev, hDt(lev));
	}

	// CheckForPhysicalBoundary(0);

	AverageDownDVDFFlux();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Wave (lev, hDt(lev));
	}

	AverageDownDVDFWave();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);
	}

	//primitive DVDF, stored in xWave
	for (int lev = 0; lev <= finest_level; ++lev) {

		Update_DVDF_Eq (lev);

		Update_DVDF_Prime (lev, Dt(), Dt());
	}

	for (int lev = 0; lev <= finest_level; ++lev) {

		Update_DVDF_Wave_SS (lev, 0.5*Dt());
	}

}

void DUGKSAmr::Advance_FD()
{
	AverageDownMacroVars();

	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {

		Update_DVDF_Eq (lev);

		Update_FD_DVDF_BarPlus (lev, 0);
		
		Update_DVDF_MesoFlux (lev, Dt());
	}

	// CheckForPhysicalBoundary(0);

	AverageDownDVDFFlux();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Wave (lev, 0);
	}

	AverageDownDVDFWave();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);
	}
}

void DUGKSAmr::Advance_DVM()
{
	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {

		Update_DVDF_Eq (lev);

		Update_DVM_BarPlus (lev);

		(this->*Update_DVDF_Grad_BarPlus) (lev);
		
#ifdef ARK_USE_MOMENTUM
		Update_DVDF_MesoFlux (lev, hDt(lev));
#endif
		Update_DVM_Wave (lev);

		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);

	}
}
void DUGKSAmr::Advance_NEQ()
{
BL_PROFILE_VAR("DUGKSAmr::Advance_NEQ_Regird()", Advance_CN_Regrid);
	
	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Eq (lev);
		Update_DVDF_NonEq (lev, Dt());
	}

	if (m_regrid_int > 0 and nstep%m_regrid_int == 0)
    {
        if (m_verbose > 100) {
        	amrex::Print() << "Regridding...\n";
      	}

        regrid(0, cur_time);
    }
    
	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Tau (lev);
		Update_DVDF_Eq (lev);
		Construct_DVDF_Wave(lev, Dt());
	}

BL_PROFILE_VAR_STOP(Advance_NEQ_Regird);

BL_PROFILE_VAR("DUGKSAmr::Advance_NEQ_DVDF()", Advance_CN_DVDF);

	for (int lev = finest_level; lev > 0; --lev) {

		int const crse_lev = lev-1;
		
		Update_DVDF_BarPlus (crse_lev, hDt(lev));
		Update_DVDF_BarPlus (lev, hDt(lev));
	}
		Update_DVDF_BarPlus (0, hDt(0));


	for (int lev = 0; lev <= finest_level; ++lev) {

		(this->*Update_DVDF_Grad_BarPlus) (lev);
		Update_DVDF_MesoFlux (lev, hDt(lev));
	}

	// CheckForPhysicalBoundary(0);

	AverageDownDVDFFlux();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Wave (lev,hDt(lev));
	}

	AverageDownDVDFWave();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);
	}

	AverageDownMacroVars();

BL_PROFILE_VAR_STOP(Advance_NEQ_DVDF);

}
