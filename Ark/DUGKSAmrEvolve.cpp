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
			WriteL2NormFile((this->*MultiLevelL2Norm)());
			// WriteSumPhaseFile(MultiLevelSum());
			WriteSumDensityFile(MultiLevelSumDensity());

			// WriteInterfacePosition(InterfacePositionTaylorBubble(finest_level));
			// WriteInterfacePosition(InterfacePositionRayleighTaylor(finest_level));
			// WriteInterfacePosition(InterfacePositionLiquidBridge(finest_level));
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

void DUGKSAmr::Update_DVDF_Wave_Lite (int const lev)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Wave_Lite()", Update_DVDF_Wave_Lite);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Update_Wave (ld.hDVDF[Ark::WAVE],
											ld.get_hDVDFnF_flux_const(),
									  	    Geom(lev).CellSize(), Dt());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Wave (ld.fDVDF[Ark::WAVE],
											   ld.get_fDVDFnF_flux_const(),
									  	       Geom(lev).CellSize(), Dt());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Update_Wave (ld.gDVDF[Ark::WAVE],
											 ld.get_gDVDFnF_flux_const(),
									   		 Geom(lev).CellSize(), Dt());
#endif
}

void DUGKSAmr::Update_DVM_Wave (int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DVMCore<Ark::DmVnPhase>::Update_Wave(ld.hDVDF, Dt(), MinL());
#endif

#ifdef ARK_USE_MOMENTUM
	// DUGKSCore<Ark::DmVnMomentum>::Update_Wave (ld.fDVDF[Ark::WAVE], ld.fDVDF[Ark::BP],
	// 										   ld.get_fDVDFnF_flux_const(),
	// 								  	       Geom(lev).CellSize(), hDt(lev), Dt());

	DVMCore<Ark::DmVnMomentum>::Update_Wave(ld.fDVDF, Dt(), MinL());
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
void DUGKSAmr::Transform_DVDF_Moment(int const lev, Ark::EnumDVDF ST)
{
	BL_PROFILE_VAR("DUGKSAmr::Transform_DVDF_Moment()", Transform_DVDF_Moment);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Transform (ld.hDVDF[ST], Ark::DmVnPhase::M);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Transform (ld.fDVDF[ST], Ark::DmVnMomentum::M);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Transform (ld.gDVDF[ST], Ark::DmVnEnergy::M);
#endif
}

void DUGKSAmr::Transform_Moment_DVDF(int const lev, Ark::EnumDVDF ST)
{
	BL_PROFILE_VAR("DUGKSAmr::Transform_Moment_DVDF()", Transform_Moment_DVDF);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Transform (ld.hDVDF[ST], Ark::DmVnPhase::IM);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Transform (ld.fDVDF[ST], Ark::DmVnMomentum::IM);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Transform (ld.gDVDF[ST], Ark::DmVnEnergy::IM);
#endif
}
void DUGKSAmr::Transform_DVDF_Moment_hDt(int const lev, Ark::EnumFaceDVDF ST)
{
	BL_PROFILE_VAR("DUGKSAmr::Transform_DVDF_Moment_hDt()", Transform_DVDF_Moment_hDt);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Transform (ld.hDVDFxF[ST], Ark::DmVnPhase::M);
	DUGKSCore<Ark::DmVnPhase>::Transform (ld.hDVDFyF[ST], Ark::DmVnPhase::M);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Transform (ld.hDVDFzF[ST], Ark::DmVnPhase::M);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Transform (ld.fDVDFxF[ST], Ark::DmVnMomentum::M);
	DUGKSCore<Ark::DmVnMomentum>::Transform (ld.fDVDFyF[ST], Ark::DmVnMomentum::M);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Transform (ld.fDVDFzF[ST], Ark::DmVnMomentum::M);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Transform (ld.gDVDFxF[ST], Ark::DmVnEnergy::M);
	DUGKSCore<Ark::DmVnEnergy>::Transform (ld.gDVDFyF[ST], Ark::DmVnEnergy::M);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Transform (ld.gDVDFzF[ST], Ark::DmVnEnergy::M);
#endif
#endif
}

void DUGKSAmr::Transform_Moment_DVDF_hDt(int const lev, Ark::EnumFaceDVDF ST)
{
	BL_PROFILE_VAR("DUGKSAmr::Transform_Moment_DVDF_hDt()", Transform_Moment_DVDF_hDt);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Transform (ld.hDVDFxF[ST], Ark::DmVnPhase::IM);
	DUGKSCore<Ark::DmVnPhase>::Transform (ld.hDVDFyF[ST], Ark::DmVnPhase::IM);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Transform (ld.hDVDFzF[ST], Ark::DmVnPhase::IM);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Transform (ld.fDVDFxF[ST], Ark::DmVnMomentum::IM);
	DUGKSCore<Ark::DmVnMomentum>::Transform (ld.fDVDFyF[ST], Ark::DmVnMomentum::IM);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Transform (ld.fDVDFzF[ST], Ark::DmVnMomentum::IM);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Transform (ld.gDVDFxF[ST], Ark::DmVnEnergy::IM);
	DUGKSCore<Ark::DmVnEnergy>::Transform (ld.gDVDFyF[ST], Ark::DmVnEnergy::IM);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Transform (ld.gDVDFzF[ST], Ark::DmVnEnergy::IM);
#endif
#endif
}


template<bool ENABLESO>
void DUGKSAmr::Construct_DVDF_Prime (int const lev, double const dTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Construct_DVDF_Prime()", Construct_DVDF_Prime);

	const double hdTime = dTime/2.0;

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	// DUGKSCore<Ark::DmVnPhase>::Construct_Prime<ENABLESO> (ld.hDVDF, upTime);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Construct_Prime<ENABLESO> (ld.fDVDF, hdTime);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Construct_Prime<ENABLESO> (ld.gDVDF, hdTime);
#endif
}

//hdTime = half of dTime
template<bool ENABLESO>
void DUGKSAmr::Construct_DVDF_NonEq (int const lev, double const dTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Construct_DVDF_NonEq()", Construct_DVDF_NonEq);

	const double hdTime = dTime/2.0;

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Construct_NonEq<ENABLESO> (ld.hDVDF, hdTime);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Construct_NonEq<ENABLESO> (ld.fDVDF, hdTime);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Construct_NonEq<ENABLESO> (ld.gDVDF, hdTime);
#endif
}

void DUGKSAmr::Construct_DVDF_WavePlus (int const lev, double const hTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Construct_DVDF_WavePlus()", Construct_DVDF_WavePlus);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Construct_WavePlus(ld.hDVDF[Ark::WAVE], ld.hDVDF[Ark::BP],
												  hTime, Dt());
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Construct_WavePlus(ld.fDVDF[Ark::WAVE], ld.fDVDF[Ark::BP],
													 hTime, Dt());
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Construct_WavePlus(ld.gDVDF[Ark::WAVE], ld.gDVDF[Ark::BP],
												   hTime, Dt());
#endif
}

template<bool ENABLESO>
void DUGKSAmr::Construct_DVDF_Wave (int const lev, const Real dTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Construct_DVDF_Wave()", Construct_DVDF_Wave);

	const double hdTime = dTime/2.0;
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Construct_Wave<ENABLESO> (ld.hDVDF, hdTime);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Construct_Wave<ENABLESO> (ld.fDVDF, hdTime);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Construct_Wave<ENABLESO> (ld.gDVDF, hdTime);
#endif
}

void DUGKSAmr::Collision_DVDF_BarPlus (int const lev, const double hTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Collision_DVDF_BarPlus()", Collision_DVDF_BarPlus);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Collision_BarPlus(ld.hDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	// Update_Tau_Momentum (ld.fDVDF[Ark::TAU], ld.msv);
	DUGKSCore<Ark::DmVnMomentum>::Collision_BarPlus(ld.fDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.fDVDF[Ark::BP], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Collision_BarPlus(ld.gDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.gDVDF[Ark::BP], &LevelData::gDVDF);
#endif
}

void DUGKSAmr::Collision_DVDF_BarPlus_Moment (int const lev, const double hTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Collision_DVDF_BarPlus_Moment()", Collision_DVDF_BarPlus_Moment);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Collision_BarPlus(ld.hDVDF, hTime, Dt());
	DUGKSCore<Ark::DmVnPhase>::Transform(ld.hDVDF[Ark::BP], Ark::DmVnPhase::IM);
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	// Update_Tau_Momentum (ld.fDVDF[Ark::TAU], ld.msv);
	DUGKSCore<Ark::DmVnMomentum>::Collision_BarPlus(ld.fDVDF, hTime, Dt());
	DUGKSCore<Ark::DmVnMomentum>::Transform(ld.fDVDF[Ark::BP], Ark::DmVnMomentum::IM);
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.fDVDF[Ark::BP], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Collision_BarPlus(ld.gDVDF, hTime, Dt());
	DUGKSCore<Ark::DmVnEnergy>::Transform(ld.gDVDF[Ark::BP], Ark::DmVnEnergy::IM);
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.gDVDF[Ark::BP], &LevelData::gDVDF);
#endif
}

void DUGKSAmr::Collision_DVDF_BarPlus_NonEq (int const lev, const double hTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Collision_DVDF_BarPlus_NonEq()", Collision_DVDF_BarPlus_NonEq);

	auto &ld = Data(lev);
	const double hhTime = hTime/2.0;

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Collision_BarPlus(ld.hDVDF, hhTime);
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	// Update_Tau_Momentum (ld.fDVDF[Ark::TAU], ld.msv);
	DUGKSCore<Ark::DmVnMomentum>::Collision_BarPlus(ld.fDVDF, hhTime);
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.fDVDF[Ark::BP], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Collision_BarPlus(ld.gDVDF, hhTime);
	DUGKSAmr::FillPatchDVDFBarPlus(lev, cur_time, ld.gDVDF[Ark::BP], &LevelData::gDVDF);
#endif
}

void DUGKSAmr::Collision_DVM_BarPlus(int const lev)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Collision_BarPlus(ld.hDVDF, Dt(), Dt());
	FillPatchDVDFBarPlus(lev, cur_time, ld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Collision_BarPlus(ld.fDVDF, Dt(), Dt());
	FillPatchDVDFBarPlus(lev, cur_time, ld.fDVDF[Ark::BP], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Collision_BarPlus(ld.gDVDF, Dt(), Dt());
	FillPatchDVDFBarPlus(lev, cur_time, ld.gDVDF[Ark::BP], &LevelData::gDVDF);
#endif
}

void DUGKSAmr::Update_FD_DVDF_BarPlus (int const lev, double const hTime)
{
	BL_PROFILE_VAR("DUGKSAmr::Update_FD_DVDF_BarPlus()", Update_FD_DVDF_BarPlus);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Collision_BarPlus(ld.hDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFWave(lev, cur_time, ld.hDVDF[Ark::WAVE], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
	// Update_Tau_Momentum(ld.fDVDF[Ark::TAU], ld.msv[Ark::PHI], ld.msv[Ark::RHO]);
	DUGKSCore<Ark::DmVnMomentum>::Collision_BarPlus(ld.fDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFWave(lev, cur_time, ld.fDVDF[Ark::WAVE], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Collision_BarPlus(ld.gDVDF, hTime, Dt());
	DUGKSAmr::FillPatchDVDFWave(lev, cur_time, ld.gDVDF[Ark::WAVE], &LevelData::gDVDF);
#endif
}

template <bool ENABLESOURCE>
void DUGKSAmr::Collision_DVDF_Flux (int const lev, double const hTime)
{

BL_PROFILE_VAR("DUGKSAmr::Collision_DVDF_Flux()", Collision_DVDF_Flux);

	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	DUGKSCore<Ark::DmVnPhase>::Collision_xFluxX<ENABLESOURCE>(ld.hDVDFxF, hTime);
	DUGKSCore<Ark::DmVnPhase>::Collision_xFluxX<ENABLESOURCE>(ld.hDVDFyF, hTime);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnPhase>::Collision_xFluxX<ENABLESOURCE>(ld.hDVDFzF, hTime);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Collision_xFluxX<ENABLESOURCE>(ld.fDVDFxF, hTime);
	DUGKSCore<Ark::DmVnMomentum>::Collision_xFluxX<ENABLESOURCE>(ld.fDVDFyF, hTime);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnMomentum>::Collision_xFluxX<ENABLESOURCE>(ld.fDVDFzF, hTime);
#endif
#endif

#ifdef ARK_USE_ENERGY
	DUGKSCore<Ark::DmVnEnergy>::Collision_xFluxX<ENABLESOURCE>(ld.gDVDFxF, hTime);
	DUGKSCore<Ark::DmVnEnergy>::Collision_xFluxX<ENABLESOURCE>(ld.gDVDFyF, hTime);
#if(AMREX_SPACEDIM == 3)
	DUGKSCore<Ark::DmVnEnergy>::Collision_xFluxX<ENABLESOURCE>(ld.gDVDFzF, hTime);
#endif
#endif
BL_PROFILE_VAR_STOP(Collision_DVDF_Flux);

BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_Flux_Boundary()", Update_DVDF_Flux_Boundary);

	FillPhysicalBoundaryDVDF(lev);

BL_PROFILE_VAR_STOP(Update_DVDF_Flux_Boundary);

}

template <bool ENABLESOURCE>
void DUGKSAmr::Update_DVDF_MesoFlux (int const lev, double const hTime)
{
BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_MesoFlux()", Update_DVDF_MesoFlux);

	(this->*Update_DVDF_Barh) (lev);
	// Update_DVDF_Barh_UW_Limiter(lev);
	Update_MacroVars_hDt (lev, 0.5*hTime);
	Update_DVDFnF_Tau (lev);
	if constexpr (ENABLESOURCE) { Update_DVDF_EqSo_hDt (lev); }
	else { Update_DVDF_Eq_hDt (lev); }
	Collision_DVDF_Flux<ENABLESOURCE> (lev, hTime);

BL_PROFILE_VAR_STOP(Update_DVDF_MesoFlux);
}

template <bool ENABLESOURCE>
void DUGKSAmr::Update_DVDF_MesoFlux_Moment (int const lev, double const hTime)
{
BL_PROFILE_VAR("DUGKSAmr::Update_DVDF_MesoFlux_Moment()", Update_DVDF_MesoFlux_Moment);

	(this->*Update_DVDF_Barh) (lev);
	// Update_DVDF_Barh_UW_Limiter(lev);
	Transform_DVDF_Moment_hDt(lev, Ark::FLUX);
	Update_MacroVars_hDt(lev, 0.5*hTime);
	Update_DVDFnF_Tau (lev);
	if constexpr (ENABLESOURCE) { Update_DVDF_EqSo_hDt (lev); }
	else { Update_DVDF_Eq_hDt (lev); }
	Collision_DVDF_Flux<ENABLESOURCE> (lev, hTime);
	Transform_Moment_DVDF_hDt(lev, Ark::FLUX);

BL_PROFILE_VAR_STOP(Update_DVDF_MesoFlux_Moment);
}

void DUGKSAmr::Update_DVDF_Wave_SS (int const lev, double const ssTime)
{
	auto &ld = Data(lev);

#ifdef ARK_USE_PHASE
	// Aborted
	// DUGKSCore<Ark::DmVnPhase>::Update_Wave_SS(ld.hDVDF[Ark::WAVE],ld.hDVDF[Ark::SO], ssTime);
#endif

	Update_MacroSourceVars(ld.msv, ssTime);

	Update_DVDF_EqSo (lev);

#ifdef ARK_USE_MOMENTUM
	DUGKSCore<Ark::DmVnMomentum>::Update_Wave_SS(ld.fDVDF[Ark::WAVE],ld.fDVDF[Ark::SO], ssTime);
#endif

#ifdef ARK_USE_ENERGY
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
		Update_DVDF_EqSo (lev);
	}

	for (int lev = finest_level; lev > 0; --lev) {

		int const crse_lev = lev-1;

		Collision_DVDF_BarPlus (crse_lev, hDt(lev));
		Collision_DVDF_BarPlus (lev, hDt(lev));
	}
		Collision_DVDF_BarPlus (0, hDt(0));


	for (int lev = 0; lev <= finest_level; ++lev) {

		(this->*Update_DVDF_Grad_BarPlus) (lev);
		Update_DVDF_MesoFlux<> (lev, hDt(lev));
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
		Update_DVDF_EqSo (lev);

		Collision_DVDF_BarPlus (lev, 0);

		(this->*Update_DVDF_Grad_BarPlus) (lev);

		Update_DVDF_MesoFlux<true> (lev, 2*hDt(lev));
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
	for (int lev = 0; lev <= finest_level; ++lev) {

		Update_DVDF_Wave_SS (lev, 0.5*Dt());

		Nullify_Source (lev);
	}

//  regridding
	AverageDownMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Construct_DVDF_NonEq<false> (lev, 0);
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
		Construct_DVDF_Wave<false> (lev, Dt());
	}

	for (int lev = finest_level; lev > 0; --lev) {

		int const crse_lev = lev-1;

		Collision_DVDF_BarPlus (crse_lev, hDt(lev));
		Collision_DVDF_BarPlus (lev, hDt(lev));
	}

		Collision_DVDF_BarPlus (0, hDt(0));

	for (int lev = 0; lev <= finest_level; ++lev) {

		(this->*Update_DVDF_Grad_BarPlus) (lev);

		Update_DVDF_MesoFlux<false> (lev, hDt(lev));
	}

	// CheckForPhysicalBoundary(0);

	AverageDownDVDFFlux();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Wave (lev, hDt(lev));
	}

	for (int lev = 0; lev <= finest_level; ++lev) {

		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);
		// FreeEnergySource(lev);
	}

	// FillPatchMacroVars();

	//Primitive distribution, stored in xWave
	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Eq (lev);
		Construct_DVDF_Prime<false> (lev, Dt());
	}

	// AverageDownMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Wave_SS (lev, 0.5*Dt());
	}
}

void DUGKSAmr::Advance_NEQ()
{
BL_PROFILE_VAR("DUGKSAmr::Advance_NEQ_Regird()", Advance_NEQ_Regrid);

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_EqSo (lev);
		Construct_DVDF_NonEq(lev, Dt());
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
		Update_DVDF_EqSo (lev);
		Construct_DVDF_Wave(lev, Dt());
	}

BL_PROFILE_VAR_STOP(Advance_NEQ_Regird);

BL_PROFILE_VAR("DUGKSAmr::Advance_NEQ_DVDF()", Advance_CN_DVDF);

	for (int lev = finest_level; lev > 0; --lev) {

		int const crse_lev = lev-1;

		Collision_DVDF_BarPlus (crse_lev, hDt(lev));
		Collision_DVDF_BarPlus (lev, hDt(lev));
	}
		Collision_DVDF_BarPlus (0, hDt(0));

	for (int lev = 0; lev <= finest_level; ++lev) {

		(this->*Update_DVDF_Grad_BarPlus) (lev);

		Update_DVDF_MesoFlux<true> (lev, hDt(lev));
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

void DUGKSAmr::Advance_FD()
{
	AverageDownMacroVars();

	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {

		Update_DVDF_EqSo (lev);

		Update_FD_DVDF_BarPlus (lev, 0);

		Update_DVDF_MesoFlux<true> (lev, Dt());
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

// void DUGKSAmr::Advance_SNEQ_S()
// {
// BL_PROFILE_VAR("DUGKSAmr::Advance_SNEQ_Regird()", Advance_CN_Regrid);

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_DVDF_Eq (lev);
// 		Construct_DVDF_NonEq (lev, Dt());
// 	}

// 	if (m_regrid_int > 0 and nstep%m_regrid_int == 0)
//     {
//         if (m_verbose > 100) {
//         	amrex::Print() << "Regridding...\n";
//       	}

//         regrid(0, cur_time);
//     }

// 	FillPatchMacroVars();

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_DVDF_Tau (lev);
// 		Update_DVDF_Eq (lev);
// 		Construct_DVDF_Wave(lev, Dt(), true);
// 	}

// BL_PROFILE_VAR_STOP(Advance_SNEQ_Regird);

// BL_PROFILE_VAR("DUGKSAmr::Advance_SNEQ_DVDF()", Advance_CN_DVDF);

// 	for (int lev = finest_level; lev > 0; --lev) {

// 		int const crse_lev = lev-1;

// 		Update_DVDF_BarPlus (crse_lev, hDt(lev));
// 		Update_DVDF_BarPlus (lev, hDt(lev));
// 	}
// 		Update_DVDF_BarPlus (0, hDt(0));

// 	for (int lev = 0; lev <= finest_level; ++lev) {

// 		(this->*Update_DVDF_Grad_BarPlus) (lev);

// 		Update_DVDF_MesoFlux (lev, hDt(lev));
// 	}

// 	// CheckForPhysicalBoundary(0);

// 	AverageDownDVDFFlux();

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_DVDF_Wave (lev,hDt(lev));
// 	}

// 	AverageDownDVDFWave();

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_MacroVars (lev);
// 		Update_DVDF_Tau (lev);
// 	}

// 	AverageDownMacroVars();

// BL_PROFILE_VAR_STOP(Advance_SNEQ_DVDF);

// }

void DUGKSAmr::Advance_DVM()
{
	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {

		Update_DVDF_EqSo (lev);

		Collision_DVM_BarPlus (lev);

		if(!Ark::tinyEqual(Dt(), MinL())) { (this->*Update_DVDF_Grad_BarPlus) (lev); }

#ifdef ARK_USE_MOMENTUM
		// Update_DVDF_MesoFlux (lev, hDt(lev));
#endif
		Update_DVM_Wave (lev);

		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);

		//source terms
		auto &ld = Data(lev);
		Update_MacroSourceVars (ld.msv, 0.5*Dt());
	}
}
/*
void DUGKSAmr::Advance_Moment_SS()
{

	for (int lev = 0; lev <= finest_level; ++lev) {
		// incorporate source effects via Strang-splitting scheme
		Update_DVDF_Wave_SS(lev, 0.5*Dt());

		Nullify_Source (lev);
	}

//  regridding
	AverageDownMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Construct_DVDF_NonEq<false> (lev, 0);
	}

	if (m_regrid_int > 0 and nstep%m_regrid_int == 0)
    {
        if (m_verbose > 100) {
        	amrex::Print() << "Regridding...\n";
      	}

        regrid(0, cur_time);
    }

	// update value of macroscopic variables in ghost cells
	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Tau (lev);
		Update_DVDF_Eq (lev);
		Construct_DVDF_Wave<false> (lev, Dt());
	}

	for (int lev = finest_level; lev > 0; --lev) {

		int const crse_lev = lev-1;

		Collision_DVDF_BarPlus_Moment (crse_lev, hDt(lev));
		Collision_DVDF_BarPlus_Moment (lev, hDt(lev));
	}

		Collision_DVDF_BarPlus_Moment (0, hDt(0));

	for (int lev = 0; lev <= finest_level; ++lev) {

		(this->*Update_DVDF_Grad_BarPlus) (lev);

		Update_DVDF_MesoFlux_Moment<false> (lev, hDt(lev));
	}

	// CheckForPhysicalBoundary(0);

	AverageDownDVDFFlux();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Transform_DVDF_Moment (lev, Ark::BP);
		Construct_DVDF_WavePlus(lev, hDt(lev));
		Transform_Moment_DVDF(lev, Ark::WAVE);
		Update_DVDF_Wave_Lite (lev);
	}

	for (int lev = 0; lev <= finest_level; ++lev) {
		Transform_DVDF_Moment(lev,Ark::WAVE);
		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);
	}

	// FillPatchMacroVars();

	//Primitive distribution, stored in xWave
	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Eq (lev);
		Construct_DVDF_Prime<false> (lev, Dt());
	}

	// AverageDownMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Wave_SS (lev, 0.5*Dt());
	}
}
*/
void DUGKSAmr::Advance_Moment_SS()
{
	FillPatchMacroVars();

	// Primitive distribution, stored in xWave
	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Eq (lev);
		Construct_DVDF_Prime<false> (lev, Dt());
	}

	// for (int lev = 0; lev <= finest_level; ++lev) {
	// 	Update_DVDF_Wave_SS (lev, 0.5*Dt());
	// }

	// for (int lev = 0; lev <= finest_level; ++lev) {
	// 	// incorporate source effects via Strang-splitting scheme
	// 	Update_DVDF_Wave_SS(lev, 0.5*Dt());
	// 	Nullify_Source (lev);
	// }

//  regridding
	AverageDownMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		// 0: f-f^{E}
		// Dt: coeff*(\tilde{f} - f^{E})
		Construct_DVDF_NonEq<false> (lev, 0);
	}

	if (m_regrid_int > 0 and nstep%m_regrid_int == 0)
    {
        if (m_verbose > 100) {
        	amrex::Print() << "Regridding...\n";
      	}

        regrid(0, cur_time);
    }

	// update value of macroscopic variables in ghost cells
	FillPatchMacroVars();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_Tau (lev);
		Update_DVDF_Eq (lev);
		Construct_DVDF_Wave<false> (lev, Dt());
	}

	for (int lev = finest_level; lev > 0; --lev) {

		int const crse_lev = lev-1;

		Collision_DVDF_BarPlus_Moment (crse_lev, hDt(lev));
		Collision_DVDF_BarPlus_Moment (lev, hDt(lev));
	}
		Collision_DVDF_BarPlus_Moment (0, hDt(0));
		// Transform_Moment_DVDF(0, Ark::WAVE);
	for (int lev = 0; lev <= finest_level; ++lev) {

		(this->*Update_DVDF_Grad_BarPlus) (lev);

		// Update_DVDF_MesoFlux<false> (lev, hDt(lev));
		Update_DVDF_MesoFlux_Moment<false> (lev, hDt(lev));
	}

	// CheckForPhysicalBoundary(0);

	AverageDownDVDFFlux();

	for (int lev = 0; lev <= finest_level; ++lev) {
		Transform_DVDF_Moment (lev, Ark::BP);
		Construct_DVDF_WavePlus(lev, hDt(lev));
		Transform_Moment_DVDF(lev, Ark::WAVE);
		Update_DVDF_Wave_Lite (lev);
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
		Transform_DVDF_Moment(lev,Ark::WAVE);
		Update_MacroVars (lev);
		Update_DVDF_Tau (lev);
	}
}

// void DUGKSAmr::Advance_Moment()
// {
// //  regridding
// 	AverageDownMacroVars();

// 	FillPatchMacroVars();


// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_DVDF_Eq (lev);
// 	}

// 	for (int lev = finest_level; lev > 0; --lev) {

// 		int const crse_lev = lev-1;

// 		Collision_DVDF_BarPlus_Moment (crse_lev, hDt(lev));
// 		Collision_DVDF_BarPlus_Moment (lev, hDt(lev));
// 	}
// 		Collision_DVDF_BarPlus_Moment (0, hDt(0));
// 		// Transform_Moment_DVDF(0, Ark::WAVE);
// 	for (int lev = 0; lev <= finest_level; ++lev) {

// 		(this->*Update_DVDF_Grad_BarPlus) (lev);

// 		// Update_DVDF_MesoFlux<false> (lev, hDt(lev));
// 		Update_DVDF_MesoFlux_Moment<false> (lev, hDt(lev));
// 	}

// 	// CheckForPhysicalBoundary(0);

// 	AverageDownDVDFFlux();

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Transform_DVDF_Moment (lev, Ark::BP);
// 		Construct_DVDF_WavePlus(lev, hDt(lev));
// 		Transform_Moment_DVDF(lev, Ark::WAVE);
// 		Update_DVDF_Wave_Lite (lev);
// 	}

// 	AverageDownDVDFWave();

// 	if (m_regrid_int > 0 and nstep%m_regrid_int == 0)
//     {
//         if (m_verbose > 100) {
//         	amrex::Print() << "Regridding...\n";
//       	}

//         regrid(0, cur_time);
//     }

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Transform_DVDF_Moment(lev,Ark::WAVE);
// 		Update_MacroVars (lev);
// 		Update_DVDF_Tau (lev);
// 	}
// 	FillPatchMacroVars();

// 	// Primitive distribution, stored in xWave
// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_DVDF_Eq (lev);
// 		Construct_DVDF_Prime<false> (lev, Dt());
// 	}

// 	AverageDownMacroVars();

// 	for (int lev = 0; lev <= finest_level; ++lev) {
// 		Update_DVDF_Wave_SS (lev, 0.5*Dt());
// 	}
// }
void DUGKSAmr::Advance_Debug()
{
BL_PROFILE_VAR("DUGKSAmr::Advance_Debug_FillMacroVars()", Advance_Debug_FillMacroVars);

	AverageDownMacroVars();

	FillPatchMacroVars();

BL_PROFILE_VAR_STOP(Advance_CN_FillMacroVars);

BL_PROFILE_VAR("DUGKSAmr::Advance_CN_DVDF()", Advance_CN_DVDF);

	for (int lev = 0; lev <= finest_level; ++lev) {
		Update_DVDF_EqSo (lev);
	}

	for (int lev = finest_level; lev > 0; --lev) {

		int const crse_lev = lev-1;

		Collision_DVDF_BarPlus (crse_lev, hDt(lev));
		Collision_DVDF_BarPlus (lev, hDt(lev));
	}
		Collision_DVDF_BarPlus (0, hDt(0));


	for (int lev = 0; lev <= finest_level; ++lev) {

		(this->*Update_DVDF_Grad_BarPlus) (lev);
		Update_DVDF_MesoFlux<> (lev, hDt(lev));
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