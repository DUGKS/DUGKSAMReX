#include <DUGKSAmr.H>
#include <DUGKSExternal.H>
#include <AMReX_FillPatchUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
#endif

using amrex::CpuBndryFuncFab;
using amrex::PhysBCFunct;
using amrex::Interpolater;

namespace {
// Interpolater* mapper = &amrex::cell_cons_interp;
Interpolater* mapper = &amrex::cell_bilinear_interp;
// Interpolater* mapper = &amrex::quadratic_interp;
// Interpolater* mapper = &amrex::quartic_interp;      //conservative_quadratic_interp
// Interpolater* mapper = &lincc_interp;
#ifdef AMREX_USE_EB
mapper = (EBFactory(0).isAllRegular()) ?
         (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#endif
}

void DUGKSAmr::FillPatchCore(int const lev, Real const time, MultiFab &mf, 
                             const Vector<MultiFab*>& cmf, const Vector<Real>& ctime,
                             const Vector<MultiFab*>& fmf, const Vector<Real>& ftime,
                             const Vector<BCRec>& bcs, amrex::UserFillBox FillExternalVariable)
{
    CpuBndryFuncFab bfunc(FillExternalVariable);

    if (lev == 0) {
        PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bfunc);
        FillPatchSingleLevel(mf, time, fmf, ftime, 0, 0, mf.nComp(), geom[lev], physbc, 0);
    } else {
        PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1], bcs, bfunc);
        PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev], bcs, bfunc);

        amrex::FillPatchTwoLevels(mf, time,
                                  cmf, ctime,
                                  fmf, ftime,
                                  0, 0, mf.nComp(), 
                                  geom[lev-1], geom[lev],
                                  cphysbc, 0, fphysbc, 0,
                                  refRatio(lev-1), mapper, bcs, 0);
    }
}
void DUGKSAmr::FillPatchPhase(int const lev, Real const time, MultiFab& phase)
{
	// if (lev == 0) {
	// 	CpuBndryFuncFab bfunc_phase(FillExternalPhase);
 //        PhysBCFunct<CpuBndryFuncFab> physbc_phase(geom[lev],get_phase_bcrec(),bfunc_phase);

 //        FillPatchSingleLevel(phi, IntVect(ng), time,
 //                             {&(m_leveldata[lev]->phase)},
 //                             {m_t_new[lev]}, 0, 0, phi.nComp(), geom[lev],
 //                             physbc_phase, 0);
 //    } else {
 //        const auto& bcrec = get_phase_bcrec();

 //        CpuBndryFuncFab bfunc_phase(FillExternalPhase);
 //        PhysBCFunct<CpuBndryFuncFab> cphysbc_phase(geom[lev-1],get_phase_bcrec(),bfunc_phase);
 //        PhysBCFunct<CpuBndryFuncFab> fphysbc_phase(geom[lev],get_phase_bcrec(),bfunc_phase);
        
 //        FillPatchTwoLevels(phi, IntVect(ng), time,
 //                           {&(m_leveldata[lev-1]->phase)}, {m_t_new[lev-1]},
 //                           {&(m_leveldata[lev]->phase)}, {m_t_new[lev]},
 //                           0, 0, 1, 
 //                           geom[lev-1], geom[lev],
 //                           cphysbc_phase, 0, fphysbc_phase, 0,
 //                           refRatio(lev-1), mapper, bcrec, 0);
 // }
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, phase,
                  {&(m_leveldata[crse_lev]->msv[Ark::PHI])}, {m_t_new[crse_lev]},
                  {&(m_leveldata[lev]->msv[Ark::PHI])}, {m_t_new[lev]},
                  get_phase_bcrec(), FillExternalPhase);
}
void DUGKSAmr::FillPatchDensity(int const lev, Real const time, MultiFab& density)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, density,
                  {&(m_leveldata[crse_lev]->msv[Ark::RHO])}, {m_t_new[crse_lev]},
                  {&(m_leveldata[lev]->msv[Ark::RHO])}, {m_t_new[lev]},
                  get_density_bcrec(), FillExternalDensity);
}
void DUGKSAmr::FillPatchVelocity(int const lev, Real const time, MultiFab& velocity)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, velocity,
                  {&(m_leveldata[crse_lev]->msv[Ark::VEL])}, {m_t_new[crse_lev]},
                  {&(m_leveldata[lev]->msv[Ark::VEL])}, {m_t_new[lev]},
                  get_velocity_bcrec(), FillExternalVelocity);
}
void DUGKSAmr::FillPatchPressure(int const lev, Real const time, MultiFab& pres)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, pres,
                  {&(m_leveldata[crse_lev]->msv[Ark::P])}, {m_t_new[crse_lev]},
                  {&(m_leveldata[lev]->msv[Ark::P])}, {m_t_new[lev]},
                  get_pressure_bcrec(), FillExternalPressure);
}
void DUGKSAmr::FillPatchTemperature(int const lev, Real const time, MultiFab& temp)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, temp,
                  {&(m_leveldata[crse_lev]->msv[Ark::T])}, {m_t_new[crse_lev]},
                  {&(m_leveldata[lev]->msv[Ark::T])}, {m_t_new[lev]},
                  get_temperature_bcrec(), FillExternalTemperature);
}
void DUGKSAmr::FillPatchGradPhase(int const lev, Real const time, MultiFab& grad_phi)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, grad_phi,
                  {&(m_leveldata[crse_lev]->msv[Ark::GradPHI])}, {m_t_new[crse_lev]},
                  {&(m_leveldata[lev]->msv[Ark::GradPHI])}, {m_t_new[lev]},
                  get_force_bcrec(), FillExternalForce);
}
void DUGKSAmr::FillPatchGradDensity(int const lev, Real const time, MultiFab& grad_rho)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, grad_rho,
                  {&(m_leveldata[crse_lev]->msv[Ark::GradRHO])}, {m_t_new[crse_lev]},
                  {&(m_leveldata[lev]->msv[Ark::GradRHO])}, {m_t_new[lev]},
                  get_force_bcrec(), FillExternalForce);
}
void DUGKSAmr::FillPatchForce(int const lev, Real const time, MultiFab& force)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, force,
                  {&(m_leveldata[crse_lev]->msv[Ark::FORCE])}, {m_t_new[crse_lev]},
                  {&(m_leveldata[lev]->msv[Ark::FORCE])}, {m_t_new[lev]},
                  get_force_bcrec(), FillExternalForce);
}
void DUGKSAmr::FillPatchDVDFWave (int const lev, Real const time, MultiFab& xWave, 
                                  Array<MultiFab, Ark::DVDFNUM> LevelData::*xDVDF)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, xWave,
                  {&(((*m_leveldata[crse_lev]).*xDVDF)[Ark::WAVE])}, {m_t_new[crse_lev]},
                  {&(((*m_leveldata[lev]).*xDVDF)[Ark::WAVE])}, {m_t_new[lev]},
                  get_DVDF_wave_bcrec(), FillExternalDVDFWave);
}
void DUGKSAmr::FillPatchDVDFBarPlus (int const lev, Real const time, MultiFab& xBP, 
                                     Array<MultiFab, Ark::DVDFNUM> LevelData::*xDVDF)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore(lev, time, xBP,
                  {&( ((*m_leveldata[crse_lev]).*xDVDF)[Ark::BP] )}, {m_t_new[crse_lev]},
                  {&( ((*m_leveldata[lev]).*xDVDF)[Ark::BP] )}, {m_t_new[lev]},
                  get_DVDF_barplus_bcrec(), FillExternalDVDFBarPlus);
}
void DUGKSAmr::FillPatchMacroVars()
{
    for(int lev = 0; lev <= finest_level; ++lev)
    {
        auto &ld = Data(lev);

        FillPatchPhase      (lev, cur_time, ld.msv[Ark::PHI]);
        FillPatchDensity    (lev, cur_time, ld.msv[Ark::RHO]);
        FillPatchVelocity   (lev, cur_time, ld.msv[Ark::VEL]);
        FillPatchPressure   (lev, cur_time, ld.msv[Ark::P]);
        FillPatchTemperature(lev, cur_time, ld.msv[Ark::T]);
        FillPatchGradPhase  (lev, cur_time, ld.msv[Ark::GradPHI]);
        FillPatchGradDensity(lev, cur_time, ld.msv[Ark::GradRHO]);
        FillPatchForce      (lev, cur_time, ld.msv[Ark::FORCE]);
    }
}
void DUGKSAmr::FillPatchNew(int const lev, Real const time, DUGKSAmr::LevelData &newld)
{
    FillPatchPhase      (lev, time, newld.msv[Ark::PHI]);
    FillPatchDensity    (lev, time, newld.msv[Ark::RHO]);
    FillPatchVelocity   (lev, time, newld.msv[Ark::VEL]);
    FillPatchPressure   (lev, time, newld.msv[Ark::P]);
    FillPatchTemperature(lev, time, newld.msv[Ark::T]);

#ifdef ARK_USE_PHASE
    FillPatchDVDFWave(lev, time, newld.hDVDF[Ark::WAVE], &LevelData::hDVDF);
    // FillPatchDVDFBarPlus(lev, time, newld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
    FillPatchDVDFWave(lev, time, newld.fDVDF[Ark::WAVE], &LevelData::fDVDF);
    // FillPatchDVDFBarPlus(lev, time, newld.fDVDF[Ark::BP], &LevelData::fDVDF);
#endif

#ifdef ARK_USE_ENERGY
    FillPatchDVDFWave(lev, time, newld.gDVDF[Ark::WAVE], &LevelData::gDVDF);
    // FillPatchDVDFBarPlus(lev, time, newld.gDVDF[Ark::BP], &LevelData::gDVDF);
#endif
}
/*fill the fine lev from the coarse lev*/
void DUGKSAmr::FillCoarsePatchCore(int const lev, Real const time, 
                                   MultiFab &mf, MultiFab const &cmf,
                                   const Vector<BCRec> &bcs, 
                                   amrex::UserFillBox FillExternalVariable)
{
    CpuBndryFuncFab bfunc(FillExternalVariable);
    PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
    PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bfunc);

    amrex::InterpFromCoarseLevel(mf, time, cmf, 0, 0, mf.nComp(),
                                 geom[lev-1], geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 refRatio(lev-1), mapper, bcs, 0);
}
void DUGKSAmr::FillCoarsePatchPhase(int const lev, Real const time, MultiFab& phase)
{
	AMREX_ASSERT(lev > 0);
  FillCoarsePatchCore(lev, time, phase, m_leveldata[lev-1]->msv[Ark::PHI],
                        get_phase_bcrec(), FillExternalPhase);
}
void DUGKSAmr::FillCoarsePatchDensity(int const lev, Real const time, MultiFab& density) 
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore(lev, time, density, m_leveldata[lev-1]->msv[Ark::RHO],
                        get_density_bcrec(), FillExternalDensity);
}

void DUGKSAmr::FillCoarsePatchVelocity(int const lev, Real const time, MultiFab& velocity)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore(lev, time, velocity, m_leveldata[lev-1]->msv[Ark::VEL],
                        get_velocity_bcrec(), FillExternalVelocity);
}
void DUGKSAmr::FillCoarsePatchPressure(int const lev, Real const time, MultiFab& pres)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore(lev, time, pres, m_leveldata[lev-1]->msv[Ark::P],
                        get_pressure_bcrec(), FillExternalPressure);
}
void DUGKSAmr::FillCoarsePatchTemperature(int const lev, Real const time, MultiFab& temp)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore(lev, time, temp, m_leveldata[lev-1]->msv[Ark::T],
                        get_temperature_bcrec(), FillExternalTemperature);
}
void DUGKSAmr::FillCoarsePatchDVDFWave (int const lev, Real const time, 
                                        MultiFab& xWave, 
                                        Array<MultiFab, Ark::DVDFNUM> LevelData:: *xDVDF)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore(lev, time, xWave, ((*m_leveldata[lev-1]).*xDVDF)[Ark::WAVE],
                        get_DVDF_wave_bcrec(), FillExternalDVDFWave);
}
void DUGKSAmr::FillCoarsePatchDVDFBarPlus (int const lev, Real const time, 
                                           MultiFab& xBP, 
                                           Array<MultiFab, Ark::DVDFNUM> LevelData:: *xDVDF)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore(lev, time, xBP, ((*m_leveldata[lev-1]).*xDVDF)[Ark::BP],
                        get_DVDF_barplus_bcrec(), FillExternalDVDFBarPlus);
}

void DUGKSAmr::FillCoarsePatchNew(int const lev, Real const time, DUGKSAmr::LevelData &newld)
{
	FillCoarsePatchPhase(lev, time, newld.msv[Ark::PHI]);
	FillCoarsePatchDensity(lev, time, newld.msv[Ark::RHO]);
  FillCoarsePatchVelocity(lev, time, newld.msv[Ark::VEL]);
  FillCoarsePatchPressure(lev, time, newld.msv[Ark::P]);
  FillCoarsePatchTemperature(lev, time, newld.msv[Ark::T]);

#ifdef ARK_USE_PHASE
    FillCoarsePatchDVDFWave(lev, time, newld.hDVDF[Ark::WAVE], &LevelData::hDVDF);
    // FillCoarsePatchDVDFBarPlus(lev, time, newld.hDVDF[Ark::BP], &LevelData::hDVDF);
#endif

#ifdef ARK_USE_MOMENTUM
    FillCoarsePatchDVDFWave(lev, time, newld.fDVDF[Ark::WAVE], &LevelData::fDVDF);
    // FillCoarsePatchDVDFBarPlus(lev, time, newld.fDVDF[Ark::BP], &LevelData::fDVDF);
 #endif

#ifdef ARK_USE_ENERGY
    FillCoarsePatchDVDFWave(lev, time, newld.gDVDF[Ark::WAVE], &LevelData::gDVDF);
    // FillCoarsePatchDVDFBarPlus(lev, cur_time, newld.gDVDF[Ark::WAVE], &LevelData::gDVDF);
#endif
}

void DUGKSAmr::FillInteriorBoundary (MultiFab &xFab, const Periodicity& period)
{
  xFab.FillBoundary(period);
}