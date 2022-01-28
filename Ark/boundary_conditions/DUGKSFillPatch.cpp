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

#ifdef AMREX_USE_EB
    Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                        (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
// Interpolater* mapper = &amrex::cell_cons_interp;
// Interpolater* mapper = &amrex::cell_bilinear_interp;
Interpolater* mapper = &amrex::quadratic_interp;
// Interpolater* mapper = &amrex::quartic_interp;      //conservative_quadratic_interp
// Interpolater* mapper = &lincc_interp;
#endif
}

// void DUGKSAmr::FillPatchCore(int const lev, Real const time, MultiFab &mf, 
//                              const Vector<MultiFab*>& cmf, const Vector<Real>& ctime,
//                              const Vector<MultiFab*>& fmf, const Vector<Real>& ftime,
//                              const Vector<BCRec>& bcs, amrex::UserFillBox FillExternalVariable)
// {
//     CpuBndryFuncFab bfunc(FillExternalVariable);

//     if (lev == 0) {
//         PhysBCFunct<CpuBndryFuncFab> physbc(geom[lev],bcs,bfunc);
//         FillPatchSingleLevel(mf, time, fmf, ftime, 0, 0, mf.nComp(), geom[lev], physbc, 0);
//     } else {
//         PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1], bcs, bfunc);
//         PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev], bcs, bfunc);

//         amrex::FillPatchTwoLevels(mf, time,
//                                   cmf, ctime,
//                                   fmf, ftime,
//                                   0, 0, mf.nComp(), 
//                                   geom[lev-1], geom[lev],
//                                   cphysbc, 0, fphysbc, 0,
//                                   refRatio(lev-1), mapper, bcs, 0);
//     }
// }

template <typename FillExternal>
void DUGKSAmr::FillPatchCore(int const lev, Real const time, MultiFab &mf, 
                             const Vector<MultiFab*>& cmf, const Vector<Real>& ctime,
                             const Vector<MultiFab*>& fmf, const Vector<Real>& ftime,
                             const Vector<BCRec>& bcs)
{
    amrex::GpuBndryFuncFab<FillExternal> bfunc(FillExternal{});

    if (lev == 0) {

        PhysBCFunct<amrex::GpuBndryFuncFab<FillExternal>> physbc(geom[lev],bcs,bfunc);
        FillPatchSingleLevel(mf, time, fmf, ftime, 0, 0, mf.nComp(), geom[lev], physbc, 0);
    } else {

        PhysBCFunct<amrex::GpuBndryFuncFab<FillExternal>> cphysbc(geom[lev-1], bcs, bfunc);
        PhysBCFunct<amrex::GpuBndryFuncFab<FillExternal>> fphysbc(geom[lev], bcs, bfunc);

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
 //                             {time}, 0, 0, phi.nComp(), geom[lev],
 //                             physbc_phase, 0);
 //    } else {
 //        const auto& bcrec = get_phase_bcrec();

 //        CpuBndryFuncFab bfunc_phase(FillExternalPhase);
 //        PhysBCFunct<CpuBndryFuncFab> cphysbc_phase(geom[lev-1],get_phase_bcrec(),bfunc_phase);
 //        PhysBCFunct<CpuBndryFuncFab> fphysbc_phase(geom[lev],get_phase_bcrec(),bfunc_phase);
        
 //        FillPatchTwoLevels(phi, IntVect(ng), time,
 //                           {&(m_leveldata[lev-1]->phase)}, {m_t_new[lev-1]},
 //                           {&(m_leveldata[lev]->phase)}, {time},
 //                           0, 0, 1, 
 //                           geom[lev-1], geom[lev],
 //                           cphysbc_phase, 0, fphysbc_phase, 0,
 //                           refRatio(lev-1), mapper, bcrec, 0);
 // }
    int const crse_lev = (lev == 0 ? lev : lev-1);
    // FillPatchCore(lev, time, phase,
    //               {&(m_leveldata[crse_lev]->msv[Ark::PHI])}, {time},
    //               {&(m_leveldata[lev]->msv[Ark::PHI])}, {time},
    //               get_phase_bcrec(), FillExternalPhase);

    FillPatchCore<FillExternalPhase>(lev, time, phase,
                  {&(m_leveldata[crse_lev]->msv[Ark::PHI])}, {time},
                  {&(m_leveldata[lev]->msv[Ark::PHI])}, {time},
                  get_phase_bcrec());
}
void DUGKSAmr::FillPatchPhasePrev(int const lev, Real const time, MultiFab& phasePrev)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);

    FillPatchCore<FillExternalPhase>(lev, time, phasePrev,
                  {&(m_leveldata[crse_lev]->phase_prev)}, {time},
                  {&(m_leveldata[lev]->phase_prev)}, {time},
                  get_phase_bcrec());
}
void DUGKSAmr::FillPatchDensity(int const lev, Real const time, MultiFab& density)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalDensity>(lev, time, density,
                  {&(m_leveldata[crse_lev]->msv[Ark::RHO])}, {time},
                  {&(m_leveldata[lev]->msv[Ark::RHO])}, {time},
                  get_density_bcrec());
}
void DUGKSAmr::FillPatchVelocity(int const lev, Real const time, MultiFab& velocity)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalVelocity>(lev, time, velocity,
                  {&(m_leveldata[crse_lev]->msv[Ark::VEL])}, {time},
                  {&(m_leveldata[lev]->msv[Ark::VEL])}, {time},
                  get_velocity_bcrec());
}
void DUGKSAmr::FillPatchVelocityPrev(int const lev, Real const time, MultiFab& velocityPrev)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalVelocity>(lev, time, velocityPrev,
                  {&(m_leveldata[crse_lev]->velocity_prev)}, {time},
                  {&(m_leveldata[lev]->velocity_prev)}, {time},
                  get_velocity_bcrec());
}
void DUGKSAmr::FillPatchPressure(int const lev, Real const time, MultiFab& pres)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalPressure>(lev, time, pres,
                  {&(m_leveldata[crse_lev]->msv[Ark::P])}, {time},
                  {&(m_leveldata[lev]->msv[Ark::P])}, {time},
                  get_pressure_bcrec());
}
void DUGKSAmr::FillPatchTemperature(int const lev, Real const time, MultiFab& temp)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalTemperature>(lev, time, temp,
                  {&(m_leveldata[crse_lev]->msv[Ark::T])}, {time},
                  {&(m_leveldata[lev]->msv[Ark::T])}, {time},
                  get_temperature_bcrec());
}
void DUGKSAmr::FillPatchChemiP(int const lev, Real const time, MultiFab& chemiP)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalEmpty>(lev, time, chemiP,
                  {&(m_leveldata[crse_lev]->chemiP)}, {time},
                  {&(m_leveldata[lev]->chemiP)}, {time},
                  get_force_bcrec());
}
void DUGKSAmr::FillPatchGradPhase(int const lev, Real const time, MultiFab& grad_phi)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalForce>(lev, time, grad_phi,
                  {&(m_leveldata[crse_lev]->msv[Ark::GradPHI])}, {time},
                  {&(m_leveldata[lev]->msv[Ark::GradPHI])}, {time},
                  get_force_bcrec());
}
void DUGKSAmr::FillPatchGradDensity(int const lev, Real const time, MultiFab& grad_rho)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalForce>(lev, time, grad_rho,
                  {&(m_leveldata[crse_lev]->msv[Ark::GradRHO])}, {time},
                  {&(m_leveldata[lev]->msv[Ark::GradRHO])}, {time},
                  get_force_bcrec());
}
void DUGKSAmr::FillPatchForce(int const lev, Real const time, MultiFab& force)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalForce>(lev, time, force,
                  {&(m_leveldata[crse_lev]->msv[Ark::FORCE])}, {time},
                  {&(m_leveldata[lev]->msv[Ark::FORCE])}, {time},
                  get_force_bcrec());
}
void DUGKSAmr::FillPatchDVDFWave (int const lev, Real const time, MultiFab& xWave, 
                                  Array<MultiFab, Ark::DVDFNUM> LevelData::*xDVDF)
{
    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalDVDFWave>(lev, time, xWave,
                  {&(((*m_leveldata[crse_lev]).*xDVDF)[Ark::WAVE])}, {time},
                  {&(((*m_leveldata[lev]).*xDVDF)[Ark::WAVE])}, {time},
                  get_DVDF_wave_bcrec());
}
void DUGKSAmr::FillPatchDVDFBarPlus (int const lev, Real const time, MultiFab& xBP, 
                                     Array<MultiFab, Ark::DVDFNUM> LevelData::*xDVDF)
{
BL_PROFILE_VAR("DUGKSAmr::FillPatchDVDFBarPlus()", FillPatchDVDFBarPlus);

    int const crse_lev = (lev == 0 ? lev : lev-1);
    FillPatchCore<FillExternalDVDFBarPlus>(lev, time, xBP,
                  {&( ((*m_leveldata[crse_lev]).*xDVDF)[Ark::BP] )}, {time},
                  {&( ((*m_leveldata[lev]).*xDVDF)[Ark::BP] )}, {time},
                  get_DVDF_barplus_bcrec());
}
void DUGKSAmr::FillPatchMacroVars()
{

BL_PROFILE_VAR("DUGKSAmr::FillPatchMacroVars()", FillPatchMacroVars);

    for(int lev = 0; lev <= finest_level; ++lev)
    {
        auto &ld = Data(lev);
#ifdef ARK_USE_PHASE
        FillPatchPhase      (lev, cur_time, ld.msv[Ark::PHI]);
        FillPatchGradPhase  (lev, cur_time, ld.msv[Ark::GradPHI]);
#endif
        FillPatchVelocity   (lev, cur_time, ld.msv[Ark::VEL]);
        FillPatchPhasePrev    (lev, cur_time, ld.phase_prev);
        FillPatchVelocityPrev (lev, cur_time, ld.velocity_prev);

#ifdef ARK_USE_MOMENTUM
        FillPatchDensity    (lev, cur_time, ld.msv[Ark::RHO]);
        FillPatchPressure   (lev, cur_time, ld.msv[Ark::P]);
        FillPatchTemperature(lev, cur_time, ld.msv[Ark::T]);
        FillPatchGradDensity(lev, cur_time, ld.msv[Ark::GradRHO]);
        FillPatchForce      (lev, cur_time, ld.msv[Ark::FORCE]);
#endif
    }
}
void DUGKSAmr::FillPatchNew(int const lev, Real const time, DUGKSAmr::LevelData &newld)
{
#ifdef ARK_USE_PHASE
        FillPatchPhase      (lev, cur_time, newld.msv[Ark::PHI]);
        FillPatchGradPhase  (lev, cur_time, newld.msv[Ark::GradPHI]);
#endif

        FillPatchVelocity     (lev, cur_time, newld.msv[Ark::VEL]);
        FillPatchPhasePrev    (lev, cur_time, newld.phase_prev);
        FillPatchVelocityPrev (lev, cur_time, newld.velocity_prev);
        FillPatchForce        (lev, cur_time, newld.msv[Ark::FORCE]);

#ifdef ARK_USE_MOMENTUM
        FillPatchDensity      (lev, cur_time, newld.msv[Ark::RHO]);
        FillPatchPressure     (lev, cur_time, newld.msv[Ark::P]);
        FillPatchTemperature  (lev, cur_time, newld.msv[Ark::T]);
        FillPatchGradDensity  (lev, cur_time, newld.msv[Ark::GradRHO]);
#endif

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
void DUGKSAmr::PostFillPatchNew (int const lev, DUGKSAmr::LevelData &newld)
{

}
/*fill the fine lev from the coarse lev*/

// void DUGKSAmr::FillCoarsePatchCore(int const lev, Real const time, 
//                                    MultiFab &mf, MultiFab const &cmf,
//                                    const Vector<BCRec> &bcs, 
//                                    amrex::UserFillBox FillExternalVariable)
// {
//     CpuBndryFuncFab bfunc(FillExternalVariable);
//     PhysBCFunct<CpuBndryFuncFab> cphysbc(geom[lev-1],bcs,bfunc);
//     PhysBCFunct<CpuBndryFuncFab> fphysbc(geom[lev],bcs,bfunc);

//     amrex::InterpFromCoarseLevel(mf, time, cmf, 0, 0, mf.nComp(),
//                                  geom[lev-1], geom[lev],
//                                  cphysbc, 0, fphysbc, 0,
//                                  refRatio(lev-1), mapper, bcs, 0);
// }

template <typename FillExternal>
void DUGKSAmr::FillCoarsePatchCore(int const lev, Real const time, 
                                   MultiFab &mf, MultiFab const &cmf,
                                   const Vector<BCRec> &bcs)
{
    amrex::GpuBndryFuncFab<FillExternal> bfunc(FillExternal{});
    PhysBCFunct<amrex::GpuBndryFuncFab<FillExternal>> cphysbc(geom[lev-1], bcs, bfunc);
    PhysBCFunct<amrex::GpuBndryFuncFab<FillExternal>> fphysbc(geom[lev  ], bcs, bfunc);

    amrex::InterpFromCoarseLevel(mf, time, cmf, 0, 0, mf.nComp(),
                                 geom[lev-1], geom[lev],
                                 cphysbc, 0, fphysbc, 0,
                                 refRatio(lev-1), mapper, bcs, 0);
}
void DUGKSAmr::FillCoarsePatchPhase(int const lev, Real const time, MultiFab& phase)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore<FillExternalPhase> (lev, time, phase, 
                                            m_leveldata[lev-1]->msv[Ark::PHI],
                                            get_phase_bcrec());
}
void DUGKSAmr::FillCoarsePatchDensity(int const lev, Real const time, MultiFab& density) 
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore<FillExternalDensity>(lev, time, density,
                                             m_leveldata[lev-1]->msv[Ark::RHO],
                                             get_density_bcrec());
}

void DUGKSAmr::FillCoarsePatchVelocity(int const lev, Real const time, MultiFab& velocity)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore<FillExternalVelocity>(lev, time, velocity, 
                                              m_leveldata[lev-1]->msv[Ark::VEL],
                                              get_velocity_bcrec());
}
void DUGKSAmr::FillCoarsePatchPressure(int const lev, Real const time, MultiFab& pres)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore<FillExternalPressure>(lev, time, pres, 
                                              m_leveldata[lev-1]->msv[Ark::P],
                                              get_pressure_bcrec());
}
void DUGKSAmr::FillCoarsePatchTemperature(int const lev, Real const time, MultiFab& temp)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore<FillExternalTemperature>(lev, time, temp, 
                                                 m_leveldata[lev-1]->msv[Ark::T],
                                                 get_temperature_bcrec());
}
void DUGKSAmr::FillCoarsePatchDVDFWave (int const lev, Real const time, 
                                        MultiFab& xWave, 
                                        Array<MultiFab, Ark::DVDFNUM> LevelData:: *xDVDF)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore<FillExternalDVDFWave>(lev, time, xWave, 
                                              ((*m_leveldata[lev-1]).*xDVDF)[Ark::WAVE],
                                              get_DVDF_wave_bcrec());
}
void DUGKSAmr::FillCoarsePatchDVDFBarPlus (int const lev, Real const time, 
                                           MultiFab& xBP, 
                                           Array<MultiFab, Ark::DVDFNUM> LevelData:: *xDVDF)
{
    AMREX_ASSERT(lev > 0);
    FillCoarsePatchCore<FillExternalDVDFBarPlus>(lev, time, xBP,
                                                 ((*m_leveldata[lev-1]).*xDVDF)[Ark::BP],
                                                 get_DVDF_barplus_bcrec());
}

void DUGKSAmr::FillCoarsePatchNew(int const lev, Real const time, DUGKSAmr::LevelData &newld)
{
#ifdef ARK_USE_PHASE
        FillPatchPhase      (lev, cur_time, newld.msv[Ark::PHI]);
        FillPatchGradPhase  (lev, cur_time, newld.msv[Ark::GradPHI]);
#endif

        FillPatchVelocity     (lev, cur_time, newld.msv[Ark::VEL]);
        FillPatchPhasePrev    (lev, cur_time, newld.phase_prev);
        FillPatchVelocityPrev (lev, cur_time, newld.velocity_prev);

#ifdef ARK_USE_MOMENTUM
        FillPatchDensity      (lev, cur_time, newld.msv[Ark::RHO]);
        FillPatchPressure     (lev, cur_time, newld.msv[Ark::P]);
        FillPatchTemperature  (lev, cur_time, newld.msv[Ark::T]);
        FillPatchGradDensity  (lev, cur_time, newld.msv[Ark::GradRHO]);
        FillPatchForce        (lev, cur_time, newld.msv[Ark::FORCE]);
#endif

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