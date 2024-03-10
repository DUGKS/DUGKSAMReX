#include <DUGKSAmr.H>

// DUGKSAmr::LevelData::LevelData (amrex::BoxArray const& ba,
//                                 amrex::DistributionMapping const& dm,
//                                 amrex::FabFactory<FArrayBox> const& fact,
//                                 int const ng)
//     : phase 	    (ba, dm, 1             , ng, MFInfo(), fact),
//       phase_old   (ba, dm, 1             , 1,  MFInfo(), fact),
//       density     (ba, dm, 1             , ng, MFInfo(), fact),
//       density_old (ba, dm, 1             , 1,  MFInfo(), fact),
//       velocity    (ba, dm, AMREX_SPACEDIM, ng, MFInfo(), fact),
//       velocity_old(ba, dm, AMREX_SPACEDIM, 1,  MFInfo(), fact),
//       pressure    (ba, dm, 1             , ng, MFInfo(), fact),
//       temperature (ba, dm, 1             , ng, MFInfo(), fact),
//       force       (ba, dm, AMREX_SPACEDIM, 1,  MFInfo(), fact),
//       grad_phase  (ba, dm, AMREX_SPACEDIM, 1,  MFInfo(), fact),
//       grad_density(ba, dm, AMREX_SPACEDIM, 1,  MFInfo(), fact)
// {
//     BoxArray ba_xF = amrex::convert(ba, IntVect::TheDimensionVector(Ark::X));
//     BoxArray ba_yF = amrex::convert(ba, IntVect::TheDimensionVector(Ark::Y));

//     phase_xF.define       (ba_xF, dm, 1             , 0, MFInfo(), fact);
//     density_xF.define     (ba_xF, dm, 1             , 0, MFInfo(), fact);
//     velocity_xF.define    (ba_xF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
//     pressure_xF.define    (ba_xF, dm, 1             , 0, MFInfo(), fact);
//     temperature_xF.define (ba_xF, dm, 1             , 0, MFInfo(), fact);
//     force_xF.define       (ba_xF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
//     gphase_xF.define      (ba_xF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
//     gdensity_xF.define    (ba_xF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);

//     phase_yF.define       (ba_yF, dm, 1             , 0, MFInfo(), fact);
//     density_yF.define     (ba_yF, dm, 1             , 0, MFInfo(), fact);
//     velocity_yF.define    (ba_yF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
//     pressure_yF.define    (ba_yF, dm, 1             , 0, MFInfo(), fact);
//     temperature_yF.define (ba_yF, dm, 1             , 0, MFInfo(), fact);
//     force_yF.define       (ba_yF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
//     gphase_yF.define      (ba_yF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
//     gdensity_yF.define    (ba_yF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
// }

void DUGKSAmr::LevelData::ConstructDVDFMultiFabAtFace (Array<MultiFab, Ark::FACEDVDFNUM>& xDVDFnF,
                                                  amrex::BoxArray const& ba_nF,
                                                  amrex::DistributionMapping const& dm,
                                                  amrex::FabFactory<FArrayBox> const& fact,
                                                  int const Q, int const ng)
{
    for(int i = Ark::TAUF;i < Ark::FACEDVDFNUM;++i)
    {
      xDVDFnF[i].define(ba_nF, dm, Q, 0, MFInfo(), fact);
    }
}

void DUGKSAmr::LevelData::ConstructDVDFMultiFabAtCell (Array<MultiFab, Ark::DVDFNUM>& xDVDF,
                                                  amrex::BoxArray const& ba,
                                                  amrex::DistributionMapping const& dm,
                                                  amrex::FabFactory<FArrayBox> const& fact,
                                                  int const Q, int const ng)
{
    //cell data
    xDVDF[Ark::TAU].define( ba, dm, Q,  0, MFInfo(), fact);
    xDVDF[Ark::WAVE].define(ba, dm, Q, ng, MFInfo(), fact);
    xDVDF[Ark::EQ].define(  ba, dm, Q,  0, MFInfo(), fact);
    xDVDF[Ark::SO].define(  ba, dm, Q,  0, MFInfo(), fact);
    xDVDF[Ark::BP].define(  ba, dm, Q, ng, MFInfo(), fact);

    // cell gradients
    for(int i = Ark::BP + 1;i < Ark::DVDFNUM; ++i)
    {
        xDVDF[i].define(ba, dm, Q, 1, MFInfo(), fact);
    }
}

void DUGKSAmr::LevelData::ConstructMSVMultiFabAtFace (Array<MultiFab, Ark::MSVNUM>& msv_nF,
                                      amrex::BoxArray const& ba_nF,
                                      amrex::DistributionMapping const& dm,
                                      amrex::FabFactory<FArrayBox> const& fact)
{
    msv_nF[Ark::PHI].define    (ba_nF, dm, 1             , 0, MFInfo(), fact);
    msv_nF[Ark::RHO].define    (ba_nF, dm, 1             , 0, MFInfo(), fact);
    msv_nF[Ark::VEL].define    (ba_nF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_nF[Ark::P].define      (ba_nF, dm, 1             , 0, MFInfo(), fact);
    msv_nF[Ark::T].define      (ba_nF, dm, 1             , 0, MFInfo(), fact);
    msv_nF[Ark::GradPHI].define(ba_nF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_nF[Ark::GradRHO].define(ba_nF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_nF[Ark::FORCE].define  (ba_nF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_nF[Ark::CHEMIP].define (ba_nF, dm, 1             , 0, MFInfo(), fact);
}

DUGKSAmr::LevelData::LevelData (amrex::BoxArray const& ba,
                                amrex::DistributionMapping const& dm,
                                amrex::FabFactory<FArrayBox> const& fact,
                                int const ng)
    : phase_prev   (ba, dm, 1             , 1, MFInfo(), fact),
      velocity_prev(ba, dm, AMREX_SPACEDIM, 1, MFInfo(), fact),
      phase_old    (ba, dm, 1             , 0, MFInfo(), fact),
      density_old  (ba, dm, 1             , 0, MFInfo(), fact),
      velocity_old (ba, dm, AMREX_SPACEDIM, 0, MFInfo(), fact)
{
    msv[Ark::PHI].define    (ba, dm, 1             , ng, MFInfo(), fact);
    msv[Ark::RHO].define    (ba, dm, 1             , ng, MFInfo(), fact);
    msv[Ark::VEL].define    (ba, dm, AMREX_SPACEDIM,  1, MFInfo(), fact);
    msv[Ark::P].define      (ba, dm, 1             , ng, MFInfo(), fact);
    msv[Ark::T].define      (ba, dm, 1             , ng, MFInfo(), fact);
    msv[Ark::GradPHI].define(ba, dm, AMREX_SPACEDIM, ng, MFInfo(), fact); //double interpolation
    msv[Ark::GradRHO].define(ba, dm, AMREX_SPACEDIM,  1, MFInfo(), fact);
    msv[Ark::FORCE].define  (ba, dm, AMREX_SPACEDIM,  1, MFInfo(), fact);
    msv[Ark::CHEMIP].define (ba, dm, 1             , ng, MFInfo(), fact);

    BoxArray ba_xF = amrex::convert(ba, IntVect::TheDimensionVector(Ark::X));
    ConstructMSVMultiFabAtFace(msv_xF, ba_xF, dm, fact);


    BoxArray ba_yF = amrex::convert(ba, IntVect::TheDimensionVector(Ark::Y));
    ConstructMSVMultiFabAtFace(msv_yF, ba_yF, dm, fact);

#if(AMREX_SPACEDIM == 3)
    BoxArray ba_zF = amrex::convert(ba, IntVect::TheDimensionVector(Ark::Z));
    ConstructMSVMultiFabAtFace(msv_zF, ba_zF, dm, fact);
#endif

#ifdef ARK_USE_PHASE

    ConstructDVDFMultiFabAtCell(hDVDF, ba, dm, fact, Ark::DVM_PHASE, ng);

    //face data
    ConstructDVDFMultiFabAtFace(hDVDFxF, ba_xF, dm, fact, Ark::DVM_PHASE, ng);
    ConstructDVDFMultiFabAtFace(hDVDFyF, ba_yF, dm, fact, Ark::DVM_PHASE, ng);
#if(AMREX_SPACEDIM == 3)
    ConstructDVDFMultiFabAtFace(hDVDFzF, ba_zF, dm, fact, Ark::DVM_PHASE, ng);
#endif

#endif

#ifdef ARK_USE_MOMENTUM

    ConstructDVDFMultiFabAtCell(fDVDF, ba, dm, fact, Ark::DVM_MOMENTUM, ng);

    //face data
    ConstructDVDFMultiFabAtFace(fDVDFxF, ba_xF, dm, fact, Ark::DVM_MOMENTUM, ng);
    ConstructDVDFMultiFabAtFace(fDVDFyF, ba_yF, dm, fact, Ark::DVM_MOMENTUM, ng);
#if(AMREX_SPACEDIM == 3)
    ConstructDVDFMultiFabAtFace(fDVDFzF, ba_zF, dm, fact, Ark::DVM_MOMENTUM, ng);
#endif

#endif

#ifdef ARK_USE_ENERGY

    ConstructDVDFMultiFabAtCell(gDVDF, ba, dm, fact, Ark::DVM_ENERGY, ng);

    //face data
    ConstructDVDFMultiFabAtFace(gDVDFxF, ba_xF, dm, fact, Ark::DVM_ENERGY, ng);
    ConstructDVDFMultiFabAtFace(gDVDFyF, ba_yF, dm, fact, Ark::DVM_ENERGY, ng);
#if(AMREX_SPACEDIM == 3)
    ConstructDVDFMultiFabAtFace(gDVDFzF, ba_zF, dm, fact, Ark::DVM_ENERGY, ng);
#endif

#endif
}