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

DUGKSAmr::LevelData::LevelData (amrex::BoxArray const& ba,
                                amrex::DistributionMapping const& dm,
                                amrex::FabFactory<FArrayBox> const& fact,
                                int const ng)
    : phase_prev   (ba, dm, 1             , 1,  MFInfo(), fact),
      velocity_prev(ba, dm, AMREX_SPACEDIM, 1,  MFInfo(), fact),
      phase_old    (ba, dm, 1             , 1,  MFInfo(), fact),
      density_old  (ba, dm, 1             , 1,  MFInfo(), fact),
      velocity_old (ba, dm, AMREX_SPACEDIM, 1,  MFInfo(), fact),
      chemiP       (ba, dm, 1             , ng, MFInfo(), fact)
{
    msv[Ark::PHI].define    (ba, dm, 1             , ng, MFInfo(), fact);
    msv[Ark::RHO].define    (ba, dm, 1             , ng, MFInfo(), fact);
    msv[Ark::VEL].define    (ba, dm, AMREX_SPACEDIM, 1,  MFInfo(), fact);
    msv[Ark::P].define      (ba, dm, 1             , ng, MFInfo(), fact);
    msv[Ark::T].define      (ba, dm, 1             , ng, MFInfo(), fact);
    msv[Ark::GradPHI].define(ba, dm, AMREX_SPACEDIM, ng, MFInfo(), fact);
    msv[Ark::GradRHO].define(ba, dm, AMREX_SPACEDIM, ng, MFInfo(), fact);
    msv[Ark::FORCE].define  (ba, dm, AMREX_SPACEDIM, ng, MFInfo(), fact);

    BoxArray ba_xF = amrex::convert(ba, IntVect::TheDimensionVector(Ark::X));

    msv_xF[Ark::PHI].define    (ba_xF, dm, 1             , 0, MFInfo(), fact);
    msv_xF[Ark::RHO].define    (ba_xF, dm, 1             , 0, MFInfo(), fact);
    msv_xF[Ark::VEL].define    (ba_xF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_xF[Ark::P].define      (ba_xF, dm, 1             , 0, MFInfo(), fact);
    msv_xF[Ark::T].define      (ba_xF, dm, 1             , 0, MFInfo(), fact);
    msv_xF[Ark::GradPHI].define(ba_xF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_xF[Ark::GradRHO].define(ba_xF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_xF[Ark::FORCE].define  (ba_xF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);


    BoxArray ba_yF = amrex::convert(ba, IntVect::TheDimensionVector(Ark::Y));

    msv_yF[Ark::PHI].define    (ba_yF, dm, 1             , 0, MFInfo(), fact);
    msv_yF[Ark::RHO].define    (ba_yF, dm, 1             , 0, MFInfo(), fact);
    msv_yF[Ark::VEL].define    (ba_yF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_yF[Ark::P].define      (ba_yF, dm, 1             , 0, MFInfo(), fact);
    msv_yF[Ark::T].define      (ba_yF, dm, 1             , 0, MFInfo(), fact);
    msv_yF[Ark::GradPHI].define(ba_yF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_yF[Ark::GradRHO].define(ba_yF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_yF[Ark::FORCE].define  (ba_yF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);

#if(AMREX_SPACEDIM == 3)
    BoxArray ba_zF = amrex::convert(ba, IntVect::TheDimensionVector(Ark::Z));

    msv_zF[Ark::PHI].define    (ba_zF, dm, 1             , 0, MFInfo(), fact);
    msv_zF[Ark::RHO].define    (ba_zF, dm, 1             , 0, MFInfo(), fact);
    msv_zF[Ark::VEL].define    (ba_zF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_zF[Ark::P].define      (ba_zF, dm, 1             , 0, MFInfo(), fact);
    msv_zF[Ark::T].define      (ba_zF, dm, 1             , 0, MFInfo(), fact);
    msv_zF[Ark::GradPHI].define(ba_zF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_zF[Ark::GradRHO].define(ba_zF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
    msv_zF[Ark::FORCE].define  (ba_zF, dm, AMREX_SPACEDIM, 0, MFInfo(), fact);
#endif

#ifdef ARK_USE_PHASE

    hDVDF[Ark::TAU].define(ba, dm, 1, 0, MFInfo(), fact);
    for(int i = Ark::WAVE;i < Ark::DVDFNUM; ++i)
    {
      hDVDF[i].define(ba, dm, Ark::DVM_PHASE, ng, MFInfo(), fact);
    }

    hDVDFxF[Ark::TAUF].define(ba_xF,dm,1,0,MFInfo(),fact);
    for(int i = Ark::FLUX;i < Ark::FACEDVDFNUM; ++i)
    {
      hDVDFxF[i].define(ba_xF, dm, Ark::DVM_PHASE, 0, MFInfo(), fact);
    }

    hDVDFyF[Ark::TAUF].define(ba_yF,dm,1,0,MFInfo(),fact);
    for(int i = Ark::FLUX;i < Ark::FACEDVDFNUM; ++i)
    {
      hDVDFyF[i].define(ba_yF, dm, Ark::DVM_PHASE, 0, MFInfo(), fact);
    }

#if(AMREX_SPACEDIM == 3)
    hDVDFzF[Ark::TAUF].define(ba_zF,dm,1,0,MFInfo(),fact);
    for(int i = Ark::FLUX;i < Ark::FACEDVDFNUM; ++i)
    {
      hDVDFzF[i].define(ba_zF, dm, Ark::DVM_PHASE, 0, MFInfo(), fact);
    }
#endif

#endif

#ifdef ARK_USE_MOMENTUM

    fDVDF[Ark::TAU].define(ba, dm, 1, 0, MFInfo(), fact);
    for(int i = Ark::WAVE; i < Ark::DVDFNUM; ++i)
    {
      fDVDF[i].define(ba, dm, Ark::DVM_MOMENTUM, ng, MFInfo(), fact);
    }

    fDVDFxF[Ark::TAUF].define(ba_xF,dm,1,0,MFInfo(),fact);
    for(int i = Ark::FLUX; i < Ark::FACEDVDFNUM; ++i)
    {
      fDVDFxF[i].define(ba_xF, dm, Ark::DVM_MOMENTUM, 0, MFInfo(), fact);
    }

    fDVDFyF[Ark::TAUF].define(ba_yF,dm,1,0,MFInfo(),fact);
    for(int i = Ark::FLUX; i < Ark::FACEDVDFNUM; ++i)
    {
      fDVDFyF[i].define(ba_yF, dm, Ark::DVM_MOMENTUM, 0, MFInfo(), fact);
    }
#if(AMREX_SPACEDIM == 3)
    fDVDFzF[Ark::TAUF].define(ba_zF,dm,1,0,MFInfo(),fact);
    for(int i = Ark::FLUX; i < Ark::FACEDVDFNUM; ++i)
    {
      fDVDFzF[i].define(ba_zF, dm, Ark::DVM_MOMENTUM, 0, MFInfo(), fact);
    }
#endif

#endif

#ifdef ARK_USE_ENERGY

    gDVDF[Ark::TAU].define(ba,dm,1,0,MFInfo(),fact);
    for(int i = Ark::WAVE;i < Ark::DVDFNUM; ++i)
    {
      gDVDF[i].define(ba, dm, Ark::DVM_ENERGY, ng, MFInfo(), fact);
    }

    gDVDFxF[Ark::TAUF].define(ba_xF,dm,1,0,MFInfo(),fact);
    for(int i = Ark::FLUX; i < Ark::FACEDVDFNUM; ++i)
    {
      gDVDFxF[i].define(ba_xF, dm, Ark::DVM_ENERGY, 0, MFInfo(), fact);
    }

    gDVDFyF[Ark::TAUF].define(ba_yF,dm,1,0,MFInfo(),fact);
    for(int i = Ark::FLUX;i < Ark::FACEDVDFNUM; ++i)
    {
      gDVDFyF[i].define(ba_yF, dm, Ark::DVM_ENERGY, 0, MFInfo(), fact);
    }
#if(AMREX_SPACEDIM == 3)
    gDVDFzF[Ark::TAUF].define(ba_zF, dm, 1,0, MFInfo(), fact);
    for(int i = Ark::FLUX;i < Ark::FACEDVDFNUM; ++i)
    {
      gDVDFzF[i].define(ba_zF, dm, Ark::DVM_ENERGY, 0, MFInfo(), fact);
    }
#endif

#endif
}