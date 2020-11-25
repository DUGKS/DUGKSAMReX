#include <DUGKSAmr.H>
#include <AMReX_FillPatchUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
#endif

void DUGKSAmr::AverageDownLevelAll(int const crse_lev)
{
	AverageDownLevelPhase(crse_lev);
	AverageDownLevelDensity(crse_lev);
	AverageDownLevelVelocity(crse_lev);
	AverageDownLevelPressure(crse_lev);
	AverageDownLevelTemperature(crse_lev);
	AverageDownLevelGradPhase(crse_lev);
	AverageDownLevelGradDensity(crse_lev);
	AverageDownLevelForce(crse_lev);
	AverageDownLevelDVDFWave(crse_lev);
}
void DUGKSAmr::AverageDownLevelMacroVars(int const crse_lev)
{
	AverageDownLevelPhase(crse_lev);
	AverageDownLevelDensity(crse_lev);
	AverageDownLevelVelocity(crse_lev);
	AverageDownLevelPressure(crse_lev);
	AverageDownLevelTemperature(crse_lev);
	// AverageDownLevelGradPhase(crse_lev);
	// AverageDownLevelGradDensity(crse_lev);
	// AverageDownLevelForce(crse_lev);
}
void DUGKSAmr::AverageDownLevelPhase(int const crse_lev)
{
	amrex::average_down(m_leveldata[crse_lev+1]->msv[Ark::PHI],
						m_leveldata[crse_lev]->msv[Ark::PHI],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, 1, refRatio(crse_lev));
}
void DUGKSAmr::AverageDownLevelDensity(int const crse_lev)
{
	amrex::average_down(m_leveldata[crse_lev+1]->msv[Ark::RHO],
						m_leveldata[crse_lev]->msv[Ark::RHO],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, 1, refRatio(crse_lev));
}
void DUGKSAmr::AverageDownLevelVelocity(int const crse_lev)
{
	amrex::average_down(m_leveldata[crse_lev+1]->msv[Ark::VEL], 
						m_leveldata[crse_lev]->msv[Ark::VEL],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, AMREX_SPACEDIM, refRatio(crse_lev));
}
void DUGKSAmr::AverageDownLevelPressure(int const crse_lev)
{
	amrex::average_down(m_leveldata[crse_lev+1]->msv[Ark::P],
						m_leveldata[crse_lev]->msv[Ark::P],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, 1, refRatio(crse_lev));
}
void DUGKSAmr::AverageDownLevelTemperature(int const crse_lev)
{
	amrex::average_down(m_leveldata[crse_lev+1]->msv[Ark::T],
						m_leveldata[crse_lev]->msv[Ark::T],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, 1, refRatio(crse_lev));
}
void DUGKSAmr::AverageDownLevelGradPhase(int const crse_lev)
{
	amrex::average_down(m_leveldata[crse_lev+1]->msv[Ark::GradPHI], 
						m_leveldata[crse_lev]->msv[Ark::GradPHI],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, AMREX_SPACEDIM, refRatio(crse_lev));
}
void DUGKSAmr::AverageDownLevelGradDensity(int const crse_lev)
{
	amrex::average_down(m_leveldata[crse_lev+1]->msv[Ark::GradRHO], 
						m_leveldata[crse_lev]->msv[Ark::GradRHO],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, AMREX_SPACEDIM, refRatio(crse_lev));
}
void DUGKSAmr::AverageDownLevelForce(int const crse_lev)
{
	amrex::average_down(m_leveldata[crse_lev+1]->msv[Ark::FORCE], 
						m_leveldata[crse_lev]->msv[Ark::FORCE],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, AMREX_SPACEDIM, refRatio(crse_lev));
}
void DUGKSAmr::AverageDownLevelDVDFWave(int const crse_lev)
{
#ifdef ARK_USE_PHASE
	amrex::average_down(m_leveldata[crse_lev+1]->hDVDF[Ark::WAVE],
			            m_leveldata[crse_lev]->hDVDF[Ark::WAVE],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, Ark::DVM_PHASE, refRatio(crse_lev));
#endif

#ifdef ARK_USE_MOMENTUM
	amrex::average_down(m_leveldata[crse_lev+1]->fDVDF[Ark::WAVE],
			            m_leveldata[crse_lev]->fDVDF[Ark::WAVE],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, Ark::DVM_MOMENTUM, refRatio(crse_lev));
#endif

#ifdef ARK_USE_ENERGY
	amrex::average_down(m_leveldata[crse_lev+1]->gDVDF[Ark::WAVE],
			            m_leveldata[crse_lev]->gDVDF[Ark::WAVE],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, Ark::DVM_MOMENTUM, refRatio(crse_lev));
#endif
}

void DUGKSAmr::AverageDownLevelDVDFBarPlus(int const crse_lev)
{
#ifdef ARK_USE_PHASE
	amrex::average_down(m_leveldata[crse_lev+1]->hDVDF[Ark::BP],
			            m_leveldata[crse_lev]->hDVDF[Ark::BP],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, Ark::DVM_PHASE, refRatio(crse_lev));
#endif

#ifdef ARK_USE_MOMENTUM
	amrex::average_down(m_leveldata[crse_lev+1]->fDVDF[Ark::BP],
			            m_leveldata[crse_lev]->fDVDF[Ark::BP],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, Ark::DVM_MOMENTUM, refRatio(crse_lev));
#endif

#ifdef ARK_USE_ENERGY
	amrex::average_down(m_leveldata[crse_lev+1]->gDVDF[Ark::BP],
			            m_leveldata[crse_lev]->gDVDF[Ark::BP],
                        // geom[crse_lev+1], geom[crse_lev],
                        0, Ark::DVM_ENERGY, refRatio(crse_lev));
#endif
}
void DUGKSAmr::AverageDownLevelDVDFFlux(int const crse_lev)
{
#ifdef ARK_USE_PHASE
	amrex::average_down_faces(m_leveldata[crse_lev+1]->hDVDFxF[Ark::FLUX],
			            	  m_leveldata[crse_lev]->hDVDFxF[Ark::FLUX],
                        	  refRatio(crse_lev),0);
	//
	amrex::average_down_faces(m_leveldata[crse_lev+1]->hDVDFyF[Ark::FLUX],
			            	  m_leveldata[crse_lev]->hDVDFyF[Ark::FLUX],
                        	  refRatio(crse_lev),0);

#if(AMREX_SPACEDIM == 3)
	amrex::average_down_faces(m_leveldata[crse_lev+1]->hDVDFzF[Ark::FLUX],
			            	  m_leveldata[crse_lev]->hDVDFzF[Ark::FLUX],
                        	  refRatio(crse_lev),0);
#endif
#endif

#ifdef ARK_USE_MOMENTUM
	amrex::average_down_faces(m_leveldata[crse_lev+1]->fDVDFxF[Ark::FLUX],
			            	  m_leveldata[crse_lev]->fDVDFxF[Ark::FLUX],
                        	  refRatio(crse_lev),0);
	//
	amrex::average_down_faces(m_leveldata[crse_lev+1]->fDVDFyF[Ark::FLUX],
			            	  m_leveldata[crse_lev]->fDVDFyF[Ark::FLUX],
                        	  refRatio(crse_lev),0);
#if(AMREX_SPACEDIM == 3)
	amrex::average_down_faces(m_leveldata[crse_lev+1]->fDVDFzF[Ark::FLUX],
			            	  m_leveldata[crse_lev]->fDVDFzF[Ark::FLUX],
                        	  refRatio(crse_lev),0);
#endif
#endif	

#ifdef ARK_USE_ENERGY
	amrex::average_down_faces(m_leveldata[crse_lev+1]->gDVDFxF[Ark::FLUX],
			            	  m_leveldata[crse_lev]->gDVDFxF[Ark::FLUX],
                        	  refRatio(crse_lev),0);
	//
	amrex::average_down_faces(m_leveldata[crse_lev+1]->gDVDFyF[Ark::FLUX],
			            	  m_leveldata[crse_lev]->gDVDFyF[Ark::FLUX],
                        	  refRatio(crse_lev),0);
#if(AMREX_SPACEDIM == 3)
	amrex::average_down_faces(m_leveldata[crse_lev+1]->gDVDFzF[Ark::FLUX],
			            	  m_leveldata[crse_lev]->gDVDFzF[Ark::FLUX],
                        	  refRatio(crse_lev),0);
#endif
#endif
}
void DUGKSAmr::AverageDownAll()
{
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelAll(lev);
	}
}
void DUGKSAmr::AverageDownMacroVars()
{
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelMacroVars(lev);
	}
}
void DUGKSAmr::AverageDownPhase()
{	
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelPhase(lev);
    }
}
void DUGKSAmr::AverageDownDensity()
{	
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelDensity(lev);
    }
}
void DUGKSAmr::AverageDownVelocity()
{	
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelVelocity(lev);
    }
}
void DUGKSAmr::AverageDownPressure()
{	
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelPressure(lev);
    }
}
void DUGKSAmr::AverageDownTemperature()
{	
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelTemperature(lev);
    }
}
void DUGKSAmr::AverageDownDVDFFlux()
{
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelDVDFFlux(lev);
    }
}
void DUGKSAmr::AverageDownDVDFWave()
{
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelDVDFWave(lev);
    }
}
void DUGKSAmr::AverageDownDVDFBarPlus()
{
	for (int lev = finest_level-1; lev >= 0; --lev)
    {
		AverageDownLevelDVDFBarPlus(lev);
    }
}