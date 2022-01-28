#include <DUGKSAmr.H>

Ark::FlowParm::FlowParm()
{
	InitFlowParmBase();
	InitFlowParmPhaseField();
}
Ark::FlowParm::~FlowParm(){}

void Ark::FlowParm::InitFlowParmBase()
{
	{
	amrex::ParmParse ppCase("case");

	ppCase.get("case_index",case_index);
	ppCase.get("case_name",case_name);
	ppCase.get("advance_scheme",advance_scheme);
	ppCase.get("DmVn_name",DmVn_name);
	ppCase.get("flux_scheme",flux_scheme);
	ppCase.get("scheme_index",scheme_index);
	ppCase.get("solver_name",solver_name);
	}

	{
	amrex::ParmParse ppGeom("geometry");
	ppGeom.get("ch_num",ch_num);
	ppGeom.get("ch_length",ch_length);
	ppGeom.get("min_length",min_length);
	}

	{
	amrex::ParmParse ppBase("base");

	ppBase.get("Rho0",_Rho0);
	ppBase.get("U0",_U0);
	ppBase.get("P0",_P0);
	ppBase.get("T0",_T0);
	ppBase.get("R0",_R0);
	ppBase.get("Mu0",_Mu0);
	ppBase.get("Nu0",_Nu0);	
	ppBase.get("TauMomentum",_TauMomentum);
	}

	{
	//Ma, Re, and Kn 
	_Mach = _U0/sqrt(3*_R0*_T0);
	_Reynolds = _Rho0*_U0*ch_length/_Mu0;
	_Knudsen = 16.0*_Mu0/(5.0*_Rho0*_R0*_T0)*sqrt(_R0*_T0/(2.0*Ark::pi));
	}
}
 //----------------------------------------------------------
 //*********Initialization of Phase Field Parameters*********
 //----------------------------------------------------------
void Ark::FlowParm::InitFlowParmPhaseField()
{
	amrex::ParmParse ppPF("PF");

	ppPF.get("mobility_phi",mobility_phi);
	ppPF.get("tau_phi",tau_phi);
	ppPF.get("widthI",widthI);
	ppPF.get("sigma",sigma);

	ppPF.get("phiL",phiL);
	ppPF.get("phiV",phiV);
	ppPF.get("rhoL",rhoL);
	ppPF.get("rhoV",rhoV);
	ppPF.get("muL",muL);
	ppPF.get("muV",muV);
	ppPF.get("radiusPF",radiusPF);
	ppPF.getarr("centerPF",centerPF);
	ppPF.getarr("gravity",gravity);

	//Phase Field Calculations
	beta  = 12.0*sigma/widthI,
	kappa = 3.0*sigma*widthI/2;
	nuL = muL/rhoL;	
	nuV = muV/rhoV;
	aPhi = (rhoL - rhoV)/(phiL - phiV);
	bPhi = rhoV - aPhi*phiV;

	centerPF[Ark::X] *= ch_length;
	centerPF[Ark::Y] *= ch_length;
#if(AMREX_SPACEDIM == 3)
	centerPF[Ark::Z] *= ch_length;
#endif

	radiusPF    *= ch_length;

	_Cahn   = widthI/ch_length;
	_Peclet = _U0*ch_length/mobility_phi;
	_Atwood = (rhoL - rhoV)/(rhoL + rhoV);
	_Morton = fabs(gravity[Ark::Y])*(rhoL - rhoV)*pow(muL,4.0)/(rhoL*rhoL*pow(sigma,3));
	_Bond   = fabs(gravity[Ark::Y])*(rhoL - rhoV)*4.0*radiusPF*radiusPF/sigma;

	//rising bubble
	// _ReynoldsPF = sqrt(gravity[Ark::y]*rhoL*(rhoL - rhoV)*2*radius)*2*radius/muL,

	//Rayleigh-Taylor
	_ReynoldsPF = rhoL*ch_length*sqrt(fabs(gravity[Ark::Y]*ch_length))/muL;
	_Laplace = 0;	//_Sigma*_radius/_NuV/_NuV,
}