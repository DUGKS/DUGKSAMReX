#ifndef _ARK_FLOWPARM_H
#define _ARK_FLOWPARM_H

#include <string>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <DUGKSNameSpace.H>

// template<typename> class DUGKSCore;

namespace Ark{

class FlowParm
{
	// template <typename> friend class DUGKSCore;

public:
	//constructors
	FlowParm();

	FlowParm (FlowParm&& rhs) = default;
    FlowParm& operator= (FlowParm&& rhs) = default;

    FlowParm (const FlowParm& rhs) = default;
    FlowParm& operator= (const FlowParm& rhs) = default;

	virtual ~FlowParm();

	//case information
	const std::string& CaseName() const {return case_name;}
	const std::string& SolverName() const {return solver_name;}
	const std::string& FluxScheme() const { return flux_scheme; }
	const std::string& MomentumName() const {return momentum_name;}
	int CaseIndex() const {return case_index;}
	
	//characteristic length
	int ChNum() const {return ch_num;}
	double ChLength() const {return ch_length;}
	double MinL() const {return min_length;}

	//basic flow properties
	double Rho0() const {return _Rho0;};
	double U0() const {return _U0;}
	double P0() const {return _P0;}
	double T0() const {return _T0;}
	double R0() const {return _R0;}
	double Mu0() const {return _Mu0;}
	double Nu0() const {return _Nu0;}
	double TauMomentum() const {return _TauF;}  //the relaxation time for fDVDF
	double Ma() const {return _Mach;}
	double Re() const {return _Reynolds;}
	double Kn() const {return _Knudsen;}

	//control parameters
	int EndStep() const {return end_step;}
	double ResTol() const {return residual_tolerance;}
	double AdaptTol() const {return adapt_tolerance;}

	//phase field parameters
	double MPhi() const {return mobility_phi;}
	double TauPhase() const {return tau_phi;} //the relaxation time for hDVDF
	double WidthI() const {return widthI;}
	double Sigma() const {return sigma;}
	double Beta() const {return beta;}
	double Kappa() const {return kappa;}
	double PhiL() const {return phiL;}
	double PhiV() const {return phiV;}
	double RhoL() const {return rhoL;}
	double RhoV() const {return rhoV;}
	double MuL() const {return muL;}
	double MuV() const {return muV;}
	double NuL() const {return nuL;}
	double NuV() const {return nuV;}
	double APhi() const {return aPhi;}
	double BPhi() const {return bPhi;}
	double RadiusPF() const {return radiusPF;}
	double Gravity(int i) const {return gravity[i];}
	double CenterPF(int i) const {return centerPF[i];}

	double Cn() const {return _Cahn;}
	double Pe() const {return _Peclet;}
	double At() const {return _Atwood;}
	double Mo() const {return _Morton;}
	double Bo() const {return _Bond;}
	double La() const {return _Laplace;}
	double RePF() const {return _ReynoldsPF;}

	//other member funcs defined by users
	
private:
	std::string case_name;
	std::string solver_name;
	std::string flux_scheme;
	std::string momentum_name;
	int case_index = 0;

	int    ch_num = 0;	//characteristic cell number per side
	double ch_length = 0; //characteristic length
	double min_length = 0; //minimum characteristic length

	//physical parameters
	double _Rho0 	 = 0;
	double _U0 		 = 0;
	double _P0 		 = 0;
	double _T0 		 = 0;
	double _R0 		 = 0;
	double _Mu0      = 0;
	double _Nu0      = 0;
	double _TauF     = 0;
	double _Mach     = 0;
	double _Reynolds = 0;
	double _Knudsen  = 0;

	//control parameters
	int end_step = 0;
	double residual_tolerance = 0;
	double adapt_tolerance = 0;

	//phase field parameters
	double mobility_phi    = 0; //mobility coefficient
	double tau_phi    = 0; //relaxation time of index equation
	double widthI      = 0; //interface width
	double sigma   = 0; //surface tension coefficient
	double beta    = 0;
	double kappa   = 0;

	double phiL = 0, phiV = 0;
	double rhoL = 0, rhoV = 0;
	double muL  = 0, muV  = 0;
	double nuL  = 0, nuV  = 0;
	double aPhi = 0, bPhi = 0;
	double radiusPF = 0;
	amrex::Vector<double> gravity;
	amrex::Vector<double> centerPF;

	double _Cahn = 0;
	double _Peclet  = 0;
	double _Atwood  = 0;
	double _Morton  = 0;
	double _Bond    = 0;
	double _Laplace = 0;
	double _ReynoldsPF = 0;

	//other parameters defined by users

private:
	//initial basic parameters
	void InitFlowParmBase();

	//initial phase field parameters
	void InitFlowParmPhaseFiled();

	//other functions for the initialization of parameters provided by users
};

}

#endif