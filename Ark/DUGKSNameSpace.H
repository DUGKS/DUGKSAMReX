#ifndef ARK_DUGKS_NAMESPACE_H
#define ARK_DUGKS_NAMESPACE_H

#include <AMReX_Box.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <DUGKSMacro.H>
#include <DUGKSDmVn.H>

using amrex::Box;
using amrex::Real;
using amrex::Array;
using amrex::Array4;
using amrex::IntVect;
using amrex::Vector;
using amrex::MultiFab;



namespace Ark
{
#if(AMREX_SPACEDIM == 2)
	enum EnumDVDF {TAU,WAVE,EQ,SO,BP,BPx,BPy,BPxx,BPyy,BPxy,BPxxx,BPxxy,BPxyy,BPyyy,BPxxxx,BPyyyy,BPxxyy,DVDFNUM};
	// enum EnumDVDF {TAU,WAVE,EQ,SO,BP,BPx,BPy,BPxx,BPyy,BPxy,DVDFNUM};

	using DmVnPhase    = D2M9;

	using DmVnMomentum = D2M9;

	using DmVnEnergy   = D2M9;

#elif (AMREX_SPACEDIM == 3)
	enum EnumDVDF {TAU,WAVE,EQ,SO,BP,BPx,BPy,BPz,BPxx,BPyy,BPzz,BPxy,BPyz,BPzx,DVDFNUM};

	using DmVnPhase    = D3V7;

	using DmVnMomentum = D3V15;

	using DmVnEnergy   = D3V15;
#endif

	enum EnumFaceDVDF {TAUF,FLUX,EQF,SOF,FACEDVDFNUM};

	enum EnumMSQ {PHI,RHO,VEL,P,T,GradPHI,GradRHO,FORCE,CHEMIP,MSVNUM};

	constexpr unsigned DVM_PHASE    = DmVnPhase::Q;

	constexpr unsigned DVM_MOMENTUM = DmVnMomentum::Q;

	constexpr unsigned DVM_ENERGY   = DmVnEnergy::Q;

	constexpr int X = static_cast<int>(amrex::Direction::x);

	constexpr int Y = static_cast<int>(amrex::Direction::y);

	#if(AMREX_SPACEDIM == 3)
	constexpr int Z = static_cast<int>(amrex::Direction::z);
	#endif

	constexpr double pi = 3.1415926535897932;
	constexpr double PI = pi;

	constexpr double tiny = 1E-15;

	constexpr unsigned digitsNum = 8;

	inline
	double Square(double const u, double const v = 0, double const w = 0)
	{
		return u*u + v*v + w*w;
	}

	inline
	double Hypot(double const u, double const v, double const w)
	{
		return std::sqrt(u*u + v*v + w*w);
	}

	template <typename T>
	bool approximatelyEqual(T A, T B)
	{
		return std::fabs(A-B) < std::numeric_limits<T>::epsilon();
	}

	template <typename T>
	bool tinyEqual(const T A, const T B, const T diff = T(tiny))
	{
		return std::fabs(A-B) < diff;
	}

	template <typename T>
	int sign (const T val) { return (val > 0) - (val < 0); }

	void PrintSplitLine(const char c = ' ');
}

namespace GradFunc
{
	//first PDs
	void Update_Grad_x (const Box &bx, int const ncomp, const Real *dx,
						Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_y (const Box &bx, int const ncomp, const Real *dx,
						Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_z (const Box &bx, int const ncomp, const Real *dx,
						Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_x_5th (const Box &bx, int const ncomp, const Real *dx,
							Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_y_5th (const Box &bx, int const ncomp, const Real *dx,
							Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_z_5th (const Box &bx, int const ncomp, const Real *dx,
							Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_VanLeer_x (const Box &bx, int const ncomp, const Real *dx,
								Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_VanLeer_y (const Box &bx, int const ncomp, const Real *dx,
								Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_VanLeer_z (const Box &bx, int const ncomp, const Real *dx,
								Array4<Real const> const &src,Array4<Real> const &dst);

	//Second PDs
	void Update_Grad_xx (const Box &bx, int const ncomp, const Real *dx,
						 Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_yy (const Box &bx, int const ncomp, const Real *dx,
						 Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_zz (const Box &bx, int const ncomp, const Real *dx,
						 Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_xx_5th	(const Box &bx, int const ncomp, const Real *dx,
							 Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_yy_5th	(const Box &bx, int const ncomp, const Real *dx,
							 Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_zz_5th	(const Box &bx, int const ncomp, const Real *dx,
							 Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_xy	(const Box &bx, int const ncomp, const Real *dx,
						 Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_yz	(const Box &bx, int const ncomp, const Real *dx,
						 Array4<Real const> const &src,Array4<Real> const &dst);

	void Update_Grad_zx	(const Box &bx, int const ncomp, const Real *dx,
						 Array4<Real const> const &src,Array4<Real> const &dst);

	//third-order PDs
	void Update_Grad_xxx (const Box &bx, int const ncomp, const Real *dx,
						  Array4<Real const> const &src, Array4<Real> const &dst);

	void Update_Grad_xxy (const Box &bx, int const ncomp, const Real *dx,
						  Array4<Real const> const &src, Array4<Real> const &dst);

	void Update_Grad_xyy (const Box &bx, int const ncomp, const Real *dx,
						  Array4<Real const> const &src, Array4<Real> const &dst);

	void Update_Grad_yyy (const Box &bx, int const ncomp, const Real *dx,
						  Array4<Real const> const &src, Array4<Real> const &dst);
	//fourth-order PDs
	void Update_Grad_xxxx (const Box &bx, int const ncomp, const Real *dx,
						   Array4<Real const> const &src, Array4<Real> const &dst);

	void Update_Grad_xxyy (const Box &bx, int const ncomp, const Real *dx,
						   Array4<Real const> const &src, Array4<Real> const &dst);

	void Update_Grad_yyyy (const Box &bx, int const ncomp, const Real *dx,
						   Array4<Real const> const &src, Array4<Real> const &dst);

	//

	void gradient_operator_central (const Box &bx, Real const* dl, Array4<Real const> const &src,
						  			Array4<Real> const &dst);
	void gradient_operator_isotropy (const Box &bx, Real const* dl, Array4<Real const> const &src,
						  			 Array4<Real> const &dst);

	void gradient_operator_isotropy_biased (const Box &bx, Real const* dl, Array4<Real const> const &src,
						  			 		Array4<Real> const &dst);

	void gradient_operator_isotropy_blend (const Box &bx, Real const* dl, Array4<Real const> const &src,
						  			 	   Array4<Real> const &dst);

	void gradient_operator_5th (const Box &bx, Real const* dl, Array4<Real const> const &src,
								Array4<Real> const &dst);

	void laplacian_operator_isotropy (const Box &bx, Real const *dl, Array4<Real const> const &src,
					  	   			  Array4<Real> const &dst);

	void laplacian_operator_5th (const Box &bx, Real const *dl, Array4<Real const> const &src,
					  	   		 Array4<Real> const &dst);

	void laplacian_normal_operator_isotropy (const Box &bx, Real const* dl,
											 Array4<Real const> const &src,
								 			 Array4<Real const> const &gradsrc,
								 			 Array4<Real> const &dst);

	void laplacian_normal_operator_5th (const Box &bx, Real const* dl,
										Array4<Real const> const &src,
								 		Array4<Real const> const &gradsrc,
								 		Array4<Real> const &dst);

	void laplacian_normal_chemip_isotropy (const Box &bx, Real const* dl,
											 Array4<Real const> const &src,
								 			 Array4<Real const> const &gradsrc,
								 			 Array4<Real> const &dst);

	double laplacian_operator_isotropy (Real const *dl, Array4<Real const> const &var,
					  	   				int const i,int const j,int const k);

	// double laplacian_operator_isotropy (Real const *dl, Array4<Real const> const &var,
	// 				  	   				int const i,int const j,int const k);

	double double_interpolation (const double upstream1, const double upstream0,
                            	 const double downstream0, const double downstream1,
								 int const token);

	double triple_interpolation (const double upstream2, const double upstream1,
                                 const double upstream0, const double downstream0,
                                 const double downstream1, const double downstream2,
                                 const int token);

	double WENO_Z (const double upstream2, const double upstream1, const double upstream0,
				   const double downstream0, const double downstream1);

	double WENO_C (const double upstream2, const double upstream1,
                   const double upstream0, const double downstream0,
                   const double downstream1, const double downstream2);
}

namespace EoS {
// The chemical potential here is the bulk chemical potential, i.e. \mu_0 = d{\psi}\d{\rho}

//  Van Der Waals Equation Li2021PRE [103, 013304 (2021)]

	double PressureVdWIdeal (const double Density, const double Temp);

	double PressureVdW (const double Density, const double Temp);

	double FreeEnergyDensityVdW (const double Density, const double Temp);

	double ChemicalPotentialVdW (const double Density, const double Temp);

// Redlich-Kwong-Soave Wen2017PRE [95, 063305 (2017)]

	double AlphaRKS (const double Temp);

	double PressureRKS (const double Density, const double Temp);

	double FreeEnergyDensityRKS (const double Density, const double Temp);

	double ChemicalPotentialRKS (const double Density, const double Temp);

// Peng-Robinson Yuan2006PoF[https://doi.org/10.1063/1.2187070]

	double AlphaPR (const double Temp);

	double PressurePR (const double Density, const double Temp);

	double FreeEnergyDensityPR (const double Density, const double Temp);

	double ChemicalPotentialPR (const double Density, const double Temp);

// Carnahan-Starling  Yuan2006PoF[https://doi.org/10.1063/1.2187070]

	double PressureCS (const double Density, const double Temp);

	double FreeEnergyDensityCS (const double Density, const double Temp);

	double ChemicalPotentialCS (const double Density, const double Temp);

// Double Well Potential
	double ChemicalPotentialDW (const double order, const double orderA, const double orderB);

// Pseudo-Potential
	double PseudoPotential(const double Density, const double Tr);
}

namespace Limiter
{
	double VanAlbada2(const double vF, const double vB);

	double VanLeer(const double vF, const double vB);
}

namespace Quadrature
{
	double NInt (int const Num, double const *f, double const *c);
}

#endif