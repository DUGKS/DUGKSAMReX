#ifndef ARK_DUGKS_NAMESPACE_H
#define ARK_DUGKS_NAMESPACE_H

#include <AMReX_Box.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <DUGKSMacro.H>

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
	enum EnumDVDF {TAU,WAVE,EQ,SO,BP,BPx,BPy,BPxx,BPyy,BPxy,DVDFNUM};
#elif (AMREX_SPACEDIM == 3)
	enum EnumDVDF {TAU,WAVE,EQ,SO,BP,BPx,BPy,BPz,BPxx,BPyy,BPzz,BPxy,BPyz,BPzx,DVDFNUM};
#endif
	enum EnumFaceDVDF {TAUF,BHF,EQF,SOF,FLUX,FACEDVDFNUM};

	enum EnumMSQ {PHI,RHO,VEL,P,T,GradPHI,GradRHO,FORCE,MSVNUM};

	constexpr int DVM_PHASE    = 9;

	constexpr int DVM_MOMENTUM = 9;

	constexpr int DVM_ENERGY   = 9;

	constexpr int dir_x = static_cast<int>(amrex::Direction::x);

	constexpr int dir_y = static_cast<int>(amrex::Direction::y);

	#if(AMREX_SPACEDIM == 3)
	constexpr int dir_z = static_cast<int>(amrex::Direction::z);
	#endif

	constexpr double pi = 3.1415926535897932;

	constexpr double tiny = 1E-20;

	constexpr int digitsNum = 8;

	inline
	double Square(double const u,double const v, double const w)
	{
		return u*u + v*v + w*w;
	}
	inline
	double Square(double const u,double const v)
	{
		return u*u + v*v;
	}

	inline
	int sgn(double v) { return (v > 0) - (v < 0);}
}
namespace DVMCore
{
	// void Update_BarPlus(Vector<MultiFab> &xDVDF,const Real dTime);
	// void Update_Wave(Vector<MultiFab> &xDVDF,const Real dTime);
	// void Fill_Wave_Boundary
	// (
	// 	Vector<MultiFab> &xDVDF,const Real dTime,const amrex::Direction FaceID
	// );
}
namespace GradFunc
{
	//first PDs
	void Update_Grad_x
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	void Update_Grad_y
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	void Update_Grad_z
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	//Second PDs
	void Update_Grad_xx
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	void Update_Grad_yy
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	void Update_Grad_zz
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	void Update_Grad_xy
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	void Update_Grad_yz
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	void Update_Grad_zx
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	//third-order PDs
	void Update_Grad_xxx
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	void Update_Grad_yyy
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	//fourth-order PDs
	void Update_Grad_xxxx
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	void Update_Grad_yyyy
	(
		const Box &bx, int const ncomp, const Real *dx,
		Array4<Real const> const &src,Array4<Real> const &dst
	);
	//
	void gradient_operator_isotropy (const Box &bx, Real const* dl, 
						  			 Array4<Real const> const &src,
						  			 Array4<Real> const &dst);

	void laplacian_operator_isotropy (const Box &bx, Real const *dl, 
									  Array4<Real const> const &src, 
					  	   			  Array4<Real> const &dst);

	double laplacian_operator_isotropy (Real const *dl, Array4<Real const> const &var, 
					  	   				int const i,int const j,int const k);
}
//
namespace InitFlow
{
	void BasicTest
	(
		MultiFab &MsQ,MultiFab &MsQhxF,MultiFab &MsQhyF,const Real *dx
	);
	void Couette();
	void AC_InterfaceSmoothDeformation
	(
		MultiFab &MsQ,MultiFab &MsQhxF,MultiFab &MsQhyF,const Real *dx
	);
	void AC_InterfaceElongation
	(
		MultiFab &MsQ,MultiFab &MsQhxF,MultiFab &MsQhyF,const Real *dx
	);
	void AC_ZalesakDisk
	(
		MultiFab &MsQ,MultiFab &MsQhxF,MultiFab &MsQhyF,const Real *dx
	);	
	void AC_TransDrop
	(
		MultiFab &MsQ,MultiFab &MsQhxF,MultiFab &MsQhyF,const Real *dx
	);
	void AC_Drop
	(
		MultiFab &MsQ,MultiFab &MsQhxF,MultiFab &MsQhyF,const Real *dx
	);
	void Uniform(MultiFab &MsQ);
	void Shock(MultiFab &MsQ);
}

#endif