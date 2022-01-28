#include <DUGKSAmr.H>

#if 0
---------------------------------------------------------------------
          CONVERSIONS FROM PHYSICAL TO MATHEMATICAL BCs


PHYS BC        |   Velocity      |    Temperature   |    Scalars     |
----------------------------------------------------------------------
INTERIOR       |                 |                  |                |
   or          |    INT_DIR      |     INT_DIR      |    INT_DIR     |
PERIODIC       |                 |                  |                |
----------------------------------------------------------------------
               |                 |                  |                |
INFLOW         |    EXT_DIR      |     EXT_DIR      |    EXT_DIR     |
               |                 |                  |                |
----------------------------------------------------------------------
               |                 |                  |                |
OUTFLOW        |    FOEXTRAP     |     FOEXTRAP     |    FOEXTRAP    |
               |                 |                  |                |
----------------------------------------------------------------------
NO SLIP WALL   |                 |                  |                |
  with         |    EXT_DIR      |  REFLECT_EVEN    |    HOEXTRAP    |
ADIABATIC TEMP |    U = V = 0    |    dT/dn = 0     |                |
----------------------------------------------------------------------
NO SLIP WALL   |                 |                  |                |
  with         |    EXT_DIR      |     EXT_DIR      |    HOEXTRAP    |
FIXED TEMP     |    U = V = 0    |                  |                |
----------------------------------------------------------------------
SLIP WALL      |                 |                  |                |
  with         |  Un=0 EXT_DIR   |  REFLECT_EVEN    |    HOEXTRAP    |
ADIABATIC TEMP |  Ut   HOEXTRAP  |    dT/dn = 0     |                |
----------------------------------------------------------------------
SLIP WALL      |                 |                  |                |
  with         |  Un=0 EXT_DIR   |     EXT_DIR      |    HOEXTRAP    |
FIXED TEMP     |  Ut   HOEXTRAP  |                  |                |
----------------------------------------------------------------------
               |                 |                  |                |
SYMMETRY       |  Un REFLECT_ODD |  REFLECT_EVEN    | REFLECT_EVEN   |
               |  Ut REFLECT_EVEN|                  |                |
----------------------------------------------------------------------


INT_DIR      : data taken from other grids or interpolated
EXT_DIR      : data specified on EDGE (FACE) of bndry for the linear solvers
               for application codes, the location of the EXT_DIR data will
               depend on how they do reconstruction, and may be edge or
               cell-centered.
HOEXTRAP     : higher order extrapolation to EDGE of bndry
HOEXTRAPCC   : linear extrapolation to CELL of bndry
FOEXTRAP     : first order extrapolation from last cell in interior
REFLECT_EVEN : F(-n) = F(n)   true reflection from interior cells
REFLECT_ODD  : F(-n) = -F(n)  true reflection from interior cells

namespace amrex {

namespace PhysBCType {
enum physicalBndryTypes : int {interior=0,inflow=1,outflow=2,symmetry=3,
                               slipwall=4,noslipwall=5};
}

namespace BCType {
enum mathematicalBndryTypes : int {
    bogus        = -666,
    reflect_odd  = -1,
    int_dir      =  0,
    reflect_even =  1,
    foextrap     =  2,
    ext_dir      =  3,
    hoextrap     =  4,
    hoextrapcc   =  5
};
}

}
#endif

void DUGKSAmr::InitBoundaryConditions()
{
	ReadBoundaryCondition("physical", m_bcrec_physical);
	ReadBoundaryCondition("phase", m_bcrec_phase);
	ReadBoundaryCondition("density", m_bcrec_density);
	ReadBoundaryCondition("velocity", m_bcrec_velocity, AMREX_SPACEDIM);
	ReadBoundaryCondition("pressure", m_bcrec_pressure);
	ReadBoundaryCondition("temperature", m_bcrec_temperature);
	ReadBoundaryCondition("force", m_bcrec_force, AMREX_SPACEDIM);

	ReadBoundaryCondition("DVDF_wave", m_bcrec_DVDF_wave, Ark::DVM_MOMENTUM);
	ReadBoundaryCondition("DVDF_barplus", m_bcrec_DVDF_barplus, Ark::DVM_MOMENTUM);

	ReadDirichletBoundaryCondition("dirichlet_phase", m_dirichlet_phase);
	ReadDirichletBoundaryCondition("dirichlet_density", m_dirichlet_density);
	ReadDirichletBoundaryCondition("dirichlet_velocity_x", m_dirichlet_velocity_x);
	ReadDirichletBoundaryCondition("dirichlet_velocity_y", m_dirichlet_velocity_y);
	ReadDirichletBoundaryCondition("dirichlet_velocity_z", m_dirichlet_velocity_z);
	ReadDirichletBoundaryCondition("dirichlet_pressure", m_dirichlet_pressure);
	ReadDirichletBoundaryCondition("dirichlet_temperature", m_dirichlet_temperature);

	SetPeriodicBoundaryCondition();

}
void DUGKSAmr::ReadBoundaryCondition(std::string const &var, 
									 Vector<BCRec> &m_bcrec, int const ncomp)
{
	{
		ParmParse pp("bc");

		std::string const var_lo = var + "_lo"; 
		std::string const var_hi = var + "_hi";

		Vector<int> bc_lo(AMREX_SPACEDIM, amrex::BCType::bogus);
    	Vector<int> bc_hi(AMREX_SPACEDIM, amrex::BCType::bogus);

    	pp.queryarr(var_lo.c_str(), bc_lo, 0, AMREX_SPACEDIM);
    	BL_ASSERT(bc_lo.size() == AMREX_SPACEDIM);

    	pp.queryarr(var_hi.c_str(), bc_hi, 0, AMREX_SPACEDIM);
    	BL_ASSERT(bc_hi.size() == AMREX_SPACEDIM);

    	m_bcrec.resize(ncomp);

    	for(int n = 0; n < ncomp; ++n)
    	{
    		for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    		{
    			m_bcrec[n].setLo(dir,bc_lo[dir]);
    			m_bcrec[n].setHi(dir,bc_hi[dir]);
    		}
    	}
	}
}
void DUGKSAmr::ReadDirichletBoundaryCondition(std::string const &var, BCValue &m_dirichlet)
{
	ParmParse pp("bc");

	std::string const var_lo = var + "_lo"; 
	std::string const var_hi = var + "_hi";

	Vector<Real> bc_value_lo(AMREX_SPACEDIM, 0.);
    Vector<Real> bc_value_hi(AMREX_SPACEDIM, 0.);

    pp.queryarr(var_lo.c_str(), bc_value_lo, 0, AMREX_SPACEDIM);
    pp.queryarr(var_hi.c_str(), bc_value_hi, 0, AMREX_SPACEDIM);

    for(int dir = 0; dir < AMREX_SPACEDIM; ++dir)
    {
    	m_dirichlet[dir]                = bc_value_lo[dir];
    	m_dirichlet[AMREX_SPACEDIM+dir] = bc_value_hi[dir];
    }

}
void DUGKSAmr::SetPeriodicityOverwrite(int const dir, Vector<BCRec> &m_bcrec_v)
{
    for(auto &bcrec : m_bcrec_v)
    {
    	bcrec.setLo(dir, amrex::BCType::int_dir);
    	bcrec.setHi(dir, amrex::BCType::int_dir);
    }
}
void DUGKSAmr::SetPeriodicityOverwrite(int const dir)
{
	SetPeriodicityOverwrite(dir, m_bcrec_physical);
	SetPeriodicityOverwrite(dir, m_bcrec_phase);
	SetPeriodicityOverwrite(dir, m_bcrec_density);
	SetPeriodicityOverwrite(dir, m_bcrec_velocity);
	SetPeriodicityOverwrite(dir, m_bcrec_pressure);
	SetPeriodicityOverwrite(dir, m_bcrec_temperature);
	SetPeriodicityOverwrite(dir, m_bcrec_force);
	SetPeriodicityOverwrite(dir, m_bcrec_DVDF_wave);
	SetPeriodicityOverwrite(dir, m_bcrec_DVDF_barplus);
}
void DUGKSAmr::SetPeriodicBoundaryCondition()
{
	if(Geom(0).isPeriodic(Ark::X))
	{
		SetPeriodicityOverwrite(Ark::X);
	}

	if(Geom(0).isPeriodic(Ark::Y))
	{
		SetPeriodicityOverwrite(Ark::Y);
	}

#if(AMREX_SPACEDIM == 3)
	if(Geom(0).isPeriodic(Ark::Z))
	{
		SetPeriodicityOverwrite(Ark::Z);
	}
#endif	
}
void DUGKSAmr::average_scalar_x (Array4<Real> const &dst, Array4<Real const> const &src, 
				  				 const Dim3 &lo, const Dim3 &hi, int const ix)
{
	LoopPS_YZ(lo,hi) 
	{
		dst(ix,j,k)  = 0.5*( src(ix - 1,j,k)  + src(ix,j,k) );
#ifdef AMREX_DEBUG
		if(dst(ix,j,k) != src(ix,j,k) || dst(ix,j,k) != src(ix-1,j,k)) {
			amrex::AllPrint().SetPrecision(12)
							<< dst(ix,j,k) <<"    "
							<< src(ix,j,k) <<"    "
						   	<< src(ix-1,j,k) <<"\n";
		}
#endif
	}
}
void DUGKSAmr::average_scalar_y (Array4<Real> const &dst, Array4<Real const> const &src, 
				  				 const Dim3 &lo, const Dim3 &hi, int const jy)
{
	LoopPS_XZ(lo,hi) 
	{
		dst(i,jy,k)  = 0.5*( src(i,jy - 1,k)  + src(i,jy,k) );
#ifdef AMREX_DEBUG
		if(dst(i,jy,k) != src(i,jy,k) || dst(i,jy,k) != src(i,jy-1,k)) 
		{
			amrex::AllPrint().SetPrecision(12)  
							<< dst(i,jy,k)<<"    "
							<< src(i,jy,k) <<"    "
							<< src(i,jy-1,k) <<"\n";
		}
#endif
	}
}
void DUGKSAmr::average_scalar_z (Array4<Real> const &dst, Array4<Real const> const &src, 
				  				 const Dim3 &lo, const Dim3 &hi, int const kz)
{
	LoopPS_XY(lo,hi) 
	{
		dst(i,j,kz)  = 0.5*( src(i,j,kz - 1)  + src(i,j,kz) );
	}
}

void DUGKSAmr::set_dirichlet_phase_xF (Array4<Real> const &phi, 
										const Dim3 &lo, const Dim3 &hi, 
							  			int const ix, const Orientation face)
{
	LoopPS_YZ(lo,hi)
	{
		phi(ix,j,k) = get_dirichlet_phase()[int(face)];
	}
}
void DUGKSAmr::set_dirichlet_phase_yF (Array4<Real> const &phi, 
										const Dim3 &lo, const Dim3 &hi, 
							  			int const jy, Orientation face)
{
	LoopPS_XZ(lo,hi)
	{
		phi(i,jy,k) = get_dirichlet_phase()[int(face)];
	}
}
void DUGKSAmr::set_dirichlet_phase_zF (Array4<Real> const &phi, 
										 const Dim3 &lo, const Dim3 &hi, 
							  			 int const kz, Orientation face)
{
	LoopPS_XY(lo,hi) 
	{
		phi(i,j,kz) = get_dirichlet_phase()[int(face)];
	}
}

void DUGKSAmr::set_dirichlet_density_xF (Array4<Real> const &rho, 
										 const Dim3 &lo, const Dim3 &hi, 
							  			 int const ix, const Orientation face)
{
	LoopPS_YZ(lo,hi)
	{
		rho(ix,j,k) = get_dirichlet_density()[int(face)];
	}
}
void DUGKSAmr::set_dirichlet_density_yF (Array4<Real> const &rho, 
										const Dim3 &lo, const Dim3 &hi, 
							  			int const jy, Orientation face)
{
	LoopPS_XZ(lo,hi)
	{
		rho(i,jy,k) = get_dirichlet_density()[int(face)];
	}
}
void DUGKSAmr::set_dirichlet_density_zF (Array4<Real> const &rho, 
										 const Dim3 &lo, const Dim3 &hi, 
							  			 int const kz, Orientation face)
{
	LoopPS_XY(lo,hi) 
	{
		rho(i,j,kz) = get_dirichlet_density()[int(face)];
	}
}

void DUGKSAmr::set_dirichlet_pressure_xF (Array4<Real> const &pres, 
										 const Dim3 &lo, const Dim3 &hi, 
							  			 int const ix, const Orientation face)
{
	LoopPS_YZ(lo,hi)
	{
		pres(ix,j,k) = get_dirichlet_pressure()[int(face)];
	}
}
void DUGKSAmr::set_dirichlet_pressure_yF (Array4<Real> const &pres, 
										 const Dim3 &lo, const Dim3 &hi, 
							  			 int const jy, Orientation face)
{
	LoopPS_XZ(lo,hi)
	{
		pres(i,jy,k) = get_dirichlet_pressure()[int(face)];		
	}
}
void DUGKSAmr::set_dirichlet_pressure_zF (Array4<Real> const &pres, 
										 const Dim3 &lo, const Dim3 &hi, 
							  			 int const kz, Orientation face)
{
	LoopPS_XY(lo,hi) 
	{
		pres(i,j,kz) = get_dirichlet_pressure()[int(face)];
	}
}

void DUGKSAmr::set_dirichlet_velocity_xF (Array4<Real> const &vel, 
										 const Dim3 &lo, const Dim3 &hi, 
							  			 int const ix, const Orientation face)
{
	LoopPS_YZ(lo,hi)
	{
		vel(ix,j,k,Ark::X) = get_dirichlet_velocity_x()[int(face)];
		vel(ix,j,k,Ark::Y) = get_dirichlet_velocity_y()[int(face)];
#if(AMREX_SPACEDIM == 3)
		vel(ix,j,k,Ark::Z) = get_dirichlet_velocity_z()[int(face)];
#endif
	}
}
void DUGKSAmr::set_dirichlet_velocity_yF (Array4<Real> const &vel, 
										 const Dim3 &lo, const Dim3 &hi, 
							  			 int const jy, Orientation face)
{
	LoopPS_XZ(lo,hi)
	{
		vel(i,jy,k,Ark::X) = get_dirichlet_velocity_x()[int(face)];
		vel(i,jy,k,Ark::Y) = get_dirichlet_velocity_y()[int(face)];
#if(AMREX_SPACEDIM == 3)
		vel(i,jy,k,Ark::Z) = get_dirichlet_velocity_z()[int(face)];
#endif		
	}
}
void DUGKSAmr::set_dirichlet_velocity_zF (Array4<Real> const &vel, 
										 const Dim3 &lo, const Dim3 &hi, 
							  			 int const kz, Orientation face)
{
	LoopPS_XY(lo,hi) 
	{
		vel(i,j,kz,Ark::X) = get_dirichlet_velocity_x()[int(face)];
		vel(i,j,kz,Ark::Y) = get_dirichlet_velocity_y()[int(face)];
#if(AMREX_SPACEDIM == 3)
		vel(i,j,kz,Ark::Z) = get_dirichlet_velocity_z()[int(face)];
#endif
	}
}