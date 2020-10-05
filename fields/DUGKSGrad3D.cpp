#include <DUGKSAmr.H>

#if(AMREX_SPACEDIM == 3)
void GradFunc::Update_Grad_x (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::dir_x];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i+1,j,k,n) - src(i-1,j,k,n))/(2.0*dlx);
	});
}
void GradFunc::Update_Grad_y (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::dir_y];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i,j+1,k,n) - src(i,j-1,k,n))/(2.0*dly);
	});
}
void GradFunc::Update_Grad_z (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::dir_z];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i,j,k+1,n) - src(i,j,k-1,n))/(2.0*dlz);
	});
}
void GradFunc::Update_Grad_xx (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::dir_x];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i-1,j,k,n) - 2.0*src(i,j,k,n) + src(i+1,j,k,n)))/(dlx*dlx);
	});
}
void GradFunc::Update_Grad_yy (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::dir_y];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i,j-1,k,n) - 2.0*src(i,j,k,n) + src(i,j+1,k,n)))/(dly*dly);
	});
}
void GradFunc::Update_Grad_zz (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::dir_z];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i,j,k-1,n) - 2.0*src(i,j,k,n) + src(i,j,k+1,n)))/(dlz*dlz);
	});
}
void GradFunc::Update_Grad_xy (const Box &bx, int const ncomp,Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dlx = dx[Ark::dir_x];
	const Real dly = dx[Ark::dir_y];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			src(i+1,j+1,k,n) + src(i-1,j-1,k,n)
		-	src(i+1,j-1,k,n) - src(i-1,j+1,k,n)
		)/(4.0*dlx*dly);
	});
}
void GradFunc::Update_Grad_yz (const Box &bx, int const ncomp,Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dly = dx[Ark::dir_y];
	const Real dlz = dx[Ark::dir_z];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			src(i,j+1,k+1,n) + src(i,j-1,k-1,n)
		-	src(i,j+1,k-1,n) - src(i,j-1,k+1,n)
		)/(4.0*dly*dlz);
	});
}
void GradFunc::Update_Grad_zx (const Box &bx, int const ncomp,Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dlz = dx[Ark::dir_z];
	const Real dlx = dx[Ark::dir_x];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			src(i+1,j,k+1,n) + src(i-1,j,k-1,n)
		-	src(i-1,j,k+1,n) - src(i+1,j,k-1,n)
		)/(4.0*dlz*dlx);
	});
}
void GradFunc::gradient_operator_isotropy (const Box &bx, Real const* dl, 
								 		   Array4<Real const> const &src,
								 		   Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::dir_x];
	const Real dly = dl[Ark::dir_y];
	const Real dlz = dl[Ark::dir_z];

	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dlx == dlz);

	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::dir_x) = 0;
	});
	//
	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::dir_y) =	0;
	});
	//
	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::dir_z) =	0;
	});
}
double GradFunc::laplacian_operator_isotropy (Real const *dl, Array4<Real const> const &var, 
					  	   					  int const i,int const j,int const k)
{
	const Real dlx = dl[Ark::dir_x];
	const Real dly = dl[Ark::dir_y];
	const Real dlz = dl[Ark::dir_z];
	
	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dly == dlz);
//
	return 0;
}					  		   
void GradFunc::laplacian_operator_isotropy (const Box &bx, Real const *dl, 
										    Array4<Real const> const &src, 
					  	   					Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::dir_x];
	const Real dly = dl[Ark::dir_y];
	const Real dlz = dl[Ark::dir_z];
	
	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dly == dlz);

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k) = 0;
	}
	);
}
#endif