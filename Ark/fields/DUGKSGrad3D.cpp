#include <DUGKSAmr.H>

#if(AMREX_SPACEDIM == 3)
void GradFunc::Update_Grad_x (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i+1,j,k,n) - src(i-1,j,k,n))/(2.0*dlx);
	});
}
void GradFunc::Update_Grad_y (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i,j+1,k,n) - src(i,j-1,k,n))/(2.0*dly);
	});
}
void GradFunc::Update_Grad_z (const Box &bx, int const ncomp, Real const *dx,
							  Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::Z];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = (src(i,j,k+1,n) - src(i,j,k-1,n))/(2.0*dlz);
	});
}
void GradFunc::Update_Grad_xx (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlx = dx[Ark::X];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i-1,j,k,n) - 2.0*src(i,j,k,n) + src(i+1,j,k,n)))/(dlx*dlx);
	});
}
void GradFunc::Update_Grad_yy (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dly = dx[Ark::Y];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i,j-1,k,n) - 2.0*src(i,j,k,n) + src(i,j+1,k,n)))/(dly*dly);
	});
}
void GradFunc::Update_Grad_zz (const Box &bx, int const ncomp, Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	double const dlz = dx[Ark::Z];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = ((src(i,j,k-1,n) - 2.0*src(i,j,k,n) + src(i,j,k+1,n)))/(dlz*dlz);
	});
}
void GradFunc::Update_Grad_xy (const Box &bx, int const ncomp,Real const *dx,
							   Array4<Real const> const &src, Array4<Real> const &dst)
{
	const Real dlx = dx[Ark::X];
	const Real dly = dx[Ark::Y];

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
	const Real dly = dx[Ark::Y];
	const Real dlz = dx[Ark::Z];

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
	const Real dlz = dx[Ark::Z];
	const Real dlx = dx[Ark::X];

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
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];

	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dlx == dlz);

	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::X) = 
		(
			(src(i+1,j,k) - src(i-1,j,k))*16

		+	(
				src(i+1,j+1,k) + src(i+1,j-1,k) + src(i+1,j,k+1) + src(i+1,j,k-1) 
			-	src(i-1,j+1,k) - src(i-1,j-1,k) - src(i-1,j,k+1) - src(i-1,j,k-1) 
			)*4

		+	(
			  	src(i+1,j+1,k+1) + src(i+1,j-1,k+1) + src(i+1,j+1,k-1) + src(i+1,j-1,k-1)
			- 	src(i-1,j+1,k+1) - src(i-1,j-1,k+1) - src(i-1,j+1,k-1) - src(i-1,j-1,k-1)
			)
		)/(72*dlx);
	});
	//
	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Y) =	
		(
			(src(i,j+1,k) - src(i,j-1,k))*16

		+	(
				src(i+1,j+1,k) + src(i-1,j+1,k) + src(i,j+1,k+1) + src(i,j+1,k-1) 
			-	src(i+1,j-1,k) - src(i-1,j-1,k) - src(i,j-1,k+1) - src(i,j-1,k-1) 
			)*4

		+	(
			  	src(i+1,j+1,k+1) + src(i-1,j+1,k+1) + src(i+1,j+1,k-1) + src(i-1,j+1,k-1)
			- 	src(i+1,j-1,k+1) - src(i-1,j-1,k+1) - src(i+1,j-1,k-1) - src(i-1,j-1,k-1)
			)
		)/(72*dly);
	});
	//
	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Z) =
		(
			(src(i,j,k+1) - src(i,j,k-1))*16

		+	(
				src(i,j+1,k+1) + src(i,j-1,k+1) + src(i+1,j,k+1) + src(i-1,j,k+1) 
			-	src(i,j+1,k-1) - src(i,j-1,k-1) - src(i+1,j,k-1) - src(i-1,j,k-1) 
			)*4

		+	(
			  	src(i+1,j+1,k+1) + src(i+1,j-1,k+1) + src(i-1,j+1,k+1) + src(i-1,j-1,k+1)
			- 	src(i+1,j+1,k-1) - src(i+1,j-1,k-1) - src(i-1,j+1,k-1) - src(i-1,j-1,k-1)
			)
		)/(72*dlz);
	});
}
double GradFunc::laplacian_operator_isotropy (Real const *dl, Array4<Real const> const &var, 
					  	   					  int const i,int const j,int const k)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];
	
	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dly == dlz);

	const Real dln = dlx;
//	
	Real tmp =
			(
				(
				  	var(i+1,j  ,k)   + var(i-1,j  ,k)
				+ 	var(i,  j+1,k)   + var(i,  j-1,k) 
				+ 	var(i,  j  ,k+1) + var(i,  j  ,k-1)
				)*20.
	//
			+	(
				  	var(i+1,j+1,k) + var(i-1,j+1,k) + var(i+1,j-1,k) + var(i-1,j-1,k)
				+ 	var(i+1,j,k+1) + var(i-1,j,k+1) + var(i+1,j,k-1) + var(i-1,j,k-1)
				+ 	var(i,j+1,k+1) + var(i,j-1,k+1) + var(i,j+1,k-1) + var(i,j-1,k-1)
				)*6.
	//
			+	(
					var(i+1,j+1,k+1) + var(i+1,j-1,k+1) + var(i+1,j+1,k-1) + var(i+1,j-1,k-1)
				+	var(i-1,j+1,k+1) + var(i-1,j-1,k+1) + var(i-1,j+1,k-1) + var(i-1,j-1,k-1)
				)
			-   var(i,j,k)*200
			)/(48*dln*dln);

	return tmp;
}				
//doi:10.1016/j.jcp.2004.05.005	  		   
void GradFunc::laplacian_operator_isotropy (const Box &bx, Real const *dl, 
										    Array4<Real const> const &src, 
					  	   					Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	const Real dlz = dl[Ark::Z];
	
	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dly == dlz);

	const Real dln = dlx;

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k) = 
		(
			(
			  	src(i+1,j  ,k)   + src(i-1,j  ,k)
			+ 	src(i,  j+1,k)   + src(i,  j-1,k) 
			+ 	src(i,  j  ,k+1) + src(i,  j  ,k-1)
			)*20.
//
		+	(
			  	src(i+1,j+1,k) + src(i-1,j+1,k) + src(i+1,j-1,k) + src(i-1,j-1,k)
			+ 	src(i+1,j,k+1) + src(i-1,j,k+1) + src(i+1,j,k-1) + src(i-1,j,k-1)
			+ 	src(i,j+1,k+1) + src(i,j-1,k+1) + src(i,j+1,k-1) + src(i,j-1,k-1)
			)*6.
//
		+	(
				src(i+1,j+1,k+1) + src(i+1,j-1,k+1) + src(i+1,j+1,k-1) + src(i+1,j-1,k-1)
			+	src(i-1,j+1,k+1) + src(i-1,j-1,k+1) + src(i-1,j+1,k-1) + src(i-1,j-1,k-1)
			)
		-   src(i,j,k)*200
		)/(48*dln*dln);
	});
}
#endif

double GradFunc::central_interpolation (const double upstream1, const double upstream0,
                            		    const double downstream0, const double downstream1)
{
    // return ( (upstream2 + downstream2)*3 
    //        - (upstream1 + downstream1)*25.0
    //        + (upstream0 + downstream0)*150.0)/256.0;

    return ( (upstream0 + downstream0)*9.0 - (upstream1 + downstream1))/16.0;
}