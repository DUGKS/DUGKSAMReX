#include <DUGKSAmr.H>

#if(AMREX_SPACEDIM == 2)
void GradFunc::Update_Grad_x
(
	const Box &bx, int const ncomp, Real const *dx,
	Array4<Real const> const &src, Array4<Real> const &dst
)
{
	double const dlx = dx[Ark::X];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i+1,j,  k,n) - src(i-1,j,  k,n))*4.0
		+	(src(i+1,j+1,k,n) - src(i-1,j+1,k,n))
		+	(src(i+1,j-1,k,n) - src(i-1,j-1,k,n))
		)/(12.0*dlx);
	}
	);
}
void GradFunc::Update_Grad_y
(
	const Box &bx, int const ncomp, Real const *dx,
	Array4<Real const> const &src, Array4<Real> const &dst
)
{
	double const dly = dx[Ark::Y];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i,  j+1,k,n) - src(i,  j-1,k,n))*4.0
		+	(src(i+1,j+1,k,n) - src(i+1,j-1,k,n))
		+	(src(i-1,j+1,k,n) - src(i-1,j-1,k,n))
		)/(12.0*dly);
	}
	);
}
void GradFunc::Update_Grad_xx
(
	const Box &bx, int const ncomp, Real const *dx,
	Array4<Real const> const &src, Array4<Real> const &dst
)
{
	double const dlx = dx[Ark::X];

	ParallelFor(bx, ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i-1,j,  k,n) - 2.0*src(i,j,  k,n) + src(i+1,j,  k,n))*10.0
		+	(src(i-1,j+1,k,n) - 2.0*src(i,j+1,k,n) + src(i+1,j+1,k,n))
		+	(src(i-1,j-1,k,n) - 2.0*src(i,j-1,k,n) + src(i+1,j-1,k,n))
		)/(12.0*dlx*dlx);
	}
	);
}
void GradFunc::Update_Grad_yy
(
	const Box &bx, int const ncomp,Real const *dx,
	Array4<Real const> const &src, Array4<Real> const &dst
)
{
	double const dly = dx[Ark::Y];

	ParallelFor(bx,ncomp,
	[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
	{
		dst(i,j,k,n) = 
		(
			(src(i,  j-1,k,n) - 2.0*src(i,  j,k,n) + src(i,  j+1,k,n))*10.0
		+	(src(i+1,j-1,k,n) - 2.0*src(i+1,j,k,n) + src(i+1,j+1,k,n))
		+	(src(i-1,j-1,k,n) - 2.0*src(i-1,j,k,n) + src(i-1,j+1,k,n))
		)/(12.0*dly*dly);
	}
	);
}
void GradFunc::Update_Grad_xy
(
	const Box &bx, int const ncomp,Real const *dx,
	Array4<Real const> const &src, Array4<Real> const &dst
)
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
	}
	);
}
void GradFunc::gradient_operator_isotropy (const Box &bx, Real const* dl, 
								 		   Array4<Real const> const &src,
								 		   Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	AMREX_ASSERT(dlx == dly);
	AMREX_ASSERT(dst.nComp() == AMREX_SPACEDIM);

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::X) = 
		(
			(src(i+1,j  ,k) - src(i-1,j  ,k))*4.0
		+	(src(i+1,j-1,k) - src(i-1,j-1,k))
		+	(src(i+1,j+1,k) - src(i-1,j+1,k))
		)/(12.0*dlx);
	});
	//
	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k,Ark::Y) =
		(
			(src(i,  j+1,k) - src(i,  j-1,k))*4.0
		+	(src(i-1,j+1,k) - src(i-1,j-1,k))
		+	(src(i+1,j+1,k) - src(i+1,j-1,k))
		)/(12.0*dly);
	});
}
double GradFunc::laplacian_operator_isotropy (Real const *dl, Array4<Real const> const &var, 
					  	   					  int const i,int const j,int const k)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];
	
	AMREX_ASSERT(dlx == dly);
//
	return
	(	
		(var(i+1,j,  k) + var(i-1,j,  k) + var(i,j+1,  k) + var(i,j-1,  k))*4.0
	+	(var(i+1,j+1,k) + var(i-1,j+1,k) + var(i-1,j-1,k) + var(i+1,j-1,k))
	-	var(i,j,k)*20
	)/6.0*dlx*dly;

}					  		   
void GradFunc::laplacian_operator_isotropy (const Box &bx, Real const *dl, 
										    Array4<Real const> const &src, 
					  	   					Array4<Real> const &dst)
{
	const Real dlx = dl[Ark::X];
	const Real dly = dl[Ark::Y];

	AMREX_ASSERT(dlx == dly);

	ParallelFor(bx,
	[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
	{
		dst(i,j,k) = 
		(
			(	
				src(i+1,j,  k) + src(i-1,j,  k)
			+	src(i,  j+1,k) + src(i,  j-1,k)
			)*4.0
		+	
			(	
				src(i+1,j+1,k) + src(i-1,j+1,k)
			+	src(i-1,j-1,k) + src(i+1,j-1,k)
			)
		-	src(i,j,k)*20
		)/(6.0*dlx*dly);
	}
	);
}
#endif