#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <DUGKSAmr.H>

int main(int argc, char* argv[])
{
	amrex::Initialize(argc,argv,true,MPI_COMM_WORLD);
	{

	BL_PROFILE("main()");

	// Issue an error if input file is not given 
	if(argc < 2) amrex::Abort("Input file must be given as command-line argument.");

	Real start_time = amrex::ParallelDescriptor::second();

	DUGKSAmr dugks_amr;

	dugks_amr.InitData();

	dugks_amr.Evolve();

	// Time spent in total
	Real end_time = amrex::ParallelDescriptor::second() - start_time;

	// Print timing results
	amrex::Print() << " Time spent in Evolve():      " << end_time << " " <<"seconds\n";
	Ark::PrintSplitLine();
	}
	amrex::Finalize();
}