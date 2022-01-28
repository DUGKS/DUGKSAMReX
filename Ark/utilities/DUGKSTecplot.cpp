#include <DUGKSAmr.H>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace::amrex;

using std::string;
using std::ofstream;
using std::ostringstream;
using std::ios;
using std::setprecision;

int const Out_precision = 12;

void TecplotGridHeader (ofstream &OutFile,const string &caseName,const string &varList,
					    int const NLx,int const NLy)
{
	OutFile << "FILETYPE = GRID\n"
			<< "TITLE = " << caseName << "\n"
			<< "VARIABLES = " << varList << "\n"
			<< "ZONE T = Grid\n"
			<< "NODES = " << (NLx+1)*(NLy+1)
			<< ", ELEMENTS = " << NLx*NLy
			<< ", DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n";
}

void TecplotSolutionHeader (ofstream &OutFile,const string &caseName,
						  const string &varList,const string &varNum,
						  int const step,int const NLx,int const NLy)
{
	OutFile << "FILETYPE = SOLUTION\n"
			<< "TITLE = Data\n"
			<< "VARIABLES = "<<varList<<"\n"
			<< "ZONE T=step"<<step<<"\n"
			<< "StrandID=1, SOLUTIONTIME="<<step<<"\n"
			<< "NODES="<<(NLx+1)*(NLy+1)
			<< ",ELEMENTS = "<< NLx*NLy
			<< ",DATAPACKING=BLOCK,ZONETYPE=FEQUADRILATERAL\n"
			<< "VARLOCATION=("<<varNum<<"=CELLCENTERED)\n";
}

// void oss_XXX(ostringstream& oss_r,const string &folder,const string &suffix,double const &t)
// {
// 	oss_r << "../FlowField/"<<folder<<"/Time" <<t<<"."<<suffix;
// }
void FileOpen(ofstream &OutFile_XXX, ostringstream &oss_XXX,string const &s)
{
	OutFile_XXX.open(oss_XXX.str().c_str());
	if(!OutFile_XXX)
	{
		amrex::Abort(("OutFile_" + s + " open failed"));
		return;
	}
	OutFile_XXX << setiosflags(ios::scientific) << setprecision(Out_precision);
}
void FileOpenAppend(ofstream &OutFile_XXX, ostringstream &oss_XXX, string const &errorinfo)
{
	OutFile_XXX.open(oss_XXX.str().c_str(),ios::app);
	if(!OutFile_XXX)
	{
		amrex::Abort(("OutFile_" + errorinfo + " open failed"));
		return;
	}
	OutFile_XXX << setiosflags(ios::scientific) << setprecision(Out_precision);
}
// void oss_XXXSet(ostringstream &oss_XXX,const string &fileName,const string &caseName)
// {
// 	oss_XXX<<"./"<<fileName<<_ARK_FILENAME
// 		   <<"_NL"<<DB.ChNumL()<<"_"<<caseName<<".dat";
// }

void writeSolutionData(ofstream& OutFile, const MultiFab &mf)
{
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(mf, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		//- region
		const amrex::Box& bx = mfi.validbox();
		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		//source
		Array4<Real const> const &mfArray = mf.const_array(mfi);

		LoopPS(lo,hi) {
			OutFile << mfArray(i,j,k) << std::endl;
		}
	}
}

void DUGKSAmr::WriteTecplotFlowField(int const lev) 
{
	int const NLx = Geom(lev).Domain().length(Ark::X);
	int const NLy = Geom(lev).Domain().length(Ark::Y);

	const MultiFab &PhiFab = Data(lev).msv[Ark::PHI];

	ostringstream oss_FlowField;

	oss_FlowField << "step" << nstep
            	  <<"_Pe"   << Pe()
            	  <<"_M"    << MPhi()
            	  <<"_wI"   << WidthI()
            	  <<"_NL"   << ChNum()
            	  <<"_CFL"  << CFL()
            	  <<"_Flux" << FluxScheme()
            	  <<"_"<<CaseName()<<".dat";

	ofstream OutFile_FlowField;
	FileOpen(OutFile_FlowField, oss_FlowField, "FlowField");

	// string varList = "Phi,Phix,Phiy,U,V,P,Rhox,Rhoy,Fx,Fy";
	string varList = "Phi";
	string varNum  = "[1]";

	TecplotSolutionHeader(OutFile_FlowField,CaseName(),varList,varNum,nstep,NLx,NLy);

	writeSolutionData(OutFile_FlowField,PhiFab);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::RHO);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::PHIx);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::PHIy);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::U);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::V);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::P);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::RHOx);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::RHOy);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::FORCEx);
	// writeSolutionData(OutFile_FlowField,MsQ,Ark::FORCEy);
}
void DUGKSAmr::WriteTecplotGrid(int const lev)
{
	int const NLx = Geom(lev).Domain().length(Ark::X);
	int const NLy = Geom(lev).Domain().length(Ark::Y);
	int const NLxPlusUnit = NLx + 1;

	double const real_lo_x = Geom(lev).ProbDomain().lo(Ark::X);
	double const real_lo_y = Geom(lev).ProbDomain().lo(Ark::Y);

	const Real *dx = Geom(lev).CellSize();

	const MultiFab & RhoFab = Data(lev).msv[Ark::RHO];

	ostringstream oss_Grid;

	oss_Grid <<"Nx"  << NLx
			 <<"_Ny" << NLy
			 <<"_"<< CaseName()
			 <<"_AMReX.dat";

	ofstream OutFile_Grid;	
	FileOpen(OutFile_Grid,oss_Grid,"Grid");

	TecplotGridHeader(OutFile_Grid,CaseName(),"X,Y",NLx,NLy);
	
	for(int j = 0;j <= NLy;++j)
	for(int i = 0;i <= NLx;++i)
	{
		OutFile_Grid<<i*dx[0]<<std::endl;
	}
	
	for(int j = 0;j <= NLy;++j)
	for(int i = 0;i <= NLx;++i)
	{
		OutFile_Grid<<j*dx[1]<<std::endl;
	}

	#ifdef _OPENMP
	#pragma omp parallel
	#endif
	for(amrex::MFIter mfi(RhoFab, amrex::TilingIfNotGPU());mfi.isValid();++mfi)
	{
		//- region
		const amrex::Box& bx = mfi.validbox();

		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		LoopPS(lo,hi)
		{
			OutFile_Grid <<j*NLxPlusUnit+i+1<<"  "<<j*NLxPlusUnit+i+2<<"  "
						 <<(j+1)*NLxPlusUnit+i+2<<"  "<<(j+1)*NLxPlusUnit+i+1<<std::endl;
		}
	}
}