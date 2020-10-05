//!The definition can be moved to the Header file if C++17 is supported
#include <DUGKSDmVn.H>
using namespace Ark;

//!definition of D2V9
const double D2V9::RT = 1.0/3.0;

const double D2V9::Lambda = 1.5; // = 1/(2RT)

const double D2V9::MaxU = 1.0;

const double D2V9::MaSpan = 1.0;//sqrt(3.0*RT)

const double D2V9::A = 1.0/4.0;
const double D2V9::B = 1.0/6.0;
const double D2V9::C = 1.0/9.0;
const double D2V9::D = 1.0/12.0;
const double D2V9::E = 1.0/18.0;
const double D2V9::F = 1.0/36.0;

const double D2V9::omega[Q] = {4.0/9.0,C,C,C,C,F,F,F,F};

const double D2V9::xi_u[Q] = {0,1,0,-1,0,1,-1,-1,1};

const double D2V9::xi_v[Q] = {0,0,1,0,-1,1,1,-1,-1};

const int    D2V9::BB[Q] = {0,3,4,1,2,7,8,5,6};