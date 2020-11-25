//!The definition can be moved to the Header file if C++17 is supported
#include <DUGKSDmVn.H>
using namespace Ark;

//#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
//#                 D2V5                  #
//#.......................................#
const double D2V5::RT = 1.0/3.0;

const double D2V5::Lambda = 1.5; // = 1/(2RT)

const double D2V5::MaxU = 1.0;

const double D2V5::MaSpan = 1.0;//sqrt(3.0*RT)

const double D2V5::A = 1.0/3.0;
const double D2V5::B = 1.0/6.0;


const double D2V5::omega[Q] = {A,B,B,B,B};

const double D2V5::xi_u[Q] = {0,1,0,-1,0};

const double D2V5::xi_v[Q] = {0,0,1,0,-1};

const unsigned D2V5::BB[Q] = {0,3,4,1,2};


//#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
//#                 D2V9                  #
//#.......................................#
const double D2V9::RT = 1.0/3.0;

const double D2V9::Lambda = 1.5; // = 1/(2RT)

const double D2V9::MaxU = 1.0;

const double D2V9::MaSpan = 1.0;//sqrt(3.0*RT)

const double D2V9::A = 4.0/9.0;
const double D2V9::B = 1.0/9.0;
const double D2V9::C = 1.0/36.0;

const double D2V9::omega[Q] = {A,B,B,B,B,C,C,C,C};

const double D2V9::xi_u[Q] = {0,1,0,-1,0,1,-1,-1,1};

const double D2V9::xi_v[Q] = {0,0,1,0,-1,1,1,-1,-1};

const unsigned D2V9::BB[Q] = {0,3,4,1,2,7,8,5,6};

//#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
//#                 D3V7                  #
//#.......................................#
const double D3V7::RT = 1.0/4.0;

const double D3V7::Lambda = 1.0/(2.*D3V7::RT);

const double D3V7::MaxU = 1.0;

const double D3V7::MaSpan = 1.0;//sqrt(3.0*RT)

const double D3V7::A = 1.0/4.0;
const double D3V7::B = 1.0/8.0;

const double D3V7::omega[Q] = {A,B,B,B,B,B,B};

const double D3V7::xi_u[Q] = {0, 1, 0, 0, -1,  0,  0};

const double D3V7::xi_v[Q] = {0, 0, 1, 0,  0, -1,  0};

const double D3V7::xi_w[Q] = {0, 0, 0, 1,  0,  0, -1};

const unsigned D3V7::BB[Q] = {0,4,5,6,1,2,3};

//#¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨#
//#                 D3V15                 #
//#.......................................#
const double D3V15::RT = 1.0/3.0;

const double D3V15::Lambda = 1.5; // = 1/(2RT)

const double D3V15::MaxU = 1.0;

const double D3V15::MaSpan = 1.0;//sqrt(3.0*RT)

const double D3V15::A = 2.0/9.0;
const double D3V15::B = 1.0/9.0;
const double D3V15::C = 1.0/72.0;

const double D3V15::omega[Q] = {A,B,B,B,B,B,B,C,C,C,C,C,C,C,C};

const double D3V15::xi_u[Q] = {0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1};

const double D3V15::xi_v[Q] = {0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1};

const double D3V15::xi_w[Q] = {0, 0, 0, 0, 0, 1,-1, 1, 1, 1, 1,-1,-1,-1,-1};

const unsigned D3V15::BB[Q] = {0,2,1,4,3,6,5,14,13,12,11,10,9,8,7};