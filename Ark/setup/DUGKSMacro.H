#ifndef DUGKS_DVDF_MACRO_H
#define DUGKS_DVDF_MACRO_H

#define ARK_PRINT_ERROR_MESSAGE(OS,FILE,LINE,FUNC) {\
OS << "\n****************Error****************\n"\
<< "file : "<< FILE << '\n'						 \
<< "line : "<< LINE << '\n'						 \
<< "func : "<< FUNC << '\n'						 \
<< "detail : ";}

#define LoopPS(A,B) \
for(int k = A.z;k <= B.z;++k) \
for(int j = A.y;j <= B.y;++j) \
for(int i = A.x;i <= B.x;++i)


#define LoopPS_XY(A,B) \
for(int j = A.y;j <= B.y;++j) \
for(int i = A.x;i <= B.x;++i)

#define LoopPS_YZ(A,B) \
for(int k = A.z;k <= B.z;++k) \
for(int j = A.y;j <= B.y;++j)

#define LoopPS_XZ(A,B) \
for(int k = A.z;k <= B.z;++k) \
for(int i = A.x;i <= B.x;++i)

#define LoopVS(NUM) for(size_t n = 0;n < NUM;++n)

#endif
