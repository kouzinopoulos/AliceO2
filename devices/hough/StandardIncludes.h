/// \file StandardIncludes.h

#ifndef ALICEO2_HOUGH_STANDARDINCLUDES_H_
#define ALICEO2_HOUGH_STANDARDINCLUDES_H_

#include <fstream>
#include <iostream>

#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <vector>

using namespace std;
/* Use these only if absolutely necessary
eg. in inline functions defined in header files */
#define STDCOUT std::cout
#define STDCERR std::cerr
#define STDENDL std::endl
#define STDIF std::ifstream
#define STDOF std::ofstream

typedef char Char_t; ///< Signed Character 1 byte
typedef unsigned char UChar_t; ///< Unsigned Character 1 byte
typedef short Short_t; ///< Signed Short integer 2 bytes
typedef unsigned short UShort_t; ///< Unsigned Short integer 2 bytes
#ifdef R__INT16
typedef long Int_t; ///< Signed integer 4 bytes
typedef unsigned long UInt_t; ///< Unsigned integer 4 bytes
#else
typedef int Int_t; ///< Signed integer 4 bytes
typedef unsigned int UInt_t; ///< Unsigned integer 4 bytes
#endif
#ifdef R__B64
typedef int Seek_t; ///< File pointer
typedef long Long_t; ///< Signed long integer 4 bytes
typedef unsigned long ULong_t; ///< Unsigned long integer 4 bytes
#else
typedef int Seek_t; ///< File pointer
typedef long Long_t; ///< Signed long integer 8 bytes
typedef unsigned long ULong_t; ///< Unsigned long integer 8 bytes
#endif
typedef float Float_t; ///< Float 4 bytes
typedef double Double_t; ///< Float 8 bytes
typedef char Text_t; ///< General string
typedef bool Bool_t; ///< Boolean (0=false, 1=true) (bool)
typedef unsigned char Byte_t; ///< Byte (8 bits)
typedef short Version_t; ///< Class version identifier
typedef const char Option_t; ///< Option string
typedef int Ssiz_t; ///< String size
typedef float Real_t; ///< TVector and TMatrix element type
typedef long long Long64_t; ///< Portable signed long integer 8 bytes
typedef unsigned long long ULong64_t; ///< Portable unsigned long integer 8 bytes

#endif
