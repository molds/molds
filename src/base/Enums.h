#ifndef INCLUDED_ENUMS
#define INCLUDED_ENUMS
#ifndef RENUMSTR_BODY
#define RENUMSTR_BODY 0
#endif
#include "../third_parties/rEnumStr/rEnumStr.h"

namespace MolDS_base{

RENUMSTR_BEGIN( TheoryType, TheoryTypeStr )
   RENUMSTR( CNDO2,  "CNDO/2" )
   RENUMSTR( INDO,   "INDO" )
   RENUMSTR( ZINDOS, "ZINDO/S" )
   RENUMSTR( TheoryType_end,  "TheoryType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( ShellType, ShellTypeStr )
   RENUMSTR( k,  "k" )
   RENUMSTR( l,  "l" )
   RENUMSTR( m,  "m" )
   RENUMSTR( ShellType_end,  "ShellType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( OrbitalType, OrbitalTypeStr )
   RENUMSTR( s,  "s" )
   RENUMSTR( py,  "py" )
   RENUMSTR( pz,  "pz" )
   RENUMSTR( px,  "px" )
   RENUMSTR( dxy,  "dxy" )
   RENUMSTR( dyz,  "dyz" )
   RENUMSTR( dzz,  "dzz" )
   RENUMSTR( dzx,  "dzx" )
   RENUMSTR( dxxyy,  "dxxyy" )
   RENUMSTR( OrbitalType_end,  "OrbitalType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( AtomType, AtomTypeStr )
   RENUMSTR( H,  "H" )
   RENUMSTR( He, "He" )
   RENUMSTR( Li, "Li" )
   RENUMSTR( Be, "Be" )
   RENUMSTR( B,  "B" )
   RENUMSTR( C,  "C" )
   RENUMSTR( N,  "N" )
   RENUMSTR( O,  "O" )
   RENUMSTR( F,  "F" )
   RENUMSTR( Ne,  "Ne" )
   RENUMSTR( Na,  "Na" )
   RENUMSTR( Mg,  "Mg" )
   RENUMSTR( Al,  "Al" )
   RENUMSTR( Si,  "Si" )
   RENUMSTR( P,  "P" )
   RENUMSTR( S,  "S" )
   RENUMSTR( Cl,  "Cl" )
   RENUMSTR( Ar,  "Ar" )
   RENUMSTR( AtomType_end,  "AtomType_end" )
RENUMSTR_END()

}
#endif

