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
   RENUMSTR( MNDO,   "MNDO" )
   RENUMSTR( PrincipalAxes, "PrincipalAxes" )
   RENUMSTR( Translate, "Translate" )
   RENUMSTR( Rotate, "Rotate" )
   RENUMSTR( NONE, "NONE" )
   RENUMSTR( TheoryType_end,  "TheoryType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( RotatingType, RotatingTypeStr )
   RENUMSTR( Axis,  "Axis" )
   RENUMSTR( Eular,  "EularAngle" )
   RENUMSTR( RotatingType_end,  "RotatingType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( RotatedObjectType, RotatedObjectTypeStr )
   RENUMSTR( System,  "System" )
   RENUMSTR( Frame,  "Frame" )
   RENUMSTR( RotatedObjectType_end,  "RotatedObjectType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( ShellType, ShellTypeStr )
   RENUMSTR( k,  "k" )
   RENUMSTR( l,  "l" )
   RENUMSTR( m,  "m" )
   RENUMSTR( ShellType_end,  "ShellType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( CartesianType, CartesianTypeStr )
   RENUMSTR( XAxis,  "XAxis" )
   RENUMSTR( YAxis,  "YAxis" )
   RENUMSTR( ZAxis,  "ZAxis" )
   RENUMSTR( CartesianType_end,  "CartesianType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( AzimuthalType, AzimuthalTypeStr )
   RENUMSTR( sAzimuthal,  "s-azimuthal-quantum-number" )
   RENUMSTR( pAzimuthal,  "p-azimuthal-quantum-number" )
   RENUMSTR( dAzimuthal,  "d-azimuthal-quantum-number" )
   RENUMSTR( AzimuthalType_end,  "AzimuthalType_end" )
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

RENUMSTR_BEGIN( STOnGType, STOnGTypeStr )
   RENUMSTR( STO1G,  "STO1G" )
   RENUMSTR( STO2G,  "STO2G" )
   RENUMSTR( STO3G,  "STO3G" )
   RENUMSTR( STO4G,  "STO4G" )
   RENUMSTR( STO5G,  "STO5G" )
   RENUMSTR( STO6G,  "STO6G" )
   RENUMSTR( STOnGType_end,  "STOnGType_end" )
RENUMSTR_END()

// For the definition of the MultipopleType, see appendix in [DT_1977].
RENUMSTR_BEGIN( MultipoleType, MultipoleTypeStr )
   RENUMSTR( sQ,  "q(small Q)" )
   RENUMSTR( Qxx, "Qxx" )
   RENUMSTR( Qyy, "Qyy" )
   RENUMSTR( Qzz, "Qzz" )
   RENUMSTR( Qxz, "Qxz" )
   RENUMSTR( Qyz, "Qyz" )
   RENUMSTR( Qxy, "Qxy" )
   RENUMSTR( mux, "mux" )
   RENUMSTR( muy, "muy" )
   RENUMSTR( muz, "muz" )
   RENUMSTR( MultipoleType_end,  "MultipoleType_end" )
RENUMSTR_END()

}
#endif

