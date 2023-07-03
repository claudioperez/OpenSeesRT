# Tags.tcl # UniaxialMaterial tags
set coreTag   1		;# core concrete
set coverTag  2 	;# cover concrete
set steelTag  3		;# steel

set shearTag  4		;# shear limit state material 
set momTag    5         ;# mom-curv hysteretic model
set axialTag  6         ;# elastic axial force-strain model
set momDegTag 7		;# degrading moment-curv hysteretic model
set axialFailTag 8  	;# axial limit state material
set centerSlipTag 9     ;# elastic slip spring

set rigidMatTag 10
set softMatTag  11

# Section tags
set flexSec  1
set shearSec 2
set flexTopSec 3
set flexBotSec 4
set axialSec 5
set shearAxialSec 6
set axialOnlySec 7
set flexShearSec 8
set shearAxialOnlySec 9

# Limit Curve tags
set shearCurveTag 1
set axialCurveTag 2

# Element tags
set bcTag       99


