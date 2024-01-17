#####################################################################################################

# #

# procUniaxialPinching.tcl #

# procedure for activating the pinching material given its parameters in the form of list #

# created NM (nmitra@u.washington.edu) dated : Feb 2002 #

#####################################################################################################

proc procUniaxialPinching { materialTag pEnvelopeStress nEnvelopeStress pEnvelopeStrain nEnvelopeStrain rDisp rForce uForce gammaK gammaD gammaF gammaE damage} {

# add material - command: uniaxialMaterial ...... paramaters as shown

#uniaxialMaterial Pinching4 tag

#### stress1P strain1P stress2P strain2P stress3P strain3P stress4P strain4P

#### stress1N strain1N stress2N strain2N stress3N strain3N stress4N strain4N

#### rDispP rForceP uForceP rDispN rForceN uForceN

#### gammaK1 gammaK2 gammaK3 gammaK4 gammaKLimit

#### gammaD1 gammaD2 gammaD3 gammaD4 gammaDLimit

#### gammaF1 gammaF2 gammaF3 gammaF4 gammaFLimit gammaE $damage

uniaxialMaterial Pinching4 $materialTag [lindex $pEnvelopeStress 0] [lindex $pEnvelopeStrain 0] [lindex $pEnvelopeStress 1] [lindex $pEnvelopeStrain 1] [lindex $pEnvelopeStress 2] [lindex $pEnvelopeStrain 2] [lindex $pEnvelopeStress 3] [lindex $pEnvelopeStrain 3] [lindex $nEnvelopeStress 0] [lindex $nEnvelopeStrain 0] [lindex $nEnvelopeStress 1] [lindex $nEnvelopeStrain 1] [lindex $nEnvelopeStress 2] [lindex $nEnvelopeStrain 2] [lindex $nEnvelopeStress 3] [lindex $nEnvelopeStrain 3] [lindex $rDisp 0] [lindex $rForce 0] [lindex $uForce 0] [lindex $rDisp 1] [lindex $rForce 1] [lindex $uForce 1] [lindex $gammaK 0] [lindex $gammaK 1] [lindex $gammaK 2] [lindex $gammaK 3] [lindex $gammaK 4] [lindex $gammaD 0] [lindex $gammaD 1] [lindex $gammaD 2] [lindex $gammaD 3] [lindex $gammaD 4] [lindex $gammaF 0] [lindex $gammaF 1] [lindex $gammaF 2] [lindex $gammaF 3] [lindex $gammaF 4] $gammaE $damage

}