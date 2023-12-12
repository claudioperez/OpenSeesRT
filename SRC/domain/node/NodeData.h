#ifndef NodeData_H
#define NodeData_H
enum class NodeData : int {
  Unknown             = -1,
  Disp                = 1, //
  Vel                 = 2, //
  Accel               = 3, //
  IncrDisp            = 4, // 3
  IncrDeltaDisp       = 5, // 4
  Reaction            = 6, // 7
  UnbalancedLoad      = 7, // 5
       RayleighForces = 8, // 10001
  DisplNorm              , // 10000
  Pressure               , // 10002
  DisplTrial             , // 0
  VelocTrial             , // 1
  AccelTrial             , // 2
  UnbalanceInclInertia   , // 6
  ReactionInclInertia    , // 8
  ReactionInclRayleigh   , // 9
  Empty                  , // 10
  EigenVector            , // 10+
  DisplSensitivity       , // 1000+
  VelocSensitivity       , // 2000+
  AccelSensitivity       , // 3000+
};

#endif // NodeData_H
