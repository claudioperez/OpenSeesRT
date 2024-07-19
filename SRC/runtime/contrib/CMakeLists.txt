add_library(OpenSeesCommunity OBJECT)

add_subdirectory(packages)

add_subdirectory(uniaxial/Parallel)
add_subdirectory(uniaxial/ReinforcingSteel)
add_subdirectory(uniaxial/AlgebraicHysteresis)
# add_subdirectory(uniaxial/SelfCentering)

# add_subdirectory(packages/YieldSurface/)

