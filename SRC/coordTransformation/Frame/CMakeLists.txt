
target_sources(OPS_Transform
    PRIVATE
        CorotFrameTransf3d.cpp
#       CorotFrameTransfWarping3d02.cpp
        CorotFrameTransfWarping3d.cpp
        LinearFrameTransf3d.cpp
        PDeltaFrameTransf3d.cpp
    PUBLIC
        CorotFrameTransf3d.h
#       CorotFrameTransfWarping3d02.h
        CorotFrameTransfWarping3d.h
        LinearFrameTransf3d.h
        PDeltaFrameTransf3d.h
)
target_include_directories(OPS_Transform PUBLIC ${CMAKE_CURRENT_LIST_DIR})

