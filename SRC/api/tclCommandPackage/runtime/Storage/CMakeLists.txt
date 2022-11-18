add_library(G3_Table)

target_sources(G3_Table
    PRIVATE
        G3_IntMap.c
        G3_StringMap.c
        G3_Table.c
    PUBLIC
        G3_IntMap.h
        G3_StringMap.h
        G3_Table.h
)

