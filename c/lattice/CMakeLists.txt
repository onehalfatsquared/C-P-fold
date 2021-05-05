set(SOURCES
	latticeP.cpp)

add_library(lattice ${SOURCES})
target_link_libraries(lattice genetic)
