set(SOURCES
	point.cpp
	euDist.cpp
	database.cpp
	adjMat.cpp
	import.cpp
	graph.cpp)

add_library(support ${SOURCES})
target_link_libraries(support nauty)