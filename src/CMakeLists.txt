configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_library(
	evolution
	STATIC
	Cell.cpp
	FFT.cpp
	Group.cpp
	Infrared.cpp
	Init.cpp
	Lattice.cpp
	Matrix.cpp
	#MatrixExp.cpp
	Measure.cpp
	Random.cpp
	Setup.cpp
)

add_executable(
	jimwlk
	main.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
)


target_link_libraries(
	jimwlk
	PRIVATE
	evolution
	GSL::gsl
	GSL::gslcblas
	PkgConfig::FFTW
)
