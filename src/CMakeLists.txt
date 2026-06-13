configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

if(BUILD_TESTING)
	find_package(GTest QUIET)
	if(NOT GTest_FOUND)
		include(FetchContent)
		FetchContent_Declare(
			googletest
			URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
		)
		FetchContent_MakeAvailable(googletest)
	endif()

	if(TARGET GTest::gtest_main)
		set(MATRIX_TEST_MAIN_TARGET GTest::gtest_main)
	elseif(TARGET gtest_main)
		set(MATRIX_TEST_MAIN_TARGET gtest_main)
	else()
		message(FATAL_ERROR "GoogleTest main target not found")
	endif()

	include(GoogleTest)
endif()

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

if(BUILD_TESTING)
	add_executable(
		matrix_unit_tests
		MatrixTests.cpp
	)
endif()


target_link_libraries(
	jimwlk
	PRIVATE
	evolution
	GSL::gsl
	GSL::gslcblas
	PkgConfig::FFTW
)

if(BUILD_TESTING)
	target_link_libraries(
		matrix_unit_tests
		PRIVATE
		evolution
		${MATRIX_TEST_MAIN_TARGET}
		GSL::gsl
		GSL::gslcblas
		PkgConfig::FFTW
	)

	gtest_discover_tests(matrix_unit_tests)
endif()
