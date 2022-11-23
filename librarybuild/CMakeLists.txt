cmake_minimum_required(VERSION 3.11)	

project(										
	msml                       					
	VERSION 0.0                     		
	DESCRIPTION "MS Math Library"      
	LANGUAGES CXX
)

add_library(									
	msml             							
	STATIC										
	src/MSBLAS.cpp									
	src/Real_Vector.cpp
	src/Matrix.cpp             							
) 

target_compile_options(							
	msml                 						
	PRIVATE                  					
	/std:c++17              					
	/W4                   						
)

set_target_properties(
	msml
	PROPERTIES
		ARCHIVE_OUTPUT_DIRECTORY_DEBUG		"${CMAKE_SOURCE_DIR}/lib"
		ARCHIVE_OUTPUT_NAME_DEBUG			"msml_debug"
		ARCHIVE_OUTPUT_DIRECTORY_RELEASE	"${CMAKE_SOURCE_DIR}/lib"
		ARCHIVE_OUTPUT_NAME_RELEASE			"msml_release"
)