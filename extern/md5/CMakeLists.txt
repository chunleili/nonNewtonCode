if(NOT TARGET MD5)
	add_library(MD5 STATIC
		md5.cpp
		md5.h
		)

	set_target_properties(MD5 PROPERTIES FOLDER "External Dependencies")
	set_target_properties(MD5 PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
	set_target_properties(MD5 PROPERTIES RELWITHDEBINFO_POSTFIX ${CMAKE_RELWITHDEBINFO_POSTFIX})
	set_target_properties(MD5 PROPERTIES MINSIZEREL_POSTFIX ${CMAKE_MINSIZEREL_POSTFIX})
	
else()
	message(STATUS "MD5 already exists, using the existing one")
endif()