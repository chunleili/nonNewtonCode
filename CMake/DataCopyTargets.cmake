# add_custom_target(CopySPlisHSPlasHShaders
# 	${CMAKE_COMMAND} -E copy_directory
# 	${PROJECT_SOURCE_DIR}/data/shaders
# 	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/shaders
# 	COMMENT "Copying SPlisHSPlasH shaders"
# )
# set_target_properties(CopySPlisHSPlasHShaders PROPERTIES FOLDER "Data copy")

# add_custom_target(CopyEmitterModels
# 	${CMAKE_COMMAND} -E copy_directory
# 	${PROJECT_SOURCE_DIR}/data/emitter_boundary
# 	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/emitter_boundary
# 	COMMENT "Copying SPlisHSPlasH emitter models"
# )
# set_target_properties(CopyEmitterModels PROPERTIES FOLDER "Data copy")

# if(NOT TARGET CopyPBDShaders)
# 	add_custom_target(CopyPBDShaders
# 		${CMAKE_COMMAND} -E copy_directory
# 		${PBD_INCLUDE_DIR}/data/shaders
# 		${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/pbd_shaders
# 		COMMENT "Copying PBD shaders"
# 	)
# 	add_dependencies(CopyPBDShaders Ext_PBD)
# 	set_target_properties(CopyPBDShaders PROPERTIES FOLDER "Data copy")
# else()
# 	message(STATUS "CopyPBDShaders target already exists. Skipping creation of CopyPBDShaders target.")
# endif()

# add_custom_target(CopyImguiFonts
# 	${CMAKE_COMMAND} -E copy
# 	${PROJECT_SOURCE_DIR}/extern/imgui/misc/fonts/Roboto-Medium.ttf
# 	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/fonts/Roboto-Medium.ttf
# 	COMMENT "Copying fonts"
# )
# set_target_properties(CopyImguiFonts PROPERTIES FOLDER "Data copy")


add_custom_target(copy_all_SPlisHSPlasH ALL
    COMMAND ${CMAKE_COMMAND} -E echo "Copying all data files"
)


add_custom_command(
        TARGET copy_all_SPlisHSPlasH  POST_BUILD
		COMMAND 
			${CMAKE_COMMAND} -E copy_directory
			${PROJECT_SOURCE_DIR}/data/shaders
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/shaders
			COMMENT "Copying SPlisHSPlasH shaders"
        )


add_custom_command(
        TARGET copy_all_SPlisHSPlasH  POST_BUILD
		COMMAND 
			${CMAKE_COMMAND} -E copy_directory
			${PROJECT_SOURCE_DIR}/data/emitter_boundary
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/emitter_boundary
			COMMENT "Copying SPlisHSPlasH emitter models"
        )


add_custom_command(
        TARGET copy_all_SPlisHSPlasH  POST_BUILD
		COMMAND 
			${CMAKE_COMMAND} -E copy_directory
			${PBD_INCLUDE_DIR}/data/shaders
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/pbd_shaders
			COMMENT "Copying PBD shaders"
        )

add_custom_command(
        TARGET copy_all_SPlisHSPlasH  POST_BUILD
		COMMAND 
			${CMAKE_COMMAND} -E copy
			${PROJECT_SOURCE_DIR}/extern/imgui/misc/fonts/Roboto-Medium.ttf
			${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/fonts/Roboto-Medium.ttf
			COMMENT "Copying fonts"
        )