target_sources(SPlisHSPlasH
PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/SurfaceParticles/SurfaceParticles.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SurfaceParticles/SurfaceParticles.h

    ${CMAKE_CURRENT_LIST_DIR}/Viscosity_Casson/Viscosity_Casson.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Viscosity_Casson/Viscosity_Casson.h
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/Coagulation.h	
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/Coagulation.cpp
	${CMAKE_CURRENT_LIST_DIR}/Diffusion/TemperatureDiffusion.h
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/TemperatureDiffusion.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/Melting.h	
    ${CMAKE_CURRENT_LIST_DIR}/Diffusion/Melting.cpp

    # ${CMAKE_CURRENT_LIST_DIR}/MyPartioReader/MyPartioReader.cpp 
	# ${CMAKE_CURRENT_LIST_DIR}/MyPartioReader/MyPartioReader.h

    ${CMAKE_CURRENT_LIST_DIR}/ShapeMatching/ShapeMatching.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ShapeMatching/ShapeMatching.h
    
    ${CMAKE_CURRENT_LIST_DIR}/Plastic/Plastic.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Plastic/Plastic.h

    ${CMAKE_CURRENT_LIST_DIR}/MyTimeStep/MyTimeStep.h
	${CMAKE_CURRENT_LIST_DIR}/MyTimeStep/MyTimeStep.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MyTimeStep/SimulationDataMyTimeStep.h
    ${CMAKE_CURRENT_LIST_DIR}/MyTimeStep/SimulationDataMyTimeStep.cpp

    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton.h
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton.cpp

    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_XSPH.h
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_XSPH.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_standard.h
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_standard.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Bender2017.h
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Bender2017.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Peer2015.h
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Peer2015.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Peer2016.h
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Peer2016.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Weiler2018.h
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Weiler2018.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Takahashi2015.h
    ${CMAKE_CURRENT_LIST_DIR}/NonNewton/NonNewton_Takahashi2015.cpp

    ${CMAKE_CURRENT_LIST_DIR}/Utils/mathUtils.h
    ${CMAKE_CURRENT_LIST_DIR}/Utils/mathUtils.cpp
)

target_include_directories(SPlisHSPlasH PRIVATE ${CMAKE_CURRENT_LIST_DIR})