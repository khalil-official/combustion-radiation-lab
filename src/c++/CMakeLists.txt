add_library(radlib SHARED)
add_library(radlib::radlib ALIAS radlib)

#################### Compile options

#-------- create macro so that source code can access the rcslw data files

add_definitions(-DRCSLW_DATA_DIR=${INSTALLED_RCSLW_DATA_DIR})

#-------- MSVC portability: rad_rcslw.cc uses M_PI, which MSVC only exposes
#-------- in <cmath> when _USE_MATH_DEFINES is defined.

if(MSVC)
    target_compile_definitions(radlib PRIVATE _USE_MATH_DEFINES)
    target_compile_options(radlib PRIVATE /W4)
    # The headers declare no dllexport macros; export every symbol so an
    # import library is generated and consumers can link the DLL.
    set_target_properties(radlib PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    target_compile_options(radlib PRIVATE -Wall -Wextra)
endif()

#################### Local source files

target_sources(radlib 
    PRIVATE
        rad_planck_mean.cc
        rad_wsgg.cc
        rad_rcslw.cc
    PUBLIC    
    FILE_SET radlib_headers
            TYPE HEADERS
            BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" 
            FILES rad.h
                  rad_planck_mean.h
                  rad_wsgg.h
                  rad_rcslw.h
                  multilinear_interpolation.h
              )

############################ Install library

include(GNUInstallDirs)

install(TARGETS  radlib                                     # install library
        EXPORT   radlib_exp                                 # define export var for next comand
        FILE_SET radlib_headers DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")

install(EXPORT      radlib_exp                              # hello export file
        FILE        radlib_exp.cmake  
        NAMESPACE   radlib::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/radlib")

install(FILES       "radlibConfig.cmake"                    # install config files
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/radlib") # edit file radlibConfig.cmake
