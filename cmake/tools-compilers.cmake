#
# Compiler options
#

# Visual studio solution directories
set_property(GLOBAL PROPERTY USE_FOLDERS on)

# Find compiler type
if (CMAKE_CXX_COMPILER_ID MATCHES ".*MSVC.*")
    set(COMPILER_MSVC TRUE)
else ()
    set(COMPILER_MSVC FALSE)
endif ()
if (CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*")
    set(COMPILER_GCC TRUE)
else ()
    set(COMPILER_GCC FALSE)
endif ()
if (CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*")
    set(COMPILER_CLANG TRUE)
else ()
    set(COMPILER_CLANG FALSE)
endif ()

#function (COMPILER_FLAGS TARGET)
    # Enable simultaneous compilation of source files for MSVC
    if (COMPILER_MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
    endif()

    # Compiler warnings
    if (COMPILER_MSVC)
        if (CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
            string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        else ()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
        endif ()
    elseif (COMPILER_GCC OR COMPILER_CLANG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    endif ()

    # Compiler optimization
    if (COMPILER_GCC OR COMPILER_CLANG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
    endif ()
#endfunction ()
