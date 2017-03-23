#
# ZipPack - Mesh Zippering
#

# Sources
file(GLOB_RECURSE SOURCES RELATIVE ${SOURCE_DIR} *.c)
file(GLOB_RECURSE HEADERS RELATIVE ${SOURCE_DIR} *.h)

# Includes
include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ${SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}
)

# Coverage
add_compile_options(-fprofile-arcs -ftest-coverage --coverage)

# Target
set(PROGRAM_NAME zipper)
add_executable(${PROGRAM_NAME} ${SOURCES})
set_property(TARGET ${PROGRAM_NAME} PROPERTY C_STANDARD 90)

# Libraries
target_link_libraries(${PROGRAM_NAME} gcov)
target_link_libraries(${PROGRAM_NAME} m)
target_link_libraries(${PROGRAM_NAME} meschach)
