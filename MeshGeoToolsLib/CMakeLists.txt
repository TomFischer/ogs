# Source files
GET_SOURCE_FILES(SOURCES_MeshGeoToolsLib_Utils Utils)
GET_SOURCE_FILES(SOURCES_MeshGeoToolsLib)
set(SOURCES ${SOURCES_MeshGeoToolsLib} ${SOURCES_MeshGeoToolsLib_Utils})

# Create the library
add_library(MeshGeoToolsLib ${SOURCES})

target_link_libraries(MeshGeoToolsLib
    BaseLib
    MathLib
    MeshLib
    GeoLib
)

if(TARGET Boost)
    add_dependencies(MeshGeoToolsLib Boost)
endif()
