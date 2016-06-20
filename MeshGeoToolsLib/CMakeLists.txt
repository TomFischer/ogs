# Source files
GET_SOURCE_FILES(SOURCES_MeshGeoToolsLib)
GET_SOURCE_FILES(SOURCES_EDITING MeshEditing)
GET_SOURCE_FILES(SOURCES_QUALITY MeshQuality)
GET_SOURCE_FILES(SOURCES_SEARCH MeshSearch)

set(SOURCES ${SOURCES_MeshGeoToolsLib} ${SOURCES_QUALITY} ${SOURCES_SEARCH}
    ${SOURCES_EDITING})

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
