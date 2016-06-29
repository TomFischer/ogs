add_executable(ComputeNodeAreasFromSurfaceMesh ComputeNodeAreasFromSurfaceMesh.cpp)
set_target_properties(ComputeNodeAreasFromSurfaceMesh PROPERTIES FOLDER Utilities)
target_link_libraries(ComputeNodeAreasFromSurfaceMesh
    MeshLib
    ${OGS_VTK_REQUIRED_LIBS}
)

add_executable(scaleProperty scaleProperty.cpp)
set_target_properties(scaleProperty PROPERTIES FOLDER Utilities)
target_link_libraries(scaleProperty
    MeshLib
    ${OGS_VTK_REQUIRED_LIBS}
)

add_executable(createNeumannBc createNeumannBc.cpp)
set_target_properties(createNeumannBc PROPERTIES FOLDER Utilities)
target_link_libraries(createNeumannBc
    MeshLib
    ${OGS_VTK_REQUIRED_LIBS}
)

if(OGS_BUILD_METIS)
    add_subdirectory(PartitionMesh)
endif()

add_executable( ConvertPolylineToVerticalPolygon
	ConvertPolylineToVerticalPolygon.cpp
	ExtractPolygon.cpp
	ExtractPolygon.h
)

target_link_libraries( ConvertPolylineToVerticalPolygon
	BaseLib
	GeoLib
	MeshLib
	ApplicationsFileIO
	${VTK_LIBRARIES}
)
ADD_VTK_DEPENDENCY(ConvertPolylineToVerticalPolygon)
set_target_properties(ConvertPolylineToVerticalPolygon PROPERTIES FOLDER Utilities)
