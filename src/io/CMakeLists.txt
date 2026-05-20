simall_add_module(io
    SOURCES
        Project.cpp VtkXmlWriter.cpp SolutionMonitors.cpp ForcesReport.cpp
        Checkpoint.cpp EnSightGoldWriter.cpp
        MeshFormats.cpp FluentMshReader.cpp Plot3dReader.cpp
        GmshReader.cpp UnvReader.cpp StarCdReader.cpp
        FaceBasedMeshReader.cpp GambitReader.cpp TetgenReader.cpp
        VtuReader.cpp
        CgnsReader.cpp CgnsWriter.cpp TecplotWriter.cpp
        Hdf5ResultStore.cpp ProjectGraph.cpp ProjectSerializer.cpp
        SchematicJson.cpp
        MovieEncoder.cpp
    HEADERS
        Project.hpp ProjectNode.hpp VtkXmlWriter.hpp SolutionMonitors.hpp
        ForcesReport.hpp Checkpoint.hpp EnSightGoldWriter.hpp
        MeshFormats.hpp FluentMshReader.hpp Plot3dReader.hpp
        GmshReader.hpp UnvReader.hpp StarCdReader.hpp
        FaceBasedMeshReader.hpp GambitReader.hpp TetgenReader.hpp
        VtuReader.hpp
        CgnsReader.hpp CgnsWriter.hpp TecplotWriter.hpp
        Hdf5ResultStore.hpp ProjectGraph.hpp ProjectSerializer.hpp
        SchematicJson.hpp
        MovieEncoder.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Solver SimAll::Workbench)

if(SIMALL_HAVE_HDF5)
    target_compile_definitions(simall_io PUBLIC SIMALL_HAVE_HDF5=1)
    target_link_libraries(simall_io PUBLIC hdf5::hdf5)
endif()
if(SIMALL_HAVE_CGNS)
    target_compile_definitions(simall_io PUBLIC SIMALL_HAVE_CGNS=1)
    target_link_libraries(simall_io PUBLIC CGNS::cgns)
endif()
