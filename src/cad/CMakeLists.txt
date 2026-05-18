simall_add_module(cad
    SOURCES
        CadKernel.cpp
        TopologyGraph.cpp
        ShapeHandleInternal.cpp
        StepReader.cpp
        IgesReader.cpp
        StlSurfaceReader.cpp
        ObjReader.cpp
        BrepIo.cpp
        Healing.cpp
        Defeaturing.cpp
        Booleans.cpp
        FeatureEdges.cpp
        Tessellator.cpp
        PersistentIdManager.cpp
    HEADERS
        CadKernel.hpp
        TopologyGraph.hpp
        ShapeHandleInternal.hpp
        StepReader.hpp
        IgesReader.hpp
        StlSurfaceReader.hpp
        ObjReader.hpp
        BrepIo.hpp
        Healing.hpp
        Defeaturing.hpp
        Booleans.hpp
        FeatureEdges.hpp
        Tessellator.hpp
        PersistentIdManager.hpp
    PUBLIC_DEPS
        SimAll::Core SimAll::Utilities)

# OpenCASCADE link target — modern config exposes ${OpenCASCADE_LIBRARIES}.
target_include_directories(simall_cad PRIVATE ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries     (simall_cad PRIVATE ${OpenCASCADE_LIBRARIES})

# ShapeHandleInternal.hpp is an *internal* header — exposing it via the
# public include directory would let other subsystems pull OpenCASCADE
# transitively. Keep it install-excluded; only src/cad/*.cpp may include it.
set_source_files_properties(ShapeHandleInternal.hpp PROPERTIES HEADER_FILE_ONLY TRUE)
