simall_add_module(parallel
    SOURCES
        Parallel.cpp
        MpiContext.cpp
        GhostExchange.cpp
        DomainPartition.cpp
        FieldReducer.cpp
        LoadBalancer.cpp
        NumaPinning.cpp
    HEADERS
        Parallel.hpp
        MpiContext.hpp
        IFieldSynchronizer.hpp
        GhostExchange.hpp
        DomainPartition.hpp
        FieldReducer.hpp
        LoadBalancer.hpp
        NumaPinning.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities SimAll::Meshing)

if(SIMALL_ENABLE_MPI AND MPI_CXX_FOUND)
    target_compile_definitions(simall_parallel PUBLIC SIMALL_HAVE_MPI)
    target_link_libraries(simall_parallel PUBLIC MPI::MPI_CXX)
endif()
if(SIMALL_ENABLE_OPENMP AND OpenMP_CXX_FOUND)
    target_compile_definitions(simall_parallel PUBLIC SIMALL_HAVE_OPENMP)
    target_link_libraries(simall_parallel PUBLIC OpenMP::OpenMP_CXX)
endif()
if(SIMALL_ENABLE_TBB AND TBB_FOUND)
    target_compile_definitions(simall_parallel PUBLIC SIMALL_HAVE_TBB)
    target_link_libraries(simall_parallel PUBLIC TBB::tbb)
endif()
