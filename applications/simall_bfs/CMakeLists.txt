add_executable(simall_bfs main.cpp)
target_link_libraries(simall_bfs PRIVATE
    SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Materials
    SimAll::Solver SimAll::Turbulence SimAll::Parallel SimAll::Io)
target_compile_features(simall_bfs PRIVATE cxx_std_20)
