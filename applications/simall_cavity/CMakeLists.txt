add_executable(simall_cavity main.cpp)
target_link_libraries(simall_cavity PRIVATE
    SimAll::Core SimAll::Utilities SimAll::Meshing SimAll::Materials
    SimAll::Solver SimAll::Parallel SimAll::Io)
target_compile_features(simall_cavity PRIVATE cxx_std_20)
