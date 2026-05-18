# =============================================================================
# Reference plugins (Week 20).  Each is a MODULE library compiled against the
# header-only SDK in plugins/SimAllPluginSdk.hpp.  They are installed under
# ${CMAKE_INSTALL_LIBDIR}/simall/plugins so the runtime can discover them
# without touching system paths.
# =============================================================================
add_subdirectory(simall_curve_postproc)
add_subdirectory(simall_csv_writer)
add_subdirectory(simall_hello_solver)

# Public SDK header — install so out-of-tree plugin authors can include it
# from the installed package.
install(FILES
    IPlugin.hpp
    SimAllPluginSdk.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/simall/plugins
    COMPONENT   devel)
