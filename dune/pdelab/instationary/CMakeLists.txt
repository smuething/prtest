set(instationarydir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/instationary)
set(instationary_HEADERS  onestep.hh         
                       pvdwriter.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${instationary_HEADERS} DESTINATION ${instationarydir})