set(gridoperatorcommondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridoperator/common)

set(headercheck_IGNORE  assembler.hh)

set(gridoperatorcommon_HEADERS             
        assembler.hh                    
        assemblerutilities.hh           
        gridoperatorutilities.hh        
        localassemblerenginebase.hh     
        timesteppingparameterinterface.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${gridoperatorcommon_HEADERS} DESTINATION ${gridoperatorcommondir})