set(gridfunctionspacedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/pdelab/gridfunctionspace)
set(gridfunctionspace_HEADERS                      
        blockwiseordering.hh                    
        compositegridfunctionspace.hh           
        compositeorderingutilities.hh           
        constraints.hh                          
        debug.hh                                
        dofinfo.hh                              
        dynamicblockwiseordering.hh             
        genericdatahandle.hh                    
        gridfunctionspace.hh                    
        gridfunctionspaceutilities.hh           
        interpolate.hh                          
        intersectionindexset.hh                 
        leafordering.hh                         
        lexicographicordering.hh                
        localfunctionspace.hh                   
        localfunctionspacetags.hh               
        localvector.hh                          
        nonleaforderingbase.hh                  
        orderingbase.hh                         
        orderinginterface.hh                    
        powercompositegridfunctionspacebase.hh  
        powergridfunctionspace.hh               
        printvector.hh                          
        tags.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

install(FILES ${gridfunctionspace_HEADERS} DESTINATION ${gridfunctionspacedir})