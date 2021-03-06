set(INCLUDE_DIRECTORIES
  ${BASE_INCLUDE_DIRECTORIES}
  ${Boost_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIR}
  ${ALIROOT}/include
  ${SIMPATH}/include/root
  ${CMAKE_SOURCE_DIR}/devices/hough
)

include_directories(SYSTEM ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${Boost_LIBRARY_DIRS}
  ${ROOT_LIBRARY_DIR}
  ${ALIROOT}/lib
  ${SIMPATH}/lib/root
)

link_directories(${LINK_DIRECTORIES})

set(SRCS
  Accumulator.cxx
  BaseTransformer.cxx
  ClusterCollection.cxx
  Draw.cxx
  Hough.cxx
  HoughTrack.cxx
  MaxFinder.cxx
  Track.cxx
  Transform.cxx
  TransformerRow.cxx
)

#add_definitions( -Ddo_mc )

set(DEPENDENCIES
  ${DEPENDENCIES}
   AliHLTTPC HLTbase
)

set(LIBRARY_NAME ALICEHOUGH)

GENERATE_LIBRARY()

Set(Exe_Names
  runHough
)

set(Exe_Source
  runHough.cxx
)

list(LENGTH Exe_Names _length)
math(EXPR _length ${_length}-1)

ForEach(_file RANGE 0 ${_length})
  list(GET Exe_Names ${_file} _name)
  list(GET Exe_Source ${_file} _src)
  set(EXE_NAME ${_name})
  set(SRCS ${_src})
  set(DEPENDENCIES dl Core RIO Gpad Hist ALICEHOUGH boost_system boost_filesystem boost_program_options)
  GENERATE_EXECUTABLE()
EndForEach(_file RANGE 0 ${_length})

