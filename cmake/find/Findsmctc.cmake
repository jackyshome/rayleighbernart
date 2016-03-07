# Try to find the smctc library and header
# Once done this will define
#
# smctc_FOUND           - system has smctc and it can be used
# smctc_INCLUDE_DIRS    - directory where the header file can be found
# smctc_LIBRARY_DIRS    - Path where smctc required libs file can be found
#


set(smctc_FOUND FALSE)

find_path(smctc_INCLUDE_DIR smctc.hh
  ""
  )

set(smctc_INCLUDE_DIRS
    ${smctc_INCLUDE_DIR}
    )

find_path(smctc_LIBRARY_DIR libsmctc.a
  ""
  )

set(smctc_LIBRARY_DIRS
    ${smctc_LIBRARY_DIR}
    )

if(smctc_INCLUDE_DIR AND
    smctc_LIBRARY_DIR
    )
    set(smctc_FOUND TRUE)
endif()