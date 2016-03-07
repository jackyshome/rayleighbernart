# Try to find the CImg library and header
# Once done this will define
#
# CImg_FOUND           - system has CImg and it can be used
# CImg_INCLUDE_DIRS    - directory where the header file can be found
#

set(CImg_FOUND FALSE)

find_path(CImg_INCLUDE_DIR CImg.h
  ""
  )

set(CImg_INCLUDE_DIRS
    ${CImg_INCLUDE_DIR}
    ${CImg_INCLUDE_DIR}/plugins
    )

if(CImg_INCLUDE_DIR)
    set(CImg_FOUND TRUE)
endif()
