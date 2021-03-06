

FIND_PATH(HEPMC3_INCLUDE_DIR HepMC3/GenEvent.h HINTS /usr/local/include $ENV{HEPMC3DIR}/include)
#FIND_LIBRARY(HEPMC3_LIBRARY NAMES HepMC3 PATHS $ENV{HEPMC3_ROOT_DIR}/lib)
FIND_LIBRARY(HEPMC3_LIB NAMES HepMC3 PATHS /usr/lib /usr/lib/HepMC3 /usr/local/lib $ENV{HEPMC3DIR}/lib)

IF (HEPMC3_INCLUDE_DIR AND HEPMC3_LIB)
   SET(HepMC3_FOUND TRUE)
ENDIF (HEPMC3_INCLUDE_DIR AND HEPMC3_LIB)


IF (HepMC3_FOUND)
   IF (NOT HEPMC3_FIND_QUIETLY)
      MESSAGE(STATUS "Found HepMC3: ${HEPMC3_LIB}")
      MESSAGE(STATUS "Found HepMC3 include: ${HEPMC3_INCLUDE_DIR}")
   ENDIF (NOT HEPMC3_FIND_QUIETLY)
ELSE (HepMC3_FOUND)
   IF (HEPMC3_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find HepMC3. We search first in the normal library paths, then in $HEPMC3IR")
   ELSE(HEPMC3_FIND_REQUIRED)
      IF(NOT HEPMC3_FIND_QUIETLY)
	 MESSAGE(STATUS "Could not find HepMC3.  We search first in the normal library paths, then in $HEPMC3DIR")
      ENDIF(NOT HEPMC3_FIND_QUIETLY)
   ENDIF (HEPMC3_FIND_REQUIRED)
   
ENDIF (HepMC3_FOUND)

