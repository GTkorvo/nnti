include(AppendSet)

# Explode an ETI set into a variable number of named component fields
FUNCTION(TRIBITS_ETI_EXPLODE fields etivar outvar)
  LIST(LENGTH fields numfields)
  MATH(EXPR N "${numfields}-1")
  SET(result "")
  foreach (i RANGE 0 ${N})
    LIST(GET fields ${i} field)
    LIST(GET ARGV   ${i} var)
    SET(tmp "")
    STRING(REGEX REPLACE ".*${field}=([^ ]*).*" "\\1" tmp ${etivar})
    SET(${result} "${result},${tmp}")
  endforeach()
  SET(${outvar} "${result}" PARENT_SCOPE)
ENDFUNCTION()

# this macro generates the explicit instantiation macros over the cross product of a number of lists
# it is used for the backwards-compatible method of specifying ETI types
FUNCTION(TpetraExpandTypesetProduct etisetvar datascalarlist lolist golist nodelist)
  foreach(ds ${datascalarlist})
    foreach(lo ${lolist})
      foreach(go ${golist})
        foreach(n ${nodelist})
          LIST(APPEND etiset "CS=${ds} DS=${ds} GO=${go} LO=${lo} N=${n}")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
  SET(${etisetvar} "${${etisetvar}};${etiset}" PARENT_SCOPE)
ENDFUNCTION()

# this macro generates the explicit instantiation macros over the cross product of a number of lists
# it is used for the backwards-compatible method of specifying ETI types
FUNCTION(TpetraExpandTypesetDoubleScalarProduct etisetvar computescalarlist datascalarlist lolist golist nodelist)
  foreach(cs ${computescalarlist})
    foreach(ds ${datascalarlist})
      foreach(lo ${lolist})
        foreach(go ${golist})
          foreach(n ${nodelist})
            LIST(APPEND etiset "CS=${cs} DS=${ds} GO=${go} LO=${lo} N=${n}")
          endforeach()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
  SET(${etisetvar} "${${etisetvar}};${etiset}" PARENT_SCOPE)
ENDFUNCTION()

# effectively, a tupled regex, wrapped in a for loop
# given a list of processed excludes (a list of comma-separated tuples of regexes), 
# determine whether processed_inst (a comma-separated tuple of types) is matched
# if the instantiation matches one of the exclusions, result is set true
FUNCTION(TRIBITS_ETI_CHECK_EXCLUSION processed_excludes processed_inst excluded)
  STRING(REPLACE "," ";" "${processed_inst}" processed_inst)
  LIST(LENGTH processed_inst numfields)
  MATH(EXPR "${numfields}-1" NFm1)
  SET(${excluded} OFF PARENT_SCOPE)
  # check to see whether this is excluded or not
  FOREACH(excl ${processed_excludes})
    #IF (${PROJECT}_VERBOSE_CONFIGURE) 
      MESSAGE(STATUS "Checking instantiation ${processed_inst} against exclusion ${excl}")
    #ENDIF()
    STRING(REPLACE "," ";" "${excl}" excl)
    SET(lcltest ON)
    FOREACH(i RANGE 0 NFm1)
      LIST(GET processed_inst ${i} f)
      LIST(GET excl ${i} e)
      IF (NOT ${f} MATCHES ${e})
        SET(lcltest OFF)
        BREAK()
      ENDIF()
    ENDFOREACH()
    IF(lcltest) 
      SET(${excluded} ON PARENT_SCOPE)
      RETURN()
    ENDIF()
  ENDFOREACH()
ENDFUNCTION()

FUNCTION(TRIBITS_ETI_INDEX_MACRO_FIELDS etifields macrofields indexvar)
ENDFUNCTION()

FUNCTION(TRIBITS_GENERATE_ETI_MACROS etifields etisetvar etiexcludelist)
  # we make lists of tuples first, because we want to make sure they don't have duplicates
  # this algorithm is O(N^2) in the number of instantiations in etisetvar
  MATH(EXPR num_macros "(${ARGC}-3)/2")
  IF(${num_macros} EQUAL 0)
    RETURN()
  ENDIF()
  MATH(EXPR num_macros_minus_1 "${num_macros}-1")
  # process macro fields into a list of indices
  FOREACH(m RANGE 0 ${num_macros_minus_1})
    MATH(EXPR m2   "${m}*2")
    MATH(EXPR m2p1 "${m}*2+1")
    LIST(GET ARGN ${m2}   macroarg)
    LIST(GET ARGN ${m2p1} macrovar${m})
    STRING(REGEX REPLACE "^(.*)\\(.*"     "\\1" macroname${m} ${macroarg})
    STRING(REGEX REPLACE "^.*\\((.*)\\)$" "\\1" macroflds${m} ${macroarg})
    # TRIBITS_ETI_INDEX_MACRO_FIELDS(${etifields} ${macrofields${m}} macroindex${m})
    #IF(${PROJECT}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Processing macro ${macroname${m}}(${macroflds${m}}) into variable ${macrovar${m}}")
    #ENDIF()
  ENDFOREACH()

  # process the exclusions once
  FOREACH(excl ${etiexcludelist})
    TRIBITS_ETI_EXPLODE("${Tpetra_ETI_FIELDS}" "${excl}" excl)
    LIST(APPEND processed_excludes "${excl}")
  ENDFOREACH()
  LIST(LENGTH etifields numfields)
  MATH(EXPR NFm1 "${numfields}-1")
  FOREACH(inst ${etisetvar})
    # strip out the types using regex; we're not assuming that they are ordered
    TRIBITS_ETI_EXPLODE("${Tpetra_ETI_FIELDS}" "${inst}" inst)
    # check whether it is on the exclude list
    TRIBITS_ETI_CHECK_EXCLUSION("${processed_excludes}" "${inst}" excluded)
    IF(excluded) 
      #IF (${PROJECT}_VERBOSE_CONFIGURE) 
        MESSAGE(STATUS "-- instantiation excluded!")
      #ENDIF()
    ELSE()
      # append tuple to list
      #list(APPEND tuple_tslgn "${cs},${ds},${lo},${go},${n}")
      #list(APPEND tuple_slgn  "${ds},${lo},${go},${n}")
      #list(APPEND tuple_lgn   "${lo},${go},${n}")
      #list(APPEND tuple_tslg  "${cs},${ds},${lo},${go}")
      #list(APPEND tuple_slg   "${ds},${lo},${go}")
      #list(APPEND tuple_lg    "${lo},${go}")
      #list(APPEND nodelist    "${n}")
    ENDIF()
  endforeach()
  # remove duplicates from lists
  #list(REMOVE_DUPLICATES tuple_tslgn)
  #list(REMOVE_DUPLICATES tuple_slgn)
  #list(REMOVE_DUPLICATES tuple_lgn)
  #list(REMOVE_DUPLICATES tuple_tslg)
  #list(REMOVE_DUPLICATES tuple_slg)
  #list(REMOVE_DUPLICATES tuple_lg)
  #list(REMOVE_DUPLICATES nodelist)
  # build the macro strings
  FOREACH(m RANGE 0 ${NFm1})
    #TRIBITS_ETI_BUILD_MACRO_STRING(... macro)
    set(${macrovar${m}} "${macro}" PARENT_SCOPE)
  ENDFOREACH()
ENDFUNCTION()

MACRO(TRIBITS_PROCESS_ETI)
  include(${PACKAGE_SOURCE_DIR}/cmake/ExplicitInstantiationSupport.cmake)
ENDMACRO()
