/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __ZOLTAN_TYPES_H
#define __ZOLTAN_TYPES_H

#include <mpi.h>
/* to get PRIdMAX, etc */
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* The default ZOLTAN_ID_TYPE is "unsigned long" but this can be over-ridden on the compile command line.  
 *
 * The type of a Zoltan object global ID is ZOLTAN_ID_TYPE.  A pointer to it is ZOLTAN_ID_PTR.
 *
 * It's decimal type specifier: printf("%" ZOLTAN_ID_SPECIFIER "\n", global_id);
 *
 * A constant of the same type:   ZOLTAN_ID_TYPE global_id = ZOLTAN_ID_CONSTANT(0);
 *
 * We assume the local number of objects fits in a 32 bit integer, but the global number may require
 * the maximum integer width available on the machine.
 *
 * The type of a global count is: intmax_t or uintmax_t
 *
 * It's signed decimal type specifier is:    printf("%" PRIdMAX "\n",globalNum);
 * It's unsigned decimal type specifier is:  printf("%" PRIuMAX "\n",globalNum);
 *
 * A constant of that type is:   intmax_t numGlobal =  INTMAX_C(0)
 *                              uintmax_t numGlobal = UINTMAX_C(0)
 *
 * We don't assume a pointer is the same size as any size of int.  If we want to store
 * a pointer in an int we use types intptr_t or uintptr_t.
 */

#undef ZOLTAN_ID_MPI_TYPE

#ifdef ZOLTAN_ID_TYPE_SHORT

typedef short ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_SHORT
#define ZOLTAN_ID_SPECIFIER  "hd"
#define ZOLTAN_ID_CONSTANT(z)  z

#endif

#ifdef ZOLTAN_ID_TYPE_INT

typedef int ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_INT
#define ZOLTAN_ID_SPECIFIER  "d"
#define ZOLTAN_ID_CONSTANT(z)  z
#endif

#ifdef ZOLTAN_ID_TYPE_LONG

typedef long ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_LONG
#define ZOLTAN_ID_SPECIFIER  "ld"
#define ZOLTAN_ID_CONSTANT(z)  z ## L
#endif

#ifdef ZOLTAN_ID_TYPE_LONG_LONG

typedef long long ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_LONG_LONG
#define ZOLTAN_ID_SPECIFIER  "Ld"
#define ZOLTAN_ID_CONSTANT(z)  z ## LL
#endif

#ifdef ZOLTAN_ID_TYPE_UNSIGNED_SHORT

typedef unsigned short ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_UNSIGNED_SHORT
#define ZOLTAN_ID_SPECIFIER  "hu"
#define ZOLTAN_ID_CONSTANT(z)  z
#endif

#ifdef ZOLTAN_ID_TYPE_UNSIGNED_INT

typedef unsigned int ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_UNSIGNED
#define ZOLTAN_ID_SPECIFIER  "u"
#define ZOLTAN_ID_CONSTANT(z)  z ## U
#endif

#ifdef ZOLTAN_ID_TYPE_UNSIGNED_LONG

typedef unsigned long ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_UNSIGNED_LONG
#define ZOLTAN_ID_SPECIFIER  "lu"
#define ZOLTAN_ID_CONSTANT(z)  z ## UL
#endif

#ifdef ZOLTAN_ID_TYPE_UNSIGNED_LONG_LONG

typedef unsigned long long ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_UNSIGNED_LONG_LONG
#define ZOLTAN_ID_SPECIFIER  "Lu"
#define ZOLTAN_ID_CONSTANT(z)  z ## ULL
#endif

#ifndef ZOLTAN_ID_MPI_TYPE

typedef unsigned int ZOLTAN_ID_TYPE;
#define ZOLTAN_ID_MPI_TYPE  MPI_UNSIGNED
#define ZOLTAN_ID_SPECIFIER  "u"
#define ZOLTAN_ID_CONSTANT(z)  z ## U

#endif

typedef ZOLTAN_ID_TYPE     *ZOLTAN_ID_PTR;

/* 
 * The MPI_Datatype for intmax_t and uintmax_t are figured out at runtime in Zoltan_set_mpi_types.
 */

extern MPI_Datatype  _mpi_intmax_datatype;
extern MPI_Datatype  _mpi_uintmax_datatype;

#define ZOLTAN_INTMAX_MPI_TYPE  _mpi_intmax_datatype
#define ZOLTAN_UINTMAX_MPI_TYPE  _mpi_uintmax_datatype

/*****************************************************************************/
/*
 * Error codes for Zoltan library
 *   ZOLTAN_OK     - no errors
 *   ZOLTAN_WARN   - some warning occurred in Zoltan library;
 *                   application should be able to continue running
 *   ZOLTAN_FATAL  - a fatal error occurred
 *   ZOLTAN_MEMERR - memory allocation failed; with this error, it could be
 *                   possible to try a different, more memory-friendly,
 *                   algorithm.
 */
/*****************************************************************************/
#define ZOLTAN_OK     0
#define ZOLTAN_WARN   1
#define ZOLTAN_FATAL  -1
#define ZOLTAN_MEMERR -2

/*****************************************************************************/
/* Hypergraph query function types 
 */
/*****************************************************************************/
#define ZOLTAN_COMPRESSED_EDGE   1
#define ZOLTAN_COMPRESSED_VERTEX 2

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif /* !__ZOLTAN_TYPES_H */
