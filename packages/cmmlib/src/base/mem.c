/*! \file mem.c
    \brief Memory allocation/deallocation routines.
    \ingroup base
*/

#include "cmm.h"
#ifndef DOXYGEN_SKIP
#include <string.h>
#endif

static void *mem_alloc(void *base, int ln);
static void mem_free(void *p);
static int _mem_blocksize(void *buf);

#ifndef DOXYGEN_SKIP
/* Turning this on will bypass mem_ routines in favor of LIBC
#define USE_LIBMALLOC 1
*/

/* Should be on for shared memory stuff, not supported for CMMLIB
#define USE_LOCK 1
*/

/* TODO

   Need to add a m_allocz, which clears the allocated buffer.

   Could add some very fast allocators/deallocators for small
      fixed size buffers.

   Add interface to bind and unbind areas:

      mem_bind_area(ptr, l)
      mem_unbind_area(ptr)

*/

#define IS_ALIGNED(p) ((((unsigned long) (p)) & 0x07) == 0)

#define FATL cmm_fatal("Fatal memory error");

static int debug = 0;

#define DEBUG if(debug)

#ifdef USE_LOCK
#define SMPLOCK(_a) smp_lock(_a)
#define SMPUNLOCK(_a) smp_unlock(_a)
#else
#define SMPLOCK(_a)     /* smp_lock(_a) */
#define SMPUNLOCK(_a)   /* smp_unlock(_a) */
#endif

/* 

  The following will work given that pointers and longs are the same size.
  This is true for the alpha (64 bit) and the x86 (32 bit).

  Note that a Block is 64 bytes for the alpha and 32 bytes for the x86.

 */

#define ALIGNPAD 3

typedef unsigned long mem_ulong;

typedef struct
    {
      void *space;              /* pointer back to the owning space */
      void *next;
      void *prev;
      mem_ulong sz;
      mem_ulong check;
      mem_ulong pad[ALIGNPAD];
    }
Block;

#define ALIGNSIZE sizeof(Block)

typedef struct
    {
      void *nextspace;
      Block *firstfree;
      Block *firstused;
      mem_ulong sz;
      mem_ulong avail;
      mem_ulong used;
      mem_ulong peak;
      mem_ulong lfb;            /* This is valid only after a call to cmm_memcheck */
      mem_ulong flags;
    }
Space;

static Space *spacelist;

/* These flags allow the use of heap. The default is unlimited use. */
static int _mem_useheap = 1;
static int _mem_heapused = 0;

#ifdef USE_LARGE_HEAP
static mem_ulong _mem_maxheap = 8000L * 1024L * 1024L;  /* 8GB    */
#else
static mem_ulong _mem_maxheap = 2000 * 1024 * 1024;     /* 2GB    */
#endif

/*

This looks ugly, but it compiles into two
simple instructions: add,and. It produces
the next greatest or equal pointer to _p
which is aligned on a boundary of size _align.
It is expected that _align is a constant, but
this is not essential.

Note that _align must be a power of 2, and
this is NOT checked!

*/

#define alignptr(_p,_align) \
 (void *) ((  ((unsigned long) (_p)) + (_align-1)) & (~(_align-1)))

#define alignsize(_p,_align) \
      ((  ((unsigned long) (_p)) + (_align-1)) & (~(_align-1)))

#define voidptr(_p) ((void *) (_p))
#define charptr(_p) ((unsigned char *) (_p))

#define FREE_CHECK 12345
#define USED_CHECK 31415

#endif

#ifndef USE_LIBC_MALLOC

/*!
   \brief Allocate memory.
   \param l Number of bytes to allocate.

   Unlike the standard malloc, this routine will print
   an error message and exit if there is a failure.

*/
void *cmm_alloc(long l)
{
  void *p;

#ifdef USE_LIBMALLOC
  p = malloc(l);
#else
  p = mem_alloc(NULL, l);
#endif

  if (l < 0)
    {
      cmm_fatal("l negative, cant allocate %ld\n", l);
    }

  if (p == 0)
    {
      cmm_fatal("cant allocate %ld\n", l);
    }

  if (!IS_ALIGNED(p))
    {
      cmm_fatal("malloc returned non-aligned pointer %p\n", p);
    }

  return p;
}

/*!
   \brief Reallocate memory.
   \param pold  Pointer to old buffer.
   \param l Number of bytes to allocate.
*/
void *cmm_realloc(void *pold, long l)
{
  void *p;
  int l1;

  if(pold == NULL)
    return cmm_alloc(l);


#ifdef USE_LIBMALLOC
  p = realloc(pold, l);
#else
  l1 = _mem_blocksize(pold);
  p = cmm_alloc(l);
  if (p)
    {
      if (!IS_ALIGNED(p))
        {
          cmm_fatal("realloc returned non-aligned pointer %p\n", p);
        }
      memcpy(p, pold, MIN(l, l1));
      cmm_free(pold);
    }

#endif

  if (p == 0)
    {
      cmm_fatal("cant (re)allocate %ld\n", l);
    }

  return p;
}

/*!
   \brief Free allocated memory.
   \param p Buffer to deallocate.

   There are serious repercussions if memory allocated with
   cmm_alloc/realloc are freed with any other routine!

*/
void cmm_free(void *p)
{
#ifdef USE_LIBMALLOC
  free(p);
#else
  mem_free(p);
#endif

}

#endif

/* This routine initializes a space */
static void _mem_resetspace(Space * sp, mem_ulong sz, int flags)
{
  Block *bp;

  bp = (Block *) (charptr(sp) + sizeof (Space));
  bp = alignptr(bp, ALIGNSIZE);

  bp->space = sp;
  bp->next = NULL;
  bp->prev = NULL;
  bp->sz = charptr(sp) + sz - charptr(bp);
  bp->sz = alignsize(bp->sz, ALIGNSIZE) - ALIGNSIZE;
  bp->check = FREE_CHECK;

  sp->firstfree = bp;
  sp->firstused = NULL;
  sp->used = 0;
  sp->avail = bp->sz;
  sp->peak = 0;
  sp->flags = flags;
  sp->sz = bp->sz;

}

/* This routine adds a space to the spacelist, with minor checking */
static void _mem_addspace(void *buf, long sz, int flags)
{
  Space *sp = buf;

  if (buf != alignptr(buf, sizeof (r8)))
    {
      cmm_fatal("addspace %p not aligned to double boundary!!\n", buf);
    }

  if (sz < sizeof (Space) + 4096)
    {
      cmm_fatal("mem_addspace error: size too small %ld\n", sz);
    }

  sp->nextspace = spacelist;
  spacelist = sp;

  _mem_resetspace(sp, sz, flags);

}

/* print a block */
static void _mem_printblock(Block * bp)
{
  printf("bp %p\n", bp);
  printf("  space %p\n", bp->space);
  printf("  next  %p\n", bp->next);
  printf("  prev  %p\n", bp->prev);
  printf("  size  %ld\n", bp->sz);
}

/* Set the heap limit */
void mem_limitheap(int nlim)
{
  if (nlim <= 0)
    _mem_useheap = 0;
  else
    {
      _mem_useheap = 1;
      _mem_maxheap = nlim;
    }
}

void mem_resetspace(void *buf)
{
  Space *sp = buf;

  SMPLOCK(&_mem_smplock);
  _mem_resetspace(sp, sp->sz, sp->flags);
  SMPUNLOCK(&_mem_smplock);
}

void mem_addspace(void *buf, long sz, int flags)
{
  SMPLOCK(&_mem_smplock);
  _mem_addspace(buf, sz, flags);
  SMPUNLOCK(&_mem_smplock);
}


void mem_status(void *base, mem_ulong * results)
{
  Space *sp = base;

/*

   results[0]  total size
           1   avail
           2   used
           3   peak used
           4   largest free block
           5   peak largest free block


*/

  if (base == NULL)
    {
      results[0] = 0;
      results[1] = 0;
      results[2] = 0;
      results[3] = 0;

      sp = spacelist;
      while (sp)
        {
          mem_ulong tres[6];

          mem_status(sp, tres);

          results[0] += tres[0];
          results[1] += tres[1];
          results[2] += tres[2];
          results[3] += tres[3];
          sp = sp->nextspace;
        }
      return;
    }

  results[0] = sp->sz;
  results[1] = sp->avail;
  results[2] = sp->used;
  results[3] = sp->peak;
}

/*!
   \brief Preallocate memory
   \param ln Number of blocks to preallocate.

   Note: block size is 1024 bytes.

*/
void cmm_preallocate(int ln)
{
  long lnbytes;
  void *base;

  lnbytes = ln*1024;

  base = malloc(lnbytes);
  _mem_addspace(base, lnbytes, 0);

}

static void *_mem_alloc(void *base, int ln)
{
  Space *sp;
  Block *bp, *bprev, *bnext, *bp1;
  int j;
  void *ret_ptr;

  ret_ptr = NULL;

  if (base == NULL)
    {
      sp = spacelist;
      while (sp)
        {
          void *p;

          p = _mem_alloc(sp, ln);
          if (p)
            {
              ret_ptr = p;
              goto done;
            }

          sp = sp->nextspace;
        }

      if (ln > 64000)
        j = ln + 1024;

      else
        j = 68000;


        printf("Allocating from heap\n");
        printf("  ln %d space %d\n",ln,j);
        cmm_memstat(NULL);


      if (_mem_useheap && (_mem_heapused + j < _mem_maxheap))
        {
          base = malloc(j);
          if (base == NULL)
            goto done;

          _mem_heapused += j;
        }
      else
        goto done;      /* failed. */

      _mem_addspace(base, j, 0);

      ret_ptr = _mem_alloc(base, ln);
      goto done;
    }

  sp = base;

  ln = alignsize(ln, ALIGNSIZE);        /* length in units of block sz */
  ln += 2 * sizeof (Block);     /* allow some extra space */
  j = ln / sizeof (Block);      /* length in units of block size */

  bp = sp->firstfree;

/* first fit algorithm: */
  while (bp)
    {
      if (bp->sz >= ln)
        break;
      bp = bp->next;
    }

  if (bp == NULL)
    goto done;  /* allocation not possible */

  bprev = bp->prev;
  bnext = bp->next;

  if (bp->sz > ln + 4 * ALIGNSIZE)
    {
      /* split the block */

      bp1 = bp + j;     /* new free block */
      bp1->check = FREE_CHECK;
      bp1->next = bnext;
      bp1->prev = bprev;
      bp1->space = sp;
      bp1->sz = bp->sz - j * sizeof (Block);
      bp->sz = j * sizeof (Block);

      /* modify the free queue, replacing bp with bp1 */
      if (bprev)
        bprev->next = bp1;
      else
        sp->firstfree = bp1;


      if (bnext)
        bnext->prev = bp1;

    }
  else
    {
      /* keep the block as is */

      /* unlink the block from the free queue */
      if (bprev)
        bprev->next = bnext;
      else
        sp->firstfree = bnext;
      if (bnext)
        bnext->prev = bprev;

    }

/* link the block into the head of the used queue */
  bnext = sp->firstused;
  bp->next = bnext;
  bp->prev = NULL;
  if (bnext)
    bnext->prev = bp;

  sp->firstused = bp;

/*
     sp->avail -= ln;
     sp->used  += ln;
*/

  sp->avail -= bp->sz;
  sp->used += bp->sz;


  if (sp->used > sp->peak)
    sp->peak = sp->used;

  bp->check = USED_CHECK;
  ret_ptr = bp + 1;

  DEBUG printf("mem_alloc %p bp %p sz %ld ln %d used %ld\n", ret_ptr, bp, bp->sz, ln, sp->used);

done:

  return ret_ptr;
}

static void *mem_alloc(void *base, int ln)
{
  void *ret_ptr;

#ifdef USE_LOCK
  if (_mem_inited != 12345)
    {
      _mem_inited = 12345;
      smp_initlock(&_mem_smplock);
    }
#endif

  SMPLOCK(&_mem_smplock);
  ret_ptr = _mem_alloc(base, ln);
  SMPUNLOCK(&_mem_smplock);

  return ret_ptr;
}


static void _mem_merge(Space * sp, Block * bp0, Block * bp1)
{
  Block *bnext;

  DEBUG printf("mem_merge %p %p\n", bp0, bp1);

  bnext = bp1->next;
  bp0->sz += bp1->sz;
  bp0->next = bnext;
  if (bnext)
    bnext->prev = bp0;
}

static int _mem_blocksize(void *buf)
{
  Block *bp;

  bp = alignptr(buf, ALIGNSIZE);
  if (bp != buf)
    {
      cmm_fatal("pointer error 1 in mem_free %p %p\n", buf, bp);
    }

  bp--;

  return bp->sz - sizeof (Block);
}



static void _mem_free(void *buf)
{
  Block *bp, *bp0, *bprev, *bnext;
  Space *sp;

  bp = alignptr(buf, ALIGNSIZE);
  if (bp != buf)
    {
      cmm_fatal("pointer error 1 in mem_free %p %p\n", buf, bp);
    }

  bp--;

  if (bp->check != USED_CHECK)
    {
      cmm_fatal("pointer error 2 in mem_free\n");
    }

  bp->check = FREE_CHECK;

  sp = bp->space;

  DEBUG printf("mem_free %p bp %p sz %ld used %ld\n", buf, bp, bp->sz, sp->used + sp->sz);

  bprev = bp->prev;
  bnext = bp->next;

/* unlink from the used queue */
  if (bprev)
    bprev->next = bnext;
  else
    sp->firstused = bnext;

  if (bnext)
    bnext->prev = bprev;

  sp->avail += bp->sz;
  sp->used -= bp->sz;

/* link into the free queue. Keep it in sorted order */
  bp0 = sp->firstfree;

  if (bp0 == NULL)
    {
      sp->firstfree = bp;
      bp->next = NULL;
      bp->prev = NULL;
      return;
    }

  while (1)
    {
      if (bp0 > bp || bp0->next == NULL)
        {
          bprev = bp0->prev;

          if (bprev)
            bprev->next = bp;
          else
            sp->firstfree = bp;

          bp->prev = bprev;
          bp->next = bp0;

          bp0->prev = bp;
          break;
        }

      bp0 = bp0->next;
    }

  bprev = bp->prev;
  bnext = bp->next;

  if (bprev)
    {
      if (charptr(bprev) + bprev->sz == charptr(bp))
        {
          /* merge bp and bprev. Then set bp = bprev */
          _mem_merge(sp, bprev, bp);
          bp = bprev;
        }
    }

  if (bnext)
    {
      if (charptr(bp) + bp->sz == charptr(bnext))
        {
          /* merge bp and bnext */
          _mem_merge(sp, bp, bnext);
        }
    }

}

/*

   Free an allocated block of memory

*/
static void mem_free(void *p)
{
  SMPLOCK(&_mem_smplock);
  _mem_free(p);
  SMPUNLOCK(&_mem_smplock);
}

/*!

   \brief Print memory status
   \param base  memory base.

  Note: if base is NULL then all memory bases will be
  checked and printed.

*/
void cmm_memstat(void *base)
{
  Space *sp = base;
  int status;

  if (base == NULL)
    {
      sp = spacelist;
      while (sp)
        {
          cmm_memstat(sp);
          sp = sp->nextspace;
        }
      return;
    }


  status = cmm_memcheck(base);

  printf("MEM: base %p\n", sp);
  printf("MEM:   status %d\n", status);
  printf("MEM:   avail  %ld\n", sp->avail);
  printf("MEM:   used   %ld\n", sp->used);
  printf("MEM:   peak   %ld\n", sp->peak);
  printf("MEM:   lfb    %ld\n", sp->lfb);

/* other info it would be nice to know:

      largest free block
      # free blocks
      # allocated blocks

*/

}

/*

   List all spaces

*/
void mem_list(void *base);
void memlistall_()
{
  mem_list(NULL);
}

/*
    List all the blocks in a space
*/
void mem_list(void *base)
{
  Space *sp = base;
  Block *bp;

  if (base == NULL)
    {
      sp = spacelist;
      while (sp)
        {
          mem_list(sp);
          sp = sp->nextspace;
        }
      return;
    }


  printf("Base %p\n", base);
  printf("Free blocks:\n");
  bp = sp->firstfree;
  while (bp)
    {
      _mem_printblock(bp);
      bp = bp->next;
    }

  printf("Used blocks:\n");
  bp = sp->firstused;
  while (bp)
    {
      _mem_printblock(bp);
      bp = bp->next;
    }

}

/*!

   \brief Check space for consistency.
   \param base Memory base.

   Go through a space and check the allocated and free blocks for consistency.
   Nothing is printed if everything checks out.

*/
int cmm_memcheck(void *base)
{
  int lused, lfree, sz, lfb;
  Space *sp = base;
  Block *bp, *bprev;

  if (base == NULL)
    {
      sp = spacelist;
      while (sp)
        {
          int status = cmm_memcheck(sp);

          if (status)
            return status;
          sp = sp->nextspace;
        }
      return 0;
    }

  sz = sp->sz;

/* check the free list */
  bp = sp->firstfree;
  bprev = NULL;
  lfree = 0;
  lfb = 0;
  while (bp)
    {
      if (bp->space != sp)
        {
          cmm_fatal("memchk error: free block has bad base ptr\n");
        }

      if (bp->check != FREE_CHECK)
        {
          cmm_fatal("memchk error: free block has bad check\n");
        }

      if (bp->prev != bprev)
        {
          cmm_fatal("memchk error: free block has bad prev ptr\n");
        }

      if (bp->sz <= 0)
        {
          cmm_fatal("memchk error: free block has bad size\n");
        }

      if (bp->sz > lfb)
        lfb = bp->sz;

      lfree = lfree + bp->sz;

      if (lfree > sz)
        {
          cmm_fatal("memchk error: free list size is bad %d > %d\n", lfree, sz);
        }

      bprev = bp;
      bp = bp->next;

      if (bp == sp->firstfree)
        {
          cmm_fatal("memchk error: free list is circular\n");
        }
    }

  sp->lfb = lfb;        /* save this, since we measured it */

  bp = sp->firstused;
  bprev = NULL;
  lused = 0;
  while (bp)
    {
      if (bp->space != sp)
        {
          cmm_fatal("memchk error: used block has bad base ptr\n");
        }

      if (bp->check != USED_CHECK)
        {
          cmm_fatal("memchk error: used block has bad check\n");
        }


      if (bp->prev != bprev)
        {
          cmm_fatal("memchk error: used block has bad prev ptr\n");
        }

      if (bp->sz <= 0)
        {
          cmm_fatal("memchk error: used block has bad size\n");
        }

      lused = lused + bp->sz;

      if (lused > sz)
        {
          cmm_fatal("memchk error: used list size is bad\n");
        }


      bprev = bp;
      bp = bp->next;

      if (bp == sp->firstused)
        {
          cmm_fatal("memchk error: used list is circular\n");
        }
    }

  if (lfree + lused != sz)
    {
      cmm_fatal("memchk error: size %d != lfree %d + lused %d.\n", sz, lfree, lused);
    }

  return 0;     /* everything is normal */
}

/*

    This returns the size of an allocated block of memory

*/
int mem_blksize(void *p)
{
  Block *bp;

  bp = p;
  bp--;
  return bp->sz - sizeof (Block);
}
