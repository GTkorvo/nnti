/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/
/**
 * TRIOS does not have any internal threading, but it should run in a
 * multithreaded environment.  The nthread library provides a single API
 * for locks and atomic counters.
 */

#include "Trios_config.h"

#include <fcntl.h>
#include <sys/stat.h>

#include <assert.h>
#include <unistd.h>

#include <errno.h>
#include <string.h>

#include <semaphore.h>


#include "Trios_logger.h"
#include "Trios_threads.h"

#define _DEBUG_LOCKS_
#undef _DEBUG_LOCKS_

/* set to LOG_UNDEFINED -- log commands will use default level */
log_level thread_debug_level = LOG_UNDEFINED;


int nthread_lock_init(
        nthread_lock_t *lock)
{
    int rc=0;

#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_lock_init: initializing lock(%p), lock->lock(%p)\n", lock, &lock->lock);
    fflush(stderr);
#endif

    lock->name=tempnam("/tmp", "trios.");
    lock->lock=sem_open(lock->name+4, O_CREAT|O_EXCL, 0600, 1);

#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_lock_init: initialized lock(%p), lock->lock(%p)\n", lock, &lock->lock);
    fflush(stderr);
#endif

    return(rc);
}

int nthread_lock(
        nthread_lock_t *lock)
{
    int rc=0;

#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_lock: locking lock(%p), lock->lock(%p)\n", lock, &lock->lock);
    fflush(stderr);
#endif

    rc=sem_wait(lock->lock);

#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_lock: locked lock(%p), lock->lock(%p)\n", lock, &lock->lock);
    fflush(stderr);
#endif

    return(rc);
}

int nthread_unlock(
        nthread_lock_t *lock)
{
    int rc=0;

#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_unlock: unlocking lock(%p), lock->lock(%p)\n", lock, &lock->lock);
    fflush(stderr);
#endif

    rc=sem_post(lock->lock);

#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_unlock: unlocked lock(%p), lock->lock(%p)\n", lock, &lock->lock);
    fflush(stderr);
#endif

    return(rc);
}


int nthread_lock_fini(
        nthread_lock_t *lock)
{
    int rc=0;

#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_lock_fini: finalizing lock(%p), lock->lock(%p)\n", lock, &lock->lock);
    fflush(stderr);
#endif

    rc=sem_close(lock->lock);
    rc=sem_unlink(lock->name);

#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_lock_fini: finalized lock(%p), lock->lock(%p)\n", lock, &lock->lock);
    fflush(stderr);
#endif

    return(rc);
}

int nthread_counter_init(
        nthread_counter_t *c)
{
    int rc=0;
    int64_t t=0;

    log_debug(thread_debug_level, "nthread_counter_init(STUB)");
    rc=nthread_lock_init(&c->lock);
    c->value = 0;

    return rc;
}

int64_t nthread_counter_increment(
        nthread_counter_t *c)
{
    int rc=0;
    int64_t t=0;

    log_debug(thread_debug_level, "nthread_counter_increment(STUB)");

    if (nthread_lock(&c->lock) != 0) {
        log_error(thread_debug_level, "failed to lock the counter lock: %s", strerror(errno));
        goto cleanup;
    }

    t = c->value;
    c->value++;

    if (nthread_unlock(&c->lock) != 0) {
        log_error(thread_debug_level, "failed to unlock the counter lock: %s", strerror(errno));
    }

cleanup:
    return t;
}

int64_t nthread_counter_decrement(
        nthread_counter_t *c)
{
    int64_t t=0;

    log_debug(thread_debug_level, "nthread_counter_decrement(STUB)");

    if (nthread_lock(&c->lock) != 0) {
        log_error(thread_debug_level, "failed to lock the counter lock: %s", strerror(errno));
        goto cleanup;
    }

    t = c->value;
    c->value--;

    if (nthread_unlock(&c->lock) != 0) {
        log_error(thread_debug_level, "failed to unlock the counter lock: %s", strerror(errno));
    }

cleanup:
    return t;
}

int nthread_counter_fini(
        nthread_counter_t *c)
{
    int rc=0;
    int64_t t=0;

    log_debug(thread_debug_level, "nthread_counter_fini(STUB)");
    rc=nthread_lock_fini(&c->lock);
    c->value = 0;

    return rc;
}
