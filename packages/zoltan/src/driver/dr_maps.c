/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/

#include "dr_const.h"
#include "dr_maps_const.h"
#include "dr_util_const.h"
#include "dr_err_const.h"

#define MAP_ALLOC 10

/*
 *  Routines to build elemental communication maps, given a distributed mesh.
 */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int build_elem_comm_maps(int proc, MESH_INFO_PTR mesh)
{
/*
 * Build element communication maps, given a distributed mesh.
 * This routine builds initial communication maps for Chaco input
 * (for Nemesis, initial communication maps are read from the Nemesis file)
 * and rebuilds communication maps after data migration.
 *
 * One communication map per neighboring processor is built.
 * The corresponding maps on neighboring processors
 * must be sorted in the same order, so that neighboring processors do not
 * have to use ghost elements.   For each communication map's pair of
 * processors, the lower-numbered processor determines the order of the
 * elements in the communication map.  The sort key is the elements' global
 * IDs on the lower-number processor; the secondary key is the neighboring
 * elements global IDs.  The secondary key is used when a single element
 * must communicate with more than one neighbor.
 */

int i, j;
ELEM_INFO *elem;
int iadj_elem;
int iadj_proc;
int indx;
int num_alloc_maps;
int max_adj = 0;
int max_adj_per_map;
int cnt, offset;
int *sindex = NULL;
int tmp;
struct map_list_head {
  int map_alloc_size;
  int *glob_id;
  int *elem_id;
  int *side_id;
  int *neigh_id;
} *tmp_maps = NULL, *map = NULL;

  /*
   *  Free the old maps, if they exist.
   */

  if (mesh->ecmap_id != NULL) {
    safe_free((void **) &(mesh->ecmap_id));
    safe_free((void **) &(mesh->ecmap_cnt));
    safe_free((void **) &(mesh->ecmap_elemids));
    safe_free((void **) &(mesh->ecmap_sideids));
    safe_free((void **) &(mesh->ecmap_neighids));
    mesh->necmap = 0;
  }

  /*
   *  Look for off-processor adjacencies.
   *  Loop over all elements 
   */

  num_alloc_maps = MAP_ALLOC;
  mesh->ecmap_id = (int *) malloc(num_alloc_maps * sizeof(int));
  mesh->ecmap_cnt = (int *) malloc(num_alloc_maps * sizeof(int));
  tmp_maps = (struct map_list_head*) malloc(num_alloc_maps 
                                          * sizeof(struct map_list_head));

  if (mesh->ecmap_id == NULL || mesh->ecmap_cnt == NULL || tmp_maps == NULL) {
    Gen_Error(0, "Fatal:  insufficient memory");
    return 0;
  }

  for (i = 0; i < mesh->num_elems; i++) {
    elem = &(mesh->elements[i]);
    for (j = 0; j < elem->adj_len; j++) {

      /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
      if (elem->adj[j] == -1) continue;

      iadj_elem = elem->adj[j];
      iadj_proc = elem->adj_proc[j];

      if (iadj_proc != proc) {
        /* 
         * Adjacent element is off-processor.
         * Add this element to the temporary data structure for 
         * the appropriate neighboring processor.
         */
        if ((indx = in_list(iadj_proc, mesh->necmap, mesh->ecmap_id)) == -1) {
          /*
           * Start a new communication map.
           */

          if (mesh->necmap >= num_alloc_maps) {
            num_alloc_maps += MAP_ALLOC;
            mesh->ecmap_id = (int *) realloc(mesh->ecmap_id,
                                            num_alloc_maps * sizeof(int));
            mesh->ecmap_cnt = (int *) realloc(mesh->ecmap_cnt,
                                             num_alloc_maps * sizeof(int));
            tmp_maps = (struct map_list_head *) realloc(tmp_maps,
                               num_alloc_maps * sizeof(struct map_list_head));
            if (mesh->ecmap_id == NULL || mesh->ecmap_cnt == NULL || 
                tmp_maps == NULL) {
              Gen_Error(0, "Fatal:  insufficient memory");
              return 0;
            }
          }
          mesh->ecmap_id[mesh->necmap] = iadj_proc;
          mesh->ecmap_cnt[mesh->necmap] = 0;
          map = &(tmp_maps[mesh->necmap]);
          map->glob_id  = (int *) malloc(MAP_ALLOC * sizeof(int));
          map->elem_id  = (int *) malloc(MAP_ALLOC * sizeof(int));
          map->side_id  = (int *) malloc(MAP_ALLOC * sizeof(int));
          map->neigh_id = (int *) malloc(MAP_ALLOC * sizeof(int));
          if (map->glob_id == NULL || map->elem_id == NULL || 
              map->side_id == NULL || map->neigh_id == NULL) {
            Gen_Error(0, "Fatal:  insufficient memory");
            return 0;
          }
          map->map_alloc_size = MAP_ALLOC;
          indx = mesh->necmap;
          mesh->necmap++;
        }
        /* Add to map for indx. */
        map = &(tmp_maps[indx]);
        if (mesh->ecmap_cnt[indx] >= map->map_alloc_size) {
          map->map_alloc_size += MAP_ALLOC;
          map->glob_id  = (int *) realloc(map->glob_id, 
                                          map->map_alloc_size * sizeof(int));
          map->elem_id  = (int *) realloc(map->elem_id, 
                                          map->map_alloc_size * sizeof(int));
          map->side_id  = (int *) realloc(map->side_id, 
                                          map->map_alloc_size * sizeof(int));
          map->neigh_id = (int *) realloc(map->neigh_id, 
                                          map->map_alloc_size * sizeof(int));
          if (map->glob_id == NULL || map->elem_id == NULL || 
              map->side_id == NULL || map->neigh_id == NULL) {
            Gen_Error(0, "Fatal:  insufficient memory");
            return 0;
          }
        }       
        tmp = mesh->ecmap_cnt[indx];
        map->glob_id[tmp] = elem->globalID;
        map->elem_id[tmp] = i;
        map->side_id[tmp] = j+1;  /* side is determined by position in
                                          adj array (+1 since not 0-based). */
        map->neigh_id[tmp] = iadj_elem;
        mesh->ecmap_cnt[indx]++;
        max_adj++;
      }
    }
  }

  /* 
   * If no communication maps, don't need to do anything else. 
   */

  if (mesh->necmap > 0) {

    /*
     * Allocate data structure for element communication map arrays.
     */

    mesh->ecmap_elemids  = (int *) malloc(max_adj * sizeof(int));
    mesh->ecmap_sideids  = (int *) malloc(max_adj * sizeof(int));
    mesh->ecmap_neighids = (int *) malloc(max_adj * sizeof(int));


    /*
     * Allocate temporary memory for sort index.
     */
    max_adj_per_map = 0;
    for (i = 0; i < mesh->necmap; i++)
      if (mesh->ecmap_cnt[i] > max_adj_per_map)
        max_adj_per_map = mesh->ecmap_cnt[i];
    sindex = (int *) malloc(max_adj_per_map * sizeof(int));

    cnt = 0;
    for (i = 0; i < mesh->necmap; i++) {

      map = &(tmp_maps[i]);
      for (j = 0; j < mesh->ecmap_cnt[i]; j++)
        sindex[j] = j;

      /*
       * Sort the map so that adjacent processors have the same ordering
       * for the communication.  
       * Assume the ordering of the lower-numbered processor in the pair
       * of communicating processors.
       */

      if (proc < mesh->ecmap_id[i]) 
        sort2_index(mesh->ecmap_cnt[i], map->glob_id, map->neigh_id, sindex);
      else
        sort2_index(mesh->ecmap_cnt[i], map->neigh_id, map->glob_id, sindex);

      /*
       * Copy sorted data into elem map arrays. 
       */

      offset = cnt;
      for (j = 0; j < mesh->ecmap_cnt[i]; j++) {
        mesh->ecmap_elemids[offset]  = map->elem_id[sindex[j]];
        mesh->ecmap_sideids[offset]  = map->side_id[sindex[j]];
        mesh->ecmap_neighids[offset] = map->neigh_id[sindex[j]];
        offset++;
      }

      cnt += mesh->ecmap_cnt[i];
    }
  }

  /* Free temporary data structure. */
  for (i = 0; i < mesh->necmap; i++) {
    safe_free((void **) &(tmp_maps[i].glob_id));
    safe_free((void **) &(tmp_maps[i].elem_id));
    safe_free((void **) &(tmp_maps[i].side_id));
    safe_free((void **) &(tmp_maps[i].neigh_id));
  }
  safe_free((void **) &tmp_maps);
  safe_free((void **) &sindex);

  return 1;
}

