#ifndef SKELETON_PUB_H
#define SKELETON_PUB_H

#include "solver.h"
#include "mal_vec.h"

#define MAL_SKEL_MAX_BONES 16

struct Skeleton {
    int parent[MAL_SKEL_MAX_BONES];
    mal_v2 local[MAL_SKEL_MAX_BONES];
    float local_ang[MAL_SKEL_MAX_BONES];
    mal_v2 world[MAL_SKEL_MAX_BONES];
    float world_ang[MAL_SKEL_MAX_BONES];
    int count;
};

struct mal_skel_add_req {
    int parent;
    float x, y, ang;
};

struct mal_skel_add_resp {
    int index;
};

/* Stateful ops bundle the instance with their inputs. */
struct mal_skel_update_req { mal_skeleton* s; };
struct mal_skel_add_wrap { mal_skeleton* s; struct mal_skel_add_req in; };

typedef struct {
    mal_v2 world[MAL_SKEL_MAX_BONES];
    float ang[MAL_SKEL_MAX_BONES];
    int count;
} mal_skel_pose;

/* Lifetime / direct wrappers (also reachable through mal_solve). */
void mal_skeleton_init(mal_skeleton* s);
int mal_skeleton_capacity(void);
mal_status mal_skeleton_add(mal_skeleton* s, const mal_skel_add_req* req,
                            mal_skel_add_resp* resp, char* why, int whycap);
mal_status mal_skeleton_update(mal_skeleton* s, char* why, int whycap);
mal_status mal_skeleton_snapshot(const mal_skeleton* s, mal_skel_pose* out,
                                 char* why, int whycap);

#endif
