#ifndef IK_PUB_H
#define IK_PUB_H

#include "solver.h"
#include "mal_vec.h"

struct mal_ik_req {
    mal_v2 root;
    mal_v2 target;
    mal_v2 pole;
    float l1, l2;
};

struct mal_ik_resp {
    int ok;               /* 0 = unreachable / too close (status stays MAL_OK) */
    mal_v2 joint;
    float a1, a2;         /* world angles of bone 1 and bone 2 */
};

mal_status mal_ik_solve(const mal_ik_req* req, mal_ik_resp* resp,
                        char* why, int whycap);

#endif
