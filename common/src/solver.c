/*
 * Unified constraint-solving entry point.
 *
 * mal_solve() owns the cross-cutting checks (null request, unknown module,
 * null payload). The per-module dispatch entries registered in
 * mal_dispatch_table own every domain-specific validation and write
 * module-prefixed error messages. Because each binary defines its own
 * table, standalone module binaries link without the other modules.
 */
#include "solver.h"
#include "mal_err.h"

const char* mal_module_name(mal_module m) {
    switch (m) {
        case MAL_MODULE_SKELETON:   return "skeleton";
        case MAL_MODULE_IK:         return "ik";
        case MAL_MODULE_CURVE:      return "curve";
        case MAL_MODULE_CONSTRAINT: return "constraint";
        default:                    return "unknown";
    }
}

mal_status mal_solve(const mal_request* r, char* why, int whycap) {
    if (!r)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "solver: null request");
    if (r->module < 0 || r->module >= MAL_MODULE__COUNT)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "solver: module id %d out of range", (int)r->module);
    if (!r->req)
        return mal_fail(why, whycap, MAL_ERR_PARAM,
                        "%s: null request payload",
                        mal_module_name(r->module));

    const mal_dispatch_entry* e = mal_dispatch_table;
    while (e->fn) {
        if (e->module == r->module)
            return e->fn(r->op, r->req, r->resp, why, whycap);
        e++;
    }
    return mal_fail(why, whycap, MAL_ERR_UNSUPPORTED,
                    "solver: module %s not linked into this binary",
                    mal_module_name(r->module));
}
