#ifndef MAL_SOLVER_H
#define MAL_SOLVER_H

/*
 * Unified constraint-solving entry point.
 *
 * skeleton / ik / curve / constraint each used to expose their own static
 * solve routines. Every solve now goes through mal_solve(), which validates
 * the request, names the owning module on failure and dispatches to the
 * module wrapper. Module headers (skeleton.h ...) document the payload
 * layouts; the wrapper structs below are what mal_solve() takes.
 *
 * Convention: stateful modules keep their instance pointer inside the
 * request payload; resp is output-only (may be NULL when the caller only
 * wants the status code).
 */

typedef enum {
    MAL_MODULE_SKELETON = 0,
    MAL_MODULE_IK = 1,
    MAL_MODULE_CURVE = 2,
    MAL_MODULE_CONSTRAINT = 3,
    MAL_MODULE__COUNT
} mal_module;

typedef enum {
    MAL_OP_SKELETON_UPDATE = 0,
    MAL_OP_SKELETON_ADD = 1,
    MAL_OP_IK_SOLVE = 0,
    MAL_OP_CURVE_EVAL = 0,
    MAL_OP_CURVE_ADD = 1,
    MAL_OP_CHAIN_STEP = 0,
    MAL_OP_CHAIN_INIT = 1
} mal_op;

typedef enum {
    MAL_OK = 0,
    MAL_ERR_PARAM = 1,       /* invalid argument, message names the module */
    MAL_ERR_UNSUPPORTED = 2, /* unknown module/op pair */
    MAL_ERR_STATE = 3        /* valid args, solver cannot satisfy them */
} mal_status;

/* Forward payload and wrapper types; complete definitions live in the
 * module headers (wrappers must be declared after their payloads). */
typedef struct Skeleton mal_skeleton;
typedef struct Curve mal_curve;
typedef struct Chain mal_chain;
typedef struct mal_ik_req mal_ik_req;
typedef struct mal_ik_resp mal_ik_resp;
typedef struct mal_skel_add_req mal_skel_add_req;
typedef struct mal_skel_add_resp mal_skel_add_resp;
typedef struct mal_skel_update_req mal_skel_update_req;
typedef struct mal_skel_add_wrap mal_skel_add_wrap;
typedef struct mal_curve_add_req mal_curve_add_req;
typedef struct mal_curve_add_resp mal_curve_add_resp;
typedef struct mal_curve_eval_req mal_curve_eval_req;
typedef struct mal_curve_eval_resp mal_curve_eval_resp;
typedef struct mal_curve_eval_wrap mal_curve_eval_wrap;
typedef struct mal_curve_add_wrap mal_curve_add_wrap;
typedef struct mal_chain_init_req mal_chain_init_req;
typedef struct mal_chain_step_req mal_chain_step_req;
typedef struct mal_chain_init_wrap mal_chain_init_wrap;
typedef struct mal_chain_step_wrap mal_chain_step_wrap;

typedef struct {
    mal_module module;
    int op;
    const void* req;   /* module-specific request wrapper, see module headers */
    void* resp;        /* output-only; NULL when only the status is needed */
} mal_request;

/*
 * Each linked module registers one dispatch entry. A standalone module
 * binary contains exactly one entry; a consolidated binary (see pose/)
 * supplies a merged table. The entry receives the already-generically-
 * validated request and owns module-level validation.
 */
typedef mal_status (*mal_dispatch_fn)(int op, const void* req, void* resp,
                                      char* why, int whycap);
typedef struct {
    mal_module module;
    mal_dispatch_fn fn;
} mal_dispatch_entry;

/* Returns MAL_OK on success; on failure fills why (if non-NULL) with a
 * message prefixed by the owning module, e.g. "constraint: iterations<0". */
mal_status mal_solve(const mal_request* r, char* why, int whycap);

/* Dispatch table for every module linked into this binary; the last entry
 * has fn == NULL. Standalone module binaries define a one-entry table;
 * the consolidated pose binary defines a four-entry table. */
extern const mal_dispatch_entry mal_dispatch_table[];

const char* mal_module_name(mal_module m);

#endif
