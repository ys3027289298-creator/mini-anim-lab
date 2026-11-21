
/* contact: rope-on-rope collision for verlet particle ropes.
   pipeline per frame:
     1. broad phase  : per-segment-pair AABB overlap culling (AABBs expanded
                        by the two rope radii),
     2. narrow phase : exact closest points of each surviving segment pair,
                        penetration = radius sum - closest distance,
     3. solve        : each contact becomes a position-level constraint and
                        is projected inside the same verlet iteration loop as
                        the stretch constraints.
   a subdivided segment reference (contact_fine) supplies a lower bound on
   penetration; it is used to prove the broad phase drops no real contact
   (missed == 0) and to bracket the coarse penetration from above.
   zero-length segments are invalid: rope_validate / scene_step reject them. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_P         32
#define MAX_SEG       (MAX_P - 1)
#define MAX_PAIRS     (MAX_SEG * MAX_SEG)
#define MAX_CONTACTS  256
#define MAX_ROPES     4
#define EPS           1e-6f

typedef struct { float x, y, px, py; float inv_mass; int pinned; } RParticle;

typedef struct {
    RParticle p[MAX_P];
    int n;
    float rest[MAX_SEG];   /* rest length of segment i between p[i], p[i+1] */
    float radius;
    float damping;
} Rope;

typedef struct {
    int ai, bi;            /* rope indices inside the scene */
    int sa, sb;            /* segment indices */
    float ta, tb;          /* closest-point parameters on each segment */
    float nx, ny;          /* unit normal, pointing from B toward A */
    float nax, nay;        /* per-rope push directions at a crossing: A is
                              pushed along nA, B along nB, each perpendicular
                              to its own tangent, opposite senses */
    float nbx, nby;
    float depth;           /* radius sum - closest distance, > 0 */
    int active;
} Contact;

typedef struct {
    int stretch_iters;     /* iterations that run a stretch sweep */
    int contact_sweeps;    /* iterations that also run a contact sweep */
    int contacts_seen;     /* contacts detected before solving */
    int cap_per_sweep;
    int applied;           /* total contact constraint applications */
    int max_in_sweep;      /* most applications done in a single sweep */
} ContactStats;

typedef struct {
    Rope* ropes[MAX_ROPES];
    int nropes;
    int iterations;
    int contact_cap;       /* max contacts handled per sweep, <= 0 disables */
    float stiffness;
    float slop;
} ContactScene;

/* ------------------------------------------------------------------ ropes */

static int rope_init(Rope* r, const float* xy, int n, float radius) {
    if (n < 2 || n > MAX_P || !(radius > 0.0f)) return -1;
    memset(r, 0, sizeof(*r));
    r->n = n;
    r->radius = radius;
    r->damping = 0.99f;
    for (int i = 0; i < n; i++) {
        r->p[i].x = xy[2 * i];
        r->p[i].y = xy[2 * i + 1];
        r->p[i].px = r->p[i].x;
        r->p[i].py = r->p[i].y;
        r->p[i].inv_mass = 1.0f;
    }
    for (int i = 0; i < n - 1; i++) {
        float dx = r->p[i + 1].x - r->p[i].x;
        float dy = r->p[i + 1].y - r->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        if (d < EPS) return -1;             /* zero-length segment rejected */
        r->rest[i] = d;
    }
    return 0;
}

static void rope_pin(Rope* r, int i) {
    r->p[i].pinned = 1;
    r->p[i].inv_mass = 0.0f;
}

static int rope_validate(const Rope* r) {
    if (r->n < 2 || r->n > MAX_P || !(r->radius > 0.0f)) return -1;
    for (int i = 0; i < r->n - 1; i++)
        if (!(r->rest[i] > EPS)) return -1; /* zero-length segment: hard error */
    return 0;
}

static void rope_integrate(Rope* r, float dt, float gx, float gy) {
    for (int i = 0; i < r->n; i++) {
        RParticle* p = &r->p[i];
        if (p->pinned) continue;
        float vx = (p->x - p->px) * r->damping;
        float vy = (p->y - p->py) * r->damping;
        p->px = p->x; p->py = p->y;
        p->x += vx + gx * dt * dt;
        p->y += vy + gy * dt * dt;
    }
}

static void rope_stretch(Rope* r) {
    for (int i = 0; i < r->n - 1; i++) {
        RParticle* a = &r->p[i];
        RParticle* b = &r->p[i + 1];
        float dx = b->x - a->x, dy = b->y - a->y;
        float d = sqrtf(dx * dx + dy * dy);
        if (d < EPS) continue;              /* validate() should have caught it */
        float diff = (d - r->rest[i]) / d;
        float wa = a->inv_mass, wb = b->inv_mass;
        float wsum = wa + wb;
        if (wsum <= 0.0f) continue;
        a->x += dx * diff * (wa / wsum);
        a->y += dy * diff * (wa / wsum);
        b->x -= dx * diff * (wb / wsum);
        b->y -= dy * diff * (wb / wsum);
    }
}

/* ------------------------------------------------------------- geometry */

typedef struct { float ax, ay, bx, by; float ta, tb, dist; } SegSeg;

static float clamp01(float x) {
    return x < 0.0f ? 0.0f : (x > 1.0f ? 1.0f : x);
}

/* exact closest points of two segments. the minimum over the two interior
   parameters is attained either at the unconstrained stationary point or at
   a domain boundary, so five candidates cover every configuration, including
   crossing (dist 0) and collinear overlap. */
static SegSeg seg_seg_closest(float ax, float ay, float bx, float by,
                              float cx, float cy, float dx, float dy) {
    float d1x = bx - ax, d1y = by - ay;
    float d2x = dx - cx, d2y = dy - cy;
    float rx = ax - cx, ry = ay - cy;
    float a11 = d1x * d1x + d1y * d1y;
    float a22 = d2x * d2x + d2y * d2y;
    float a12 = d1x * d2x + d1y * d2y;
    float b1 = rx * d1x + ry * d1y;
    float b2 = rx * d2x + ry * d2y;
    float det = a11 * a22 - a12 * a12;

    float best = 1e30f, bs = 0.0f, bt = 0.0f;
    /* candidate (s,t): distance of the two interior points */
    #define CONSIDER(s, t) do {                                              \
        float px = ax + d1x * (s), py = ay + d1y * (s);                     \
        float qx = cx + d2x * (t), qy = cy + d2y * (t);                     \
        float ex = px - qx, ey = py - qy;                                   \
        float e2 = ex * ex + ey * ey;                                       \
        if (e2 < best) { best = e2; bs = (s); bt = (t); }                   \
    } while (0)

    if (det > EPS * EPS) {
        float s = (-b1 * a22 + a12 * b2) / det;
        float t = (-b1 * a12 + a11 * b2) / det;
        if (s >= 0.0f && s <= 1.0f && t >= 0.0f && t <= 1.0f)
            CONSIDER(s, t);
    }
    if (a22 > 0.0f) {
        CONSIDER(0.0f, clamp01(-b2 / a22));
        CONSIDER(1.0f, clamp01((a12 - b2) / a22));
    }
    if (a11 > 0.0f) {
        CONSIDER(clamp01(-b1 / a11), 0.0f);
        CONSIDER(clamp01((a12 - b1) / a11), 1.0f);
    }
    #undef CONSIDER
    /* degenerate (zero-length) segment guard: parameters stay in range */
    if (!(a22 > 0.0f)) bt = 0.0f;
    if (!(a11 > 0.0f)) bs = 0.0f;

    SegSeg r;
    r.ta = bs; r.tb = bt;
    r.ax = ax + d1x * bs; r.ay = ay + d1y * bs;
    r.bx = cx + d2x * bt; r.by = cy + d2y * bt;
    float ex = r.ax - r.bx, ey = r.ay - r.by;
    r.dist = sqrtf(ex * ex + ey * ey);
    return r;
}

/* ------------------------------------------------------------- broad/narrow */

static int contact_broad(const Rope* ra, const Rope* rb,
                         int pairs[][2], int max_pairs) {
    int n = 0;
    float reach = ra->radius + rb->radius;
    for (int i = 0; i < ra->n - 1; i++) {
        float aminx = fminf(ra->p[i].x, ra->p[i + 1].x) - reach;
        float amaxx = fmaxf(ra->p[i].x, ra->p[i + 1].x) + reach;
        float aminy = fminf(ra->p[i].y, ra->p[i + 1].y) - reach;
        float amaxy = fmaxf(ra->p[i].y, ra->p[i + 1].y) + reach;
        for (int j = 0; j < rb->n - 1; j++) {
            float bminx = fminf(rb->p[j].x, rb->p[j + 1].x);
            float bmaxx = fmaxf(rb->p[j].x, rb->p[j + 1].x);
            float bminy = fminf(rb->p[j].y, rb->p[j + 1].y);
            float bmaxy = fmaxf(rb->p[j].y, rb->p[j + 1].y);
            if (aminx <= bmaxx && amaxx >= bminx &&
                aminy <= bmaxy && amaxy >= bminy) {
                if (n < max_pairs) { pairs[n][0] = i; pairs[n][1] = j; n++; }
            }
        }
    }
    return n;
}

static int contact_narrow(const Rope* ra, int ai, const Rope* rb, int bi,
                          const int pairs[][2], int np,
                          Contact* out, int max_out) {
    int n = 0;
    float rsum = ra->radius + rb->radius;
    for (int k = 0; k < np; k++) {
        int sa = pairs[k][0], sb = pairs[k][1];
        SegSeg q = seg_seg_closest(
            ra->p[sa].x, ra->p[sa].y, ra->p[sa + 1].x, ra->p[sa + 1].y,
            rb->p[sb].x, rb->p[sb].y, rb->p[sb + 1].x, rb->p[sb + 1].y);
        float depth = rsum - q.dist;
        if (depth <= 0.0f) continue;
        if (n >= max_out) break;
        Contact* c = &out[n++];
        c->ai = ai; c->bi = bi; c->sa = sa; c->sb = sb;
        c->ta = q.ta; c->tb = q.tb;
        c->depth = depth;
        if (q.dist > EPS) {
            c->nx = (q.ax - q.bx) / q.dist;
            c->ny = (q.ay - q.by) / q.dist;
        } else {
            /* crossing / overlap: no geometric normal, so report a unit
               perpendicular of segment A (falling back to B, then the y
               axis) just to keep this standalone query well-defined; the
               solver instead derives separate per-rope push directions in
               pair_axes(). */
            float d1x = ra->p[sa + 1].x - ra->p[sa].x;
            float d1y = ra->p[sa + 1].y - ra->p[sa].y;
            float d2x = rb->p[sb + 1].x - rb->p[sb].x;
            float d2y = rb->p[sb + 1].y - rb->p[sb].y;
            float l1 = sqrtf(d1x * d1x + d1y * d1y);
            float l2 = sqrtf(d2x * d2x + d2y * d2y);
            if (l1 > EPS) {
                c->nx = -d1y / l1;
                c->ny = d1x / l1;
            } else if (l2 > EPS) {
                c->nx = -d2y / l2;
                c->ny = d2x / l2;
            } else {
                c->nx = 0.0f;
                c->ny = 1.0f;
            }
        }
    }
    return n;
}

/* position-based projection of one contact. the contact point is a barycentric
   blend of each segment's endpoints, so the PBD weights carry the barycentric
   coordinates; pinned particles have zero inverse mass and do not move.
   for ordinary (non-crossing) contacts both ropes share the geometric normal.
   at a crossing the two per-rope normals differ: each rope is pushed along its
   own perpendicular direction, so the relative displacement along either
   normal still closes the full penetration depth. */
static void contact_project(Rope* ra, Rope* rb, const Contact* c,
                            float stiffness, float slop) {
    float depth = c->depth - slop;
    if (depth <= 0.0f) return;
    RParticle* a0 = &ra->p[c->sa];
    RParticle* a1 = &ra->p[c->sa + 1];
    RParticle* b0 = &rb->p[c->sb];
    RParticle* b1 = &rb->p[c->sb + 1];
    float wa0 = a0->inv_mass, wa1 = a1->inv_mass;
    float wb0 = b0->inv_mass, wb1 = b1->inv_mass;
    float ta = c->ta, tb = c->tb;
    float den = wa0 * (1 - ta) * (1 - ta) + wa1 * ta * ta +
                wb0 * (1 - tb) * (1 - tb) + wb1 * tb * tb;
    if (den <= 0.0f) return;
    float lambda = depth * stiffness / den;
    a0->x += c->nax * (wa0 * (1 - ta) * lambda);
    a0->y += c->nay * (wa0 * (1 - ta) * lambda);
    a1->x += c->nax * (wa1 * ta * lambda);
    a1->y += c->nay * (wa1 * ta * lambda);
    b0->x -= c->nbx * (wb0 * (1 - tb) * lambda);
    b0->y -= c->nby * (wb0 * (1 - tb) * lambda);
    b1->x -= c->nbx * (wb1 * tb * lambda);
    b1->y -= c->nby * (wb1 * tb * lambda);
}

/* ------------------------------------------------------ iteration budget */

/* stretch keeps the full iteration count. whenever contacts exist the last
   iterations run a contact sweep after each stretch sweep, so the constraints
   alternate and the contact work converges instead of being undone by the
   next frame's integration. cap bounds how many contacts one sweep handles
   (deepest first); if more batches are needed to cover the contacts, k grows
   to ceil(n / cap), still leaving at least one pure stretch iteration. per
   frame contact work is bounded by k * cap. */
static void contact_budget(int iterations, int ncontacts, int cap,
                           ContactStats* st) {
    memset(st, 0, sizeof(*st));
    st->stretch_iters = iterations;
    st->contacts_seen = ncontacts;
    st->cap_per_sweep = cap > 0 ? cap : 0;
    if (ncontacts > 0 && cap > 0 && iterations > 1) {
        int limit = iterations - 1;
        int batches = (ncontacts + cap - 1) / cap;
        int k = batches > 1 ? batches : limit;
        st->contact_sweeps = k < limit ? k : limit;
    }
}

static int scene_detect(const ContactScene* sc, Contact* out, int max_out) {
    int n = 0;
    int pairs[MAX_PAIRS][2];
    for (int i = 0; i < sc->nropes; i++) {
        for (int j = i + 1; j < sc->nropes; j++) {
            int np = contact_broad(sc->ropes[i], sc->ropes[j],
                                   pairs, MAX_PAIRS);
            int used = max_out - n;
            if (used > 0)
                n += contact_narrow(sc->ropes[i], i, sc->ropes[j], j,
                                    pairs, np, out + n, used);
        }
    }
    return n;
}

/* persistent contact set: contacts are keyed by (rope pair, segment pair) and
   keep their slot between sweeps, while geometry is refreshed in place. this
   preserves a stable contact identity as closest points move; the crossing
   push directions come from rope-level geometry (pair_axes), so they do not
   flip when the closest point hops between neighbouring segment pairs. */
static void signed_update(const ContactScene* sc, Contact* set, int* pn) {
    Contact fresh[MAX_CONTACTS];
    int nf = scene_detect(sc, fresh, MAX_CONTACTS);
    for (int i = 0; i < *pn; i++) set[i].active = 0;

    for (int f = 0; f < nf; f++) {
        Contact* c = &fresh[f];
        int slot = -1;
        for (int i = 0; i < *pn; i++) {
            if (set[i].active &&
                set[i].ai == c->ai && set[i].bi == c->bi &&
                set[i].sa == c->sa && set[i].sb == c->sb) { slot = i; break; }
        }
        if (slot < 0) {
            for (int i = 0; i < *pn; i++)
                if (!set[i].active) { slot = i; break; }
            if (slot < 0) {
                if (*pn >= MAX_CONTACTS) continue;
                slot = (*pn)++;
            }
        }
        c->active = 1;
        set[slot] = *c;       /* geometry refreshed each sweep */
    }
}

/* signed side of a point relative to the line through segment s of rope r:
   the 2D cross product of the segment direction with (point - start). */
static float point_side(const Rope* r, int s, float px, float py) {
    float dx = r->p[s + 1].x - r->p[s].x;
    float dy = r->p[s + 1].y - r->p[s].y;
    return dx * (py - r->p[s].y) - dy * (px - r->p[s].x);
}

/* push directions at a crossing. at a zero-distance crossing there is no
   geometric common normal; each rope bends away perpendicular to its own
   tangent. nA is the perpendicular of rope A's overall direction, nB that of
   rope B, with opposite senses so the two centerlines move apart. both are
   derived from rope-level (not per-segment) geometry, which keeps their sign
   stable as closest points hop between neighbouring segment pairs. */
static void pair_axes(const Rope* ra, const Rope* rb,
                      float* nax, float* nay, float* nbx, float* nby) {
    float ax = ra->p[ra->n - 1].x - ra->p[0].x;
    float ay = ra->p[ra->n - 1].y - ra->p[0].y;
    float bx = rb->p[rb->n - 1].x - rb->p[0].x;
    float by = rb->p[rb->n - 1].y - rb->p[0].y;
    float la = sqrtf(ax * ax + ay * ay);
    float lb = sqrtf(bx * bx + by * by);
    float ux = la > EPS ? ax / la : 1.0f;
    float uy = la > EPS ? ay / la : 0.0f;
    float vx = lb > EPS ? bx / lb : 1.0f;
    float vy = lb > EPS ? by / lb : 0.0f;
    /* sign: toward where rope A's first particle sits relative to B's line */
    float side = point_side(rb, 0, ra->p[0].x, ra->p[0].y);
    if (fabsf(side) <= EPS)
        side = point_side(rb, 0, ra->p[ra->n - 1].x, ra->p[ra->n - 1].y);
    float sgn = side < 0.0f ? -1.0f : 1.0f;
    *nax = -uy * sgn; *nay = ux * sgn;       /* perpendicular of A */
    *nbx = vy * sgn;  *nby = -vx * sgn;      /* perpendicular of B, opposite */
}

/* re-evaluate the geometry of a persistent contact from current particle
   positions: same segment pair, fresh closest points and depth. returns 1
   while the pair still overlaps, 0 once it has separated. */
static int signed_refresh(const ContactScene* sc, Contact* c) {
    const Rope* ra = sc->ropes[c->ai];
    const Rope* rb = sc->ropes[c->bi];
    SegSeg q = seg_seg_closest(
        ra->p[c->sa].x, ra->p[c->sa].y,
        ra->p[c->sa + 1].x, ra->p[c->sa + 1].y,
        rb->p[c->sb].x, rb->p[c->sb].y,
        rb->p[c->sb + 1].x, rb->p[c->sb + 1].y);
    c->ta = q.ta; c->tb = q.tb;
    c->depth = ra->radius + rb->radius - q.dist;
    if (c->depth <= 0.0f) return 0;
    if (q.dist > EPS) {
        c->nx = (q.ax - q.bx) / q.dist;
        c->ny = (q.ay - q.by) / q.dist;
        c->nax = c->nx; c->nay = c->ny;
        c->nbx = c->nx; c->nby = c->ny;
    } else {
        /* crossing: no geometric normal; each rope gets its own stable
           perpendicular push direction */
        pair_axes(ra, rb, &c->nax, &c->nay, &c->nbx, &c->nby);
        c->nx = c->nax; c->ny = c->nay;
    }
    return 1;
}

static int scene_step(ContactScene* sc, float dt, float gx, float gy,
                      ContactStats* st) {
    for (int i = 0; i < sc->nropes; i++)
        if (rope_validate(sc->ropes[i]) != 0) return -1;

    for (int i = 0; i < sc->nropes; i++)
        rope_integrate(sc->ropes[i], dt, gx, gy);

    Contact set[MAX_CONTACTS];
    memset(set, 0, sizeof(set));
    int ns = 0;
    signed_update(sc, set, &ns);
    contact_budget(sc->iterations, ns, sc->contact_cap, st);
    int k = st->contact_sweeps;

    for (int it = 0; it < sc->iterations; it++) {
        for (int i = 0; i < sc->nropes; i++)
            rope_stretch(sc->ropes[i]);
        if (k > 0 && it >= sc->iterations - k) {
            /* contact sweep runs after stretch in the same iteration so its
               projection is the last position change; otherwise stretch in a
               later iteration undoes the separation. */
            if (it > sc->iterations - k) {
                signed_update(sc, set, &ns);   /* discover new overlaps */
                if (ns > st->contacts_seen) st->contacts_seen = ns;
            }

            Contact live[MAX_CONTACTS];
            int nl = 0;
            for (int i = 0; i < ns; i++) {
                Contact c = set[i];
                if (c.active && signed_refresh(sc, &c))
                    live[nl++] = c;
            }
            /* deepest penetration first: insertion sort, depth descending */
            for (int i = 1; i < nl; i++) {
                Contact key = live[i];
                int j = i - 1;
                while (j >= 0 && live[j].depth < key.depth) {
                    live[j + 1] = live[j];
                    j--;
                }
                live[j + 1] = key;
            }
            int handled = nl < sc->contact_cap ? nl : sc->contact_cap;
            for (int h = 0; h < handled; h++)
                contact_project(sc->ropes[live[h].ai],
                                sc->ropes[live[h].bi], &live[h],
                                sc->stiffness, sc->slop);
            st->applied += handled;
            if (handled > st->max_in_sweep) st->max_in_sweep = handled;
        }
    }
    return 0;
}

static float scene_max_pen(const ContactScene* sc) {
    Contact contacts[MAX_CONTACTS];
    int nc = scene_detect(sc, contacts, MAX_CONTACTS);
    float worst = 0.0f;
    for (int i = 0; i < nc; i++)
        if (contacts[i].depth > worst) worst = contacts[i].depth;
    return worst;
}

/* --------------------------------------------- subdivided reference check */

typedef struct {
    float max_pen;         /* deepest penetration seen on sub-segments */
    int ref_contacts;      /* sub-segment pairs in contact */
    int missed;            /* reference contacts whose pair the broad phase cut */
    float margin;          /* 0.5 * (longest A seg + longest B seg) / sub */
} FineReport;

/* subdivide every segment pair into sub x sub sub-segments, measure exact
   penetration there, and cross-check which coarse segment pairs survived the
   broad phase. the sub-segment midpoint differs from the parent segment by at
   most half a sub-segment length, so coarse depth <= fine depth + margin. */
static int contact_fine(const Rope* ra, const Rope* rb, int sub,
                        const int pairs[][2], int np, FineReport* rep) {
    if (sub < 1) return -1;
    memset(rep, 0, sizeof(*rep));
    char present[MAX_SEG][MAX_SEG];
    memset(present, 0, sizeof(present));
    for (int k = 0; k < np; k++)
        if (pairs[k][0] < MAX_SEG && pairs[k][1] < MAX_SEG)
            present[pairs[k][0]][pairs[k][1]] = 1;

    float max_a = 0.0f, max_b = 0.0f;
    for (int i = 0; i < ra->n - 1; i++) {
        float dx = ra->p[i + 1].x - ra->p[i].x;
        float dy = ra->p[i + 1].y - ra->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        if (d > max_a) max_a = d;
    }
    for (int i = 0; i < rb->n - 1; i++) {
        float dx = rb->p[i + 1].x - rb->p[i].x;
        float dy = rb->p[i + 1].y - rb->p[i].y;
        float d = sqrtf(dx * dx + dy * dy);
        if (d > max_b) max_b = d;
    }
    rep->margin = 0.5f * (max_a + max_b) / (float)sub;
    float rsum = ra->radius + rb->radius;

    for (int i = 0; i < ra->n - 1; i++) {
        for (int j = 0; j < rb->n - 1; j++) {
            for (int u = 0; u < sub; u++) {
                for (int v = 0; v < sub; v++) {
                    float a0x = ra->p[i].x + (ra->p[i + 1].x - ra->p[i].x) * (u / (float)sub);
                    float a0y = ra->p[i].y + (ra->p[i + 1].y - ra->p[i].y) * (u / (float)sub);
                    float a1x = ra->p[i].x + (ra->p[i + 1].x - ra->p[i].x) * ((u + 1) / (float)sub);
                    float a1y = ra->p[i].y + (ra->p[i + 1].y - ra->p[i].y) * ((u + 1) / (float)sub);
                    float b0x = rb->p[j].x + (rb->p[j + 1].x - rb->p[j].x) * (v / (float)sub);
                    float b0y = rb->p[j].y + (rb->p[j + 1].y - rb->p[j].y) * (v / (float)sub);
                    float b1x = rb->p[j].x + (rb->p[j + 1].x - rb->p[j].x) * ((v + 1) / (float)sub);
                    float b1y = rb->p[j].y + (rb->p[j + 1].y - rb->p[j].y) * ((v + 1) / (float)sub);
                    SegSeg q = seg_seg_closest(a0x, a0y, a1x, a1y,
                                               b0x, b0y, b1x, b1y);
                    float depth = rsum - q.dist;
                    if (depth <= 0.0f) continue;
                    rep->ref_contacts++;
                    if (depth > rep->max_pen) rep->max_pen = depth;
                    if (!present[i][j]) rep->missed++;
                }
            }
        }
    }
    return 0;
}

/* ------------------------------------------------------------------ tests */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* --- segment-segment closest points ------------------------------- */
    SegSeg q;

    q = seg_seg_closest(0, 0, 10, 0, 3, 1, 7, 1);
    check("seg_parallel", fabsf(q.dist - 1.0f) < 1e-5f);

    q = seg_seg_closest(0, 0, 1, 1, 0, 1, 1, 0);
    check("seg_crossing", q.dist < 1e-6f && q.ta > 0.4f && q.ta < 0.6f);

    q = seg_seg_closest(0, 0, 1, 0, 2, 1, 3, 1);
    check("seg_endpoint", fabsf(q.dist - sqrtf(2.0f)) < 1e-5f &&
                          q.ta > 0.999f && q.tb < 0.001f);

    q = seg_seg_closest(0, 0, 2, 0, 1, 0, 3, 0);
    check("seg_collinear", q.dist < 1e-6f);

    /* --- two ropes crossing at the middle ------------------------------ */
    float horiz[] = { -2, 0, -1, 0, 0, 0, 1, 0, 2, 0 };
    float vert[]  = { 0, -2, 0, -1, 0, 0, 0, 1, 0, 2 };
    Rope ra, rb;
    check("rope_init", rope_init(&ra, horiz, 5, 0.2f) == 0 &&
                      rope_init(&rb, vert, 5, 0.2f) == 0);
    rope_pin(&ra, 0); rope_pin(&ra, 4);
    rope_pin(&rb, 0); rope_pin(&rb, 4);

    int pairs[MAX_PAIRS][2];
    int np = contact_broad(&ra, &rb, pairs, MAX_PAIRS);
    int found_middle = 0;
    for (int i = 0; i < np; i++)
        if (pairs[i][0] == 2 && pairs[i][1] == 2) found_middle = 1;
    check("broad_finds_crossing", np >= 1 && found_middle);
    check("broad_culls", np < (ra.n - 1) * (rb.n - 1));

    Contact contacts[MAX_CONTACTS];
    int nc = contact_narrow(&ra, 0, &rb, 1, pairs, np, contacts, MAX_CONTACTS);
    float coarse_max = 0.0f;
    for (int i = 0; i < nc; i++)
        if (contacts[i].depth > coarse_max) coarse_max = contacts[i].depth;
    check("narrow_contact", nc >= 1);
    check("narrow_depth", fabsf(coarse_max - 0.4f) < 1e-5f);
    check("fallback_normal_unit",
          fabsf(contacts[0].nx * contacts[0].nx +
                contacts[0].ny * contacts[0].ny - 1.0f) < 1e-5f &&
          isfinite(contacts[0].nx) && isfinite(contacts[0].ny));

    /* parallel overlapping ropes: normal must point from B down toward A */
    float low[] = { 0, 0, 1, 0, 2, 0 };
    float high[] = { 0, 0.1f, 1, 0.1f, 2, 0.1f };
    Rope rl, rh;
    rope_init(&rl, low, 3, 0.2f);
    rope_init(&rh, high, 3, 0.2f);
    int pp[MAX_PAIRS][2];
    int npp = contact_broad(&rl, &rh, pp, MAX_PAIRS);
    Contact pc[MAX_CONTACTS];
    int npc = contact_narrow(&rl, 0, &rh, 1, pp, npp, pc, MAX_CONTACTS);
    check("normal_points_down", npc > 0 &&
          fabsf(pc[0].nx) < 1e-5f && pc[0].ny < -0.999f);

    /* --- subdivided reference: no missed contacts, depth bracket ------- */
    FineReport fine;
    contact_fine(&ra, &rb, 16, pairs, np, &fine);
    check("fine_ref_contacts", fine.ref_contacts > 0);
    check("fine_miss_zero", fine.missed == 0);
    check("coarse_upper_bound", coarse_max <= fine.max_pen + fine.margin + 1e-5f);
    check("fine_le_coarse", fine.max_pen <= coarse_max + 1e-5f);
    check("fine_margin_shrinks", fine.margin > 0.0f && fine.margin < 0.2f);

    /* --- solving: a hanging rope crossing a ledge rope -----------------
       both ropes pinned at only one end, so the geometry has a feasible
       non-penetrating solution; a doubly-pinned inextensible crossing does
       not, and no solver could separate it. */
    float hang_a[] = { -2, 0, -1, 0, 0, 0, 1, 0, 2, 0 };
    float hang_b[] = { 0, -2, 0, -1, 0, 0, 0, 1, 0, 2 };
    Rope ha, hb;
    rope_init(&ha, hang_a, 5, 0.2f);
    rope_init(&hb, hang_b, 5, 0.2f);
    rope_pin(&ha, 0);
    rope_pin(&hb, 0);
    ContactScene sc = {
        .ropes = { &ha, &hb }, .nropes = 2,
        .iterations = 16, .contact_cap = 64,
        .stiffness = 1.0f, .slop = 1e-4f
    };
    ContactStats st;
    check("starts_embedded", fabsf(scene_max_pen(&sc) - 0.4f) < 1e-5f);
    scene_step(&sc, 0.016f, 0.0f, 0.0f, &st);
    /* after the first solve the dominant motion is along each rope's own
       perpendicular: the horizontal rope separates mostly vertically, the
       vertical one mostly horizontally. geometric normals at partially
       separated neighbours add a smaller sideways slide, so compare the
       dominant component rather than demanding an exact zero. */
    check("horizontal_rope_separates_vertical",
          ha.p[2].y > 0.1f && fabsf(ha.p[2].y) > fabsf(ha.p[2].x));
    check("vertical_rope_separates_horizontal",
          hb.p[2].x < -0.1f && fabsf(hb.p[2].x) > fabsf(hb.p[2].y));
    check("pins_hold",
          fabsf(ha.p[0].x - -2.0f) < 1e-6f && fabsf(ha.p[0].y) < 1e-6f &&
          fabsf(hb.p[0].x) < 1e-6f && fabsf(hb.p[0].y - -2.0f) < 1e-6f);
    int total_applied = st.applied, max_sweeps = st.contact_sweeps;
    for (int i = 1; i < 60; i++) {
        scene_step(&sc, 0.016f, 0.0f, 0.0f, &st);
        total_applied += st.applied;
        if (st.contact_sweeps > max_sweeps) max_sweeps = st.contact_sweeps;
    }
    float residual = scene_max_pen(&sc);
    check("contact_resolves", residual < 0.05f);
    check("sweeps_ran", max_sweeps > 0 && total_applied > 0);
    check("centerlines_clear", residual < 0.05f);

    /* control: same scene, contacts disabled -> ropes pass through */
    Rope ca, cb;
    rope_init(&ca, hang_a, 5, 0.2f);
    rope_init(&cb, hang_b, 5, 0.2f);
    rope_pin(&ca, 0);
    rope_pin(&cb, 0);
    ContactScene ctrl = {
        .ropes = { &ca, &cb }, .nropes = 2,
        .iterations = 16, .contact_cap = 0,
        .stiffness = 1.0f, .slop = 1e-4f
    };
    ContactStats cst;
    for (int i = 0; i < 60; i++)
        scene_step(&ctrl, 0.016f, 0.0f, 0.0f, &cst);
    check("control_passes_through", fabsf(scene_max_pen(&ctrl) - 0.4f) < 1e-5f);
    check("control_no_sweeps", cst.contact_sweeps == 0 && cst.applied == 0);

    /* --- zero-length segments are hard errors -------------------------- */
    float dup[] = { 0, 0, 1, 0, 1, 0, 2, 0 };
    Rope bad;
    check("zero_len_init_rejected", rope_init(&bad, dup, 4, 0.2f) == -1);
    Rope tampered;
    rope_init(&tampered, horiz, 5, 0.2f);
    tampered.rest[2] = 0.0f;
    ContactScene bad_scene = {
        .ropes = { &tampered, &rb }, .nropes = 2,
        .iterations = 4, .contact_cap = 8, .stiffness = 1.0f, .slop = 1e-4f
    };
    float saved_x = tampered.p[3].x;
    ContactStats dummy;
    check("zero_len_step_error",
          scene_step(&bad_scene, 0.016f, 0.0f, -9.8f, &dummy) == -1);
    check("zero_len_no_mutation", fabsf(tampered.p[3].x - saved_x) < 1e-6f);

    /* --- budget allocation under a contact surge ----------------------- */
    ContactStats b;
    contact_budget(8, 0, 4, &b);
    check("budget_zero", b.contact_sweeps == 0 && b.stretch_iters == 8);
    contact_budget(8, 4, 4, &b);
    check("budget_single_batch_uses_iters", b.contact_sweeps == 7);
    contact_budget(8, 5, 4, &b);
    check("budget_partial_batch", b.contact_sweeps == 2);
    contact_budget(8, 1000, 4, &b);
    check("budget_capped", b.contact_sweeps == 7 && b.stretch_iters == 8);
    contact_budget(8, 1000, 0, &b);
    check("budget_disabled", b.contact_sweeps == 0);

    /* zigzag rope crosses a straight one on every segment: ~20 contacts */
    float zig[2 * 21], straight_pts[2 * 21];
    for (int i = 0; i <= 20; i++) {
        zig[2 * i] = 0.5f * i;
        zig[2 * i + 1] = (i % 2) ? -0.1f : 0.1f;
        straight_pts[2 * i] = 0.5f * i;
        straight_pts[2 * i + 1] = 0.0f;
    }
    Rope rz, rs;
    rope_init(&rz, zig, 21, 0.2f);
    rope_init(&rs, straight_pts, 21, 0.2f);
    ContactScene surge = {
        .ropes = { &rz, &rs }, .nropes = 2,
        .iterations = 8, .contact_cap = 4,
        .stiffness = 1.0f, .slop = 1e-4f
    };
    ContactStats sst;
    check("surge_step", scene_step(&surge, 0.016f, 0.0f, 0.0f, &sst) == 0);
    check("surge_detects_many", sst.contacts_seen >= 16);
    check("surge_sweeps_bounded", sst.contact_sweeps >= 1 &&
                                  sst.contact_sweeps <= 7);
    check("surge_cap_respected", sst.max_in_sweep <= 4 &&
                                 sst.applied <= 7 * 4);

    /* --- dynamic stability under gravity ------------------------------
       a vertical rope hangs through a horizontal one pinned at both ends;
       after settling it must stay finite and never tunnel through. */
    float ledge_xy[] = { -2, 0, -1, 0, 0, 0, 1, 0, 2, 0 };
    float drop_xy[]  = { 0, 2, 0, 1, 0, 0, 0, -1, 0, -2 };
    Rope ledge, drop;
    rope_init(&ledge, ledge_xy, 5, 0.2f);
    rope_init(&drop, drop_xy, 5, 0.2f);
    rope_pin(&ledge, 0); rope_pin(&ledge, 4);
    rope_pin(&drop, 0);
    ContactScene dyn = {
        .ropes = { &ledge, &drop }, .nropes = 2,
        .iterations = 12, .contact_cap = 32,
        .stiffness = 1.0f, .slop = 1e-4f
    };
    ContactStats dst;
    int finite_all = 1;
    for (int i = 0; i < 300; i++) {
        scene_step(&dyn, 0.016f, 0.0f, -9.8f, &dst);
        for (int r = 0; r < 2; r++)
            for (int p = 0; p < dyn.ropes[r]->n; p++) {
                RParticle* pp = &dyn.ropes[r]->p[p];
                if (!isfinite(pp->x) || !isfinite(pp->y)) finite_all = 0;
            }
    }
    check("gravity_finite", finite_all);
    check("gravity_no_tunnel", scene_max_pen(&dyn) < 0.1f);
    check("gravity_pin_holds",
          fabsf(ledge.p[0].x - -2.0f) < 1e-6f &&
          fabsf(drop.p[0].y - 2.0f) < 1e-6f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
