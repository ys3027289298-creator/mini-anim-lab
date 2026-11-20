
/* timeline: unified fixed-step timeline merging fixed/audio/network clocks.
 *
 * - master clock is an integer tick count: time = tick * dt, so float error
 *   never accumulates; an accumulator emits fixed steps from elapsed time.
 * - network clock offset/drift is estimated from round-trip samples with
 *   weighted least squares (weight 1/rtt): theta(t) = a + b*t.
 * - events from all three sources quantize onto the tick grid and come out
 *   sorted by (tick, priority, arrival seq); the seq key makes it stable.
 * - drift of an external reference clock is compensated one bounded step at
 *   a time; any clock moving backwards returns TL_ERR_BACKWARDS.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>

#define TL_EPS 1e-9
#define TL_MAX_EVENTS 128
#define TL_MAX_RTT_SAMPLES 16

enum { TL_OK = 0, TL_ERR_PARAM = -1, TL_ERR_FULL = -2, TL_ERR_BACKWARDS = -3 };

typedef enum { SRC_FIXED = 0, SRC_AUDIO = 1, SRC_NETWORK = 2, SRC_COUNT = 3 } TlSource;

typedef struct {
    double t1, t2, t3, t4; /* local send, remote recv, remote send, local recv */
} RttSample;

typedef struct {
    RttSample ring[TL_MAX_RTT_SAMPLES];
    int count;          /* valid samples, capped at TL_MAX_RTT_SAMPLES */
    int head;           /* next ring slot to overwrite */
    double last_remote; /* last t3, for remote monotonicity */
    int have_remote;
    double a, b;        /* fitted theta(t) = a + b*t */
    double min_rtt;
    int fitted;
} NetClock;

typedef struct {
    int id;
    TlSource src;
    int priority;       /* smaller value fires first */
    uint64_t seq;       /* global arrival order, stable tiebreak */
    int64_t tick;
    double local_time;  /* timestamp in the local clock domain */
} TlEvent;

typedef struct {
    double dt;
    int64_t tick;       /* grid clock: current time is tick * dt */
    double acc;         /* residual time not yet turned into ticks */
    NetClock net;
    TlEvent events[TL_MAX_EVENTS];
    int pending;
    uint64_t seq_next;
    double last_time[SRC_COUNT];
    int seen[SRC_COUNT];
    int comp_up;        /* extra ticks emitted by drift compensation */
    int comp_down;      /* delayed ticks */
} Timeline;

/* ---- network clock: round-trip offset fit ---- */

static void net_init(NetClock* n) { memset(n, 0, sizeof(*n)); }

/* Weighted least squares over the ring, weight w = 1/(rtt + eps), fitting
 * theta(t) = a + b*t.  rtt and theta follow the NTP four-timestamp form:
 *   rtt   = (t4-t1) - (t3-t2)
 *   theta = ((t2-t1) + (t3-t4)) / 2                       */
static void net_fit(NetClock* n) {
    double sw = 0, st = 0, sth = 0, stt = 0, stth = 0;
    n->min_rtt = 1e300;
    for (int i = 0; i < n->count; i++) {
        const RttSample* s = &n->ring[i];
        double rtt = (s->t4 - s->t1) - (s->t3 - s->t2);
        double theta = ((s->t2 - s->t1) + (s->t3 - s->t4)) * 0.5;
        if (rtt < n->min_rtt) n->min_rtt = rtt;
        double w = 1.0 / (rtt + TL_EPS);
        sw += w;
        st += w * s->t1;
        sth += w * theta;
        stt += w * s->t1 * s->t1;
        stth += w * s->t1 * theta;
    }
    double denom = sw * stt - st * st;
    if (fabs(denom) < TL_EPS) {
        n->b = 0.0;
        n->a = sth / sw;
    } else {
        /* WLS normal equations for theta = a + b*t:
         *   b = (Sw*S(tθ) - St*Sθ) / (Sw*S(t²) - St²),  a = (Sθ - b*St)/Sw */
        n->b = (sw * stth - st * sth) / denom;
        n->a = (sth - n->b * st) / sw;
    }
    n->fitted = 1;
}

/* Add one round-trip sample.  Rejects impossible timing and a remote clock
 * that moves backwards between samples. */
static int net_add_sample(NetClock* n, double t1, double t2, double t3, double t4) {
    if (!(t4 >= t1 - TL_EPS) || !(t3 >= t2 - TL_EPS) || t4 < 0.0 || t1 < 0.0)
        return TL_ERR_PARAM;
    double rtt = (t4 - t1) - (t3 - t2);
    if (rtt < -TL_EPS) return TL_ERR_PARAM;
    if (n->have_remote && t3 < n->last_remote - TL_EPS) return TL_ERR_BACKWARDS;
    RttSample* s = &n->ring[n->head];
    s->t1 = t1; s->t2 = t2; s->t3 = t3; s->t4 = t4;
    n->head = (n->head + 1) % TL_MAX_RTT_SAMPLES;
    if (n->count < TL_MAX_RTT_SAMPLES) n->count++;
    n->last_remote = t3;
    n->have_remote = 1;
    net_fit(n);
    return TL_OK;
}

/* local(t_remote) = t_remote - theta(t_local_ref).  Needs one sample. */
static double net_to_local(const NetClock* n, double remote_time, double local_ref) {
    if (!n->fitted) return remote_time;
    return remote_time - (n->a + n->b * local_ref);
}

/* ---- timeline ---- */

static void timeline_init(Timeline* tl, double dt) {
    memset(tl, 0, sizeof(*tl));
    tl->dt = dt;
    net_init(&tl->net);
}

static double timeline_time(const Timeline* tl) { return (double)tl->tick * tl->dt; }

static int64_t quantize_tick(double t, double dt) {
    return (int64_t)floor(t / dt + 0.5);
}

/* Merge one event.  remote_time is used for SRC_NETWORK and converted through
 * the fitted clock; local_time is the local receive instant for that event.
 * Every source's local-domain timestamp must be non-decreasing. */
static int timeline_push(Timeline* tl, TlSource src, int id, int priority,
                         double local_time, double remote_time) {
    if (tl == NULL || src < 0 || src >= SRC_COUNT || local_time < -TL_EPS)
        return TL_ERR_PARAM;
    if (local_time < 0.0) local_time = 0.0;

    double local_stamp = local_time;
    if (src == SRC_NETWORK) {
        if (remote_time < -TL_EPS) return TL_ERR_PARAM;
        local_stamp = net_to_local(&tl->net, remote_time, local_time);
        if (local_stamp < -TL_EPS) return TL_ERR_PARAM;
        if (local_stamp < 0.0) local_stamp = 0.0;
    }

    if (tl->seen[src] && local_stamp < tl->last_time[src] - TL_EPS)
        return TL_ERR_BACKWARDS;
    if (tl->pending >= TL_MAX_EVENTS) return TL_ERR_FULL;

    TlEvent* e = &tl->events[tl->pending++];
    e->id = id;
    e->src = src;
    e->priority = priority;
    e->seq = tl->seq_next++;
    e->local_time = local_stamp;
    e->tick = quantize_tick(local_stamp, tl->dt);
    tl->last_time[src] = local_stamp;
    tl->seen[src] = 1;
    return TL_OK;
}

/* Consume wall time and emit fixed steps.  Negative elapsed is a backwards
 * clock.  Returns the new tick count, or a negative error. */
static int64_t timeline_advance(Timeline* tl, double elapsed) {
    if (elapsed < -TL_EPS) return TL_ERR_BACKWARDS;
    if (elapsed < 0.0) elapsed = 0.0;
    tl->acc += elapsed;
    int64_t steps = (int64_t)floor(tl->acc / tl->dt + TL_EPS);
    if (steps > 0) {
        tl->tick += steps;
        tl->acc -= (double)steps * tl->dt;
    }
    return tl->tick;
}

/* Drift compensation against an external reference (e.g. audio hardware
 * clock).  While the reference and the grid differ by a full step or more,
 * nudge the accumulator one bounded step at a time, so the master converges
 * without ever jumping; on return |drift| is guaranteed below one dt. */
static int timeline_correct(Timeline* tl, double ref_time) {
    if (ref_time < -TL_EPS) return TL_ERR_PARAM;
    double drift = ref_time - timeline_time(tl);
    while (drift >= tl->dt) {                 /* reference ahead: emit a step */
        tl->acc += tl->dt;
        tl->comp_up++;
        drift -= tl->dt;
    }
    while (drift <= -tl->dt) {                /* reference behind: delay a step */
        if (tl->acc >= tl->dt) {
            tl->acc -= tl->dt;
        } else if (tl->tick > 0) {
            tl->tick--;
            tl->acc += tl->dt;
        } else {
            break;                            /* nothing to delay at t = 0 */
        }
        tl->comp_down++;
        drift += tl->dt;
    }
    return TL_OK;
}

/* ---- stable output ordering ---- */

static int event_less(const TlEvent* a, const TlEvent* b) {
    if (a->tick != b->tick) return a->tick < b->tick;
    if (a->priority != b->priority) return a->priority < b->priority;
    return a->seq < b->seq;
}

static void sort_events(TlEvent* ev, int n) { /* insertion sort, stable by key */
    for (int i = 1; i < n; i++) {
        TlEvent key = ev[i];
        int j = i - 1;
        while (j >= 0 && event_less(&key, &ev[j])) {
            ev[j + 1] = ev[j];
            j--;
        }
        ev[j + 1] = key;
    }
}

/* Drain events whose tick is due into out (capacity max_out), in stable order. */
static int timeline_drain(Timeline* tl, TlEvent* out, int max_out) {
    sort_events(tl->events, tl->pending);
    int n = 0;
    int i = 0;
    while (i < tl->pending && tl->events[i].tick <= tl->tick) {
        if (n < max_out) out[n++] = tl->events[i];
        i++;
    }
    int remain = tl->pending - i;
    memmove(tl->events, tl->events + i, (size_t)remain * sizeof(TlEvent));
    tl->pending = remain;
    return n;
}

static int timeline_pending(const Timeline* tl) { return tl->pending; }

/* ============================== self test ============================== */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ground-truth remote clock: theta(t) = theta0 + (r-1)*t.
 * one_way is the symmetric one-way delay; asym adds bias to the return path. */
static void make_sample(NetClock* n, double t1, double one_way,
                        double theta0, double r, double asym) {
    double t2 = (t1 + one_way) * r + theta0;                 /* remote recv */
    double t3 = t2;                                          /* instant reply */
    double theta_at_t4 = theta0 + (r - 1.0) * (t1 + 2.0 * one_way + asym);
    double t4 = t3 - theta_at_t4 + one_way + asym;           /* local recv */
    (void)n;
    net_add_sample(n, t1, t2, t3, t4);
}

int main(void) {
    /* ---- 1. network offset/drift fit ---- */
    {
        NetClock n;
        net_init(&n);
        const double theta0 = 0.07, r = 1.0005;
        double sym[][2] = {
            {0.0, 0.0005}, {1.0, 0.0005}, {2.0, 0.0005}, {3.0, 0.0005},
            {4.0, 0.0005}, {5.0, 0.0005}, {6.0, 0.0005}, {7.0, 0.0005},
        };
        for (size_t i = 0; i < sizeof(sym) / sizeof(sym[0]); i++)
            make_sample(&n, sym[i][0], sym[i][1], theta0, r, 0.0);
        /* two high-rtt asymmetric outliers: their theta is biased but the
         * 1/rtt weight must keep the fit close to the truth. */
        make_sample(&n, 8.0, 0.02, theta0, r, 0.02);
        make_sample(&n, 9.0, 0.02, theta0, r, 0.02);

        check("fit_a", fabs(n.a - theta0) < 1e-3);
        check("fit_b", fabs(n.b - (r - 1.0)) < 1e-4);
        check("min_rtt", n.min_rtt > 0.0 && n.min_rtt < 0.002);

        /* round trip: remote -> local -> remote lands where it started */
        double rem = 10.0 * r + theta0;
        double loc = net_to_local(&n, rem, 10.0);
        check("net_to_local", fabs(loc - 10.0) < 2e-3);

        /* ring overwrite keeps exactly the capacity */
        NetClock cap;
        net_init(&cap);
        for (int i = 0; i < TL_MAX_RTT_SAMPLES + 5; i++)
            net_add_sample(&cap, (double)i, (double)i + 0.1,
                           (double)i + 0.2, (double)i + 0.3);
        check("rtt_capacity", cap.count == TL_MAX_RTT_SAMPLES);
    }

    /* ---- 2. backwards clocks are errors ---- */
    {
        NetClock n;
        net_init(&n);
        check("rtt_bad_order", net_add_sample(&n, 1.0, 1.1, 1.2, 0.9) == TL_ERR_PARAM);
        net_add_sample(&n, 0.0, 1.0, 2.0, 2.1);      /* remote t3 = 2.1 */
        /* internally valid (t3 >= t2), but remote t3 1.8 is behind 2.1 */
        check("rtt_backwards", net_add_sample(&n, 1.0, 1.5, 1.8, 2.6)
                               == TL_ERR_BACKWARDS);

        Timeline tl;
        timeline_init(&tl, 0.1);
        check("push_ok", timeline_push(&tl, SRC_AUDIO, 1, 0, 0.5, 0.0) == TL_OK);
        check("push_backwards", timeline_push(&tl, SRC_AUDIO, 2, 0, 0.4, 0.0)
                               == TL_ERR_BACKWARDS);
        check("advance_backwards", timeline_advance(&tl, -0.1) == TL_ERR_BACKWARDS);
        check("bad_src", timeline_push(&tl, (TlSource)9, 1, 0, 0.1, 0.0)
                        == TL_ERR_PARAM);
    }

    /* ---- 3. fixed-step accumulator has no cumulative float error ---- */
    {
        Timeline tl;
        timeline_init(&tl, 0.1);
        int64_t t = 0;
        for (int i = 0; i < 100000; i++) t = timeline_advance(&tl, 0.1);
        check("exact_steps", t == 100000);
        check("zero_residual", fabs(tl.acc) < TL_EPS);
        check("exact_time", fabs(timeline_time(&tl) - 10000.0) < TL_EPS);

        Timeline tl2;
        timeline_init(&tl2, 1.0 / 3.0);
        for (int i = 0; i < 30; i++) timeline_advance(&tl2, 0.1);
        check("fractional_steps", tl2.tick == 9 && fabs(tl2.acc) < TL_EPS + 1e-12);
    }

    /* ---- 4. drift compensation stays bounded ---- */
    {
        Timeline tl;
        timeline_init(&tl, 0.1);
        double lag = 0.0;
        for (int frame = 1; frame <= 200; frame++) {
            timeline_advance(&tl, 0.1);
            lag += 0.102; /* reference runs 2% fast */
            timeline_correct(&tl, lag);
        }
        double residual = lag - timeline_time(&tl);
        check("comp_bounded", fabs(residual) < tl.dt + TL_EPS);
        check("comp_used", tl.comp_up > 0 && tl.comp_down == 0);

        Timeline slow;
        timeline_init(&slow, 0.1);
        double lead = 0.0;
        for (int frame = 1; frame <= 200; frame++) {
            timeline_advance(&slow, 0.1);
            lead += 0.098; /* reference runs 2% slow */
            timeline_correct(&slow, lead);
        }
        check("comp_slow_bounded",
              fabs(lead - timeline_time(&slow)) < slow.dt + TL_EPS);
        check("comp_slow_used", slow.comp_down > 0);
    }

    /* ---- 5. three-source merge: order and counts match a plain sequence ---- */
    {
        NetClock n;
        net_init(&n);
        const double theta0 = 0.05, r = 1.0;
        for (int i = 0; i < 8; i++)
            make_sample(&n, (double)i, 0.0005, theta0, r, 0.0);

        Timeline tl;
        timeline_init(&tl, 0.1);
        tl.net = n;

        /* reference: events in raw local time, built per source, then sorted */
        typedef struct {
            double t; int prio; uint64_t seq; int id; TlSource src;
        } RefEvent;
        RefEvent ref[TL_MAX_EVENTS];
        int rn = 0;
        int expect[SRC_COUNT] = {0, 0, 0};

        /* fixed: grid-aligned; audio: offset between grids; network: remote
         * timestamps converted back into the local domain by the fit. */
        double fixed_t[] = {0.0, 0.3, 0.7, 1.2};
        double audio_t[] = {0.15, 0.45, 0.95};
        double net_remote[] = {0.25 + theta0, 0.85 + theta0, 1.05 + theta0};
        double net_recv[]   = {0.26, 0.86, 1.06};

        for (int i = 0; i < 4; i++) {
            timeline_push(&tl, SRC_FIXED, 100 + i, 1, fixed_t[i], 0.0);
            ref[rn++] = (RefEvent){fixed_t[i], 1, 0, 100 + i, SRC_FIXED};
            expect[SRC_FIXED]++;
        }
        for (int i = 0; i < 3; i++) {
            timeline_push(&tl, SRC_AUDIO, 200 + i, 2, audio_t[i], 0.0);
            ref[rn++] = (RefEvent){audio_t[i], 2, 0, 200 + i, SRC_AUDIO};
            expect[SRC_AUDIO]++;
        }
        for (int i = 0; i < 3; i++) {
            timeline_push(&tl, SRC_NETWORK, 300 + i, 3, net_recv[i], net_remote[i]);
            double loc = net_to_local(&tl.net, net_remote[i], net_recv[i]);
            ref[rn++] = (RefEvent){loc, 3, 0, 300 + i, SRC_NETWORK};
            expect[SRC_NETWORK]++;
        }
        check("push_count", tl.pending == rn && rn == 10);

        /* reference ordering: quantize the same way, then tick/prio, raw order
         * as the final tiebreak (matches arrival seq within one push run). */
        for (int i = 0; i < rn; i++)
            for (int j = i + 1; j < rn; j++) {
                int64_t ti = quantize_tick(ref[i].t, 0.1);
                int64_t tj = quantize_tick(ref[j].t, 0.1);
                if (tj < ti || (tj == ti && ref[j].prio < ref[i].prio) ||
                    (tj == ti && ref[j].prio == ref[i].prio && ref[j].t < ref[i].t)) {
                    RefEvent tmp = ref[i]; ref[i] = ref[j]; ref[j] = tmp;
                }
            }

        TlEvent out[TL_MAX_EVENTS];
        tl.tick = INT64_MAX / 4; /* everything due */
        int got = timeline_drain(&tl, out, TL_MAX_EVENTS);
        check("drain_count", got == rn);
        check("drain_empty", timeline_pending(&tl) == 0);

        int order_ok = 1, id_ok = 1, got_src[SRC_COUNT] = {0, 0, 0};
        for (int i = 0; i < got; i++) {
            got_src[out[i].src]++;
            if (out[i].id != ref[i].id) id_ok = 0;
            if (i > 0) {
                const TlEvent* a = &out[i - 1];
                const TlEvent* b = &out[i];
                if (event_less(b, a)) order_ok = 0;
            }
        }
        check("merge_order", order_ok);
        check("merge_matches_ref", id_ok);
        check("merge_src_counts",
              got_src[0] == expect[0] && got_src[1] == expect[1] &&
              got_src[2] == expect[2]);
        check("net_quantized", out[0].tick >= 0);
    }

    /* ---- 6. same-tick priority ordering is stable ---- */
    {
        Timeline tl;
        timeline_init(&tl, 0.1);
        /* all land on tick 5; pushed out of priority order */
        timeline_push(&tl, SRC_FIXED, 1, 5, 0.5, 0.0);
        timeline_push(&tl, SRC_AUDIO, 2, 1, 0.5, 0.0);
        timeline_push(&tl, SRC_NETWORK, 3, 3, 0.51, 0.5);
        TlEvent out[8];
        tl.tick = 10;
        int n = timeline_drain(&tl, out, 8);
        check("same_tick_count", n == 3);
        check("same_tick_prio", out[0].id == 2 && out[1].id == 3 && out[2].id == 1);

        /* equal priority keeps arrival order */
        Timeline tl2;
        timeline_init(&tl2, 0.1);
        timeline_push(&tl2, SRC_FIXED, 7, 0, 0.2, 0.0);
        timeline_push(&tl2, SRC_FIXED, 8, 0, 0.2, 0.0);
        timeline_push(&tl2, SRC_FIXED, 9, 0, 0.2, 0.0);
        tl2.tick = 5;
        TlEvent o2[8];
        timeline_drain(&tl2, o2, 8);
        check("stable_arrival", o2[0].id == 7 && o2[1].id == 8 && o2[2].id == 9);
    }

    /* ---- 7. scaling vs no-scaling: same order and counts ---- */
    {
        Timeline scaled, plain;
        timeline_init(&scaled, 0.1);
        timeline_init(&plain, 0.1);
        double times[] = {0.05, 0.12, 0.31, 0.34, 0.7, 0.99, 1.11};
        int prios[] = {2, 1, 1, 3, 2, 1, 2};
        for (int i = 0; i < 7; i++) {
            timeline_push(&scaled, SRC_FIXED, i, prios[i], times[i], 0.0);
            /* "no scaling": same events, same quantized grid, same keys */
            timeline_push(&plain, SRC_FIXED, i, prios[i], times[i], 0.0);
        }
        TlEvent a[16], b[16];
        scaled.tick = 100;
        plain.tick = 100;
        int na = timeline_drain(&scaled, a, 16);
        int nb = timeline_drain(&plain, b, 16);
        int same = na == nb;
        for (int i = 0; same && i < na; i++)
            same = a[i].id == b[i].id && a[i].tick == b[i].tick;
        check("scale_equiv", same);

        /* advancing under two different feeding schedules converges to the
         * same ticks: fractional accumulation never reorders events. */
        Timeline s1, s2;
        timeline_init(&s1, 0.1);
        timeline_init(&s2, 0.1);
        for (int i = 0; i < 10; i++) timeline_advance(&s1, 0.1);
        for (int i = 0; i < 20; i++) timeline_advance(&s2, 0.05);
        check("schedule_converge", s1.tick == s2.tick && s1.tick == 10);
    }

    /* ---- 8. capacity ---- */
    {
        Timeline tl;
        timeline_init(&tl, 0.1);
        int rc = TL_OK;
        for (int i = 0; i < TL_MAX_EVENTS + 4; i++)
            rc = timeline_push(&tl, SRC_FIXED, i, 0, (double)i * 0.1, 0.0);
        check("event_full", rc == TL_ERR_FULL && tl.pending == TL_MAX_EVENTS);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
