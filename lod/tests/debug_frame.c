/* diagnostic harness: include the module directly, reproduce one frame */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/lod.c"

int main(int argc, char** argv) {
    int target = argc > 1 ? atoi(argv[1]) : 20;
    Skeleton sk; ChainTable ct;
    skel_build(&sk); chains_build(&ct);
    Sim sim;
    float dt = 1.0f / 30.0f;
    for (int frame = 0; frame <= target; frame++) {
        float t = frame * dt;
        float dist = sim_distance(t);
        Vec2 cam = v2(sim_camx(t), 0.0f);
        Frustum fr = frustum_make(cam, dist);
        skel_pose(&sk, t); skel_fk(&sk);
        int before = sim.lod;
        int switched = lod_step(&sim.lod, dist);
        if (frame != target) continue;

        RenderList rl, rt;
        build_render(&sk, &ct, sim.lod, 1, &rl);
        build_render(&sk, &ct, 0, 1, &rt);
        int vl[MAX_RENDER], vt[MAX_RENDER];
        cull_list(&rl, &fr, 1, vl);
        cull_list(&rt, &fr, 1, vt);

        printf("frame=%d t=%.3f dist=%.3f lod=%d->%d sw=%d\n",
               frame, t, dist, before, sim.lod, switched);
        printf("--- merged (linked) ---\n");
        for (int i = 0; i < rl.n; i++)
            printf("[%d] %s vis=%d p0=(%.3f,%.3f) p1=(%.3f,%.3f) r=%.4f "
                   "tipj=%d cover=%d\n", i,
                   rl.b[i].chain >= 0 ? ct.c[rl.b[i].chain].name : "bone",
                   vl[i], rl.b[i].p0.x, rl.b[i].p0.y,
                   rl.b[i].p1.x, rl.b[i].p1.y, rl.b[i].radius,
                   rl.b[i].tip_joint, rl.b[i].cover);
        printf("--- LOD0 truth bones (visible only) ---\n");
        for (int i = 0; i < rt.n; i++)
            if (vt[i])
                printf("[%d] p0=(%.3f,%.3f) p1=(%.3f,%.3f) r=%.4f tipj=%d\n",
                       i, rt.b[i].p0.x, rt.b[i].p0.y,
                       rt.b[i].p1.x, rt.b[i].p1.y, rt.b[i].radius,
                       rt.b[i].tip_joint);

        Framebuf fl, ft;
        fb_init(&fl, SCR_W, SCR_H); fb_init(&ft, SCR_W, SCR_H);
        render_scene(&fl, &rl, vl, &fr, sim.lod, 1, 0);
        render_scene(&ft, &rt, vt, &fr, 0, 2, 0);
        int n = 0;
        for (int y = 0; y < SCR_H; y++)
            for (int x = 0; x < SCR_W; x++) {
                int idx = (y * SCR_W + x) * 3;
                int tc = !(ft.px[idx] == 24 && ft.px[idx+1] == 24 &&
                           ft.px[idx+2] == 28);
                int lc = !(fl.px[idx] == 24 && fl.px[idx+1] == 24 &&
                           fl.px[idx+2] == 28);
                if (tc && !lc) {
                    if (n < 30) printf("missing px=(%d,%d)\n", x, y);
                    n++;
                }
            }
        printf("missing_total=%d\n", n);

        Camera cam0;
        cam0.center = fr.center;
        cam0.scale = FOCAL / fr.dist;
        /* probe the first missing pixel against every merged capsule */
        int sx = 135, sy = 120;
        double pxx = sx + 0.5, pyy = sy + 0.5;
        for (int i = 0; i < rl.n; i++) {
            double axx = SCR_W * 0.5 + cam0.scale *
                         (rl.b[i].p0.x - cam0.center.x);
            double ayy = SCR_H * 0.5 - cam0.scale *
                         (rl.b[i].p0.y - cam0.center.y);
            double bxx = SCR_W * 0.5 + cam0.scale *
                         (rl.b[i].p1.x - cam0.center.x);
            double byy = SCR_H * 0.5 - cam0.scale *
                         (rl.b[i].p1.y - cam0.center.y);
            double ddx = bxx - axx, ddy = byy - ayy;
            double ll2 = ddx * ddx + ddy * ddy;
            double tt = ((pxx - axx) * ddx + (pyy - ayy) * ddy) / ll2;
            if (tt < 0) tt = 0;
            if (tt > 1) tt = 1;
            double qxx = axx + tt * ddx, qyy = ayy + tt * ddy;
            double dd2 = (pxx - qxx) * (pxx - qxx) +
                         (pyy - qyy) * (pyy - qyy);
            double rr = rl.b[i].cover ? ceil(rl.b[i].radius * cam0.scale)
                                      : floor(rl.b[i].radius * cam0.scale);
            printf("probe cap[%d] t=%.3f d2=%.3f r=%.3f r2=%.3f fill=%d\n",
                   i, tt, dd2, rr, rr * rr, dd2 <= rr * rr);
        }
        /* what color is the linked buffer there? */
        int idx = (sy * SCR_W + sx) * 3;
        printf("linked px rgb=(%d,%d,%d) truth rgb=(%d,%d,%d)\n",
               fl.px[idx], fl.px[idx+1], fl.px[idx+2],
               ft.px[idx], ft.px[idx+1], ft.px[idx+2]);

        /* same probe against LOD0 truth bones: which one fills it? */
        for (int i = 0; i < rt.n; i++) {
            double axx = SCR_W * 0.5 + cam0.scale *
                         (rt.b[i].p0.x - cam0.center.x);
            double ayy = SCR_H * 0.5 - cam0.scale *
                         (rt.b[i].p0.y - cam0.center.y);
            double bxx = SCR_W * 0.5 + cam0.scale *
                         (rt.b[i].p1.x - cam0.center.x);
            double byy = SCR_H * 0.5 - cam0.scale *
                         (rt.b[i].p1.y - cam0.center.y);
            double ddx = bxx - axx, ddy = byy - ayy;
            double ll2 = ddx * ddx + ddy * ddy;
            double tt = ((pxx - axx) * ddx + (pyy - ayy) * ddy) / ll2;
            if (tt < 0) tt = 0;
            if (tt > 1) tt = 1;
            double qxx = axx + tt * ddx, qyy = ayy + tt * ddy;
            double dd2 = (pxx - qxx) * (pxx - qxx) +
                         (pyy - qyy) * (pyy - qyy);
            double rr = floor(rt.b[i].radius * cam0.scale);
            if (dd2 <= rr * rr)
                printf("truth bone[%d] tipj=%d fills it t=%.3f d2=%.2f r=%.1f\n",
                       i, rt.b[i].tip_joint, tt, dd2, rr);
        }
        /* merged leg capsules full geometry */
        for (int i = 4; i < rl.n; i++) {
            double axx = SCR_W * 0.5 + cam0.scale *
                         (rl.b[i].p0.x - cam0.center.x);
            double ayy = SCR_H * 0.5 - cam0.scale *
                         (rl.b[i].p0.y - cam0.center.y);
            double bxx = SCR_W * 0.5 + cam0.scale *
                         (rl.b[i].p1.x - cam0.center.x);
            double byy = SCR_H * 0.5 - cam0.scale *
                         (rl.b[i].p1.y - cam0.center.y);
            printf("merged[%d] screen (%.1f,%.1f)->(%.1f,%.1f) r=%.4f\n",
                   i, axx, ayy, bxx, byy, rl.b[i].radius);
        }
        fb_free(&fl); fb_free(&ft);
    }
    return 0;
}
