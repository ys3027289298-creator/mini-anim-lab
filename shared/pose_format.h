/*
 * pose_format.h - one shared pose export format for skeleton/ik/curve/path.
 *
 * The four modules used to write slightly different layouts.  They now all
 * emit the same line-oriented container:
 *
 *   mal-pose 1
 *   kind <skeleton|ik|curve|path>
 *   count <record count>
 *   <record-name> <key>=<value> <key>=<value> ...
 *
 * The container (magic, version, kind, count, record grammar) is identical
 * everywhere; only the record name and field set differ per module:
 *
 *   skeleton : bone index= parent= lx= ly= lang= wx= wy= ang=
 *   ik       : ik   root_x= root_y= target_x= target_y= ok=
 *                   joint_x= joint_y= a1= a2=
 *   curve    : key  t= v= tin= tout=
 *   path     : meta total=          (one record, before the points)
 *              pt   x= y= cum=
 *
 * Floats are written with %.9g.  That is FLT_DECIMAL_DIG for binary32, so
 * strtod() followed by a cast back to float reproduces the exact original
 * bit pattern; readers can therefore compare with == and get a true
 * roundtrip rather than a tolerance match.
 *
 * Header-only and static-inline so each module directory keeps building with
 * its own single-file Makefile and no extra flags.
 */
#ifndef MAL_POSE_FORMAT_H
#define MAL_POSE_FORMAT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAL_POSE_MAGIC "mal-pose"
#define MAL_POSE_VERSION 1

#define MAL_POSE_LINE 512
#define MAL_POSE_MAX_FIELDS 16
#define MAL_POSE_KEY_MAX 32
#define MAL_POSE_NAME_MAX 32

typedef struct {
    char key[MAL_POSE_KEY_MAX];
    double val;
} PfField;

typedef struct {
    char name[MAL_POSE_NAME_MAX];
    PfField field[MAL_POSE_MAX_FIELDS];
    int n;
} PfRecord;

/* --------------------------------------------------------------- writing */

static inline void pf_write_header(FILE* f, const char* kind, int count) {
    fprintf(f, "%s %d\nkind %s\ncount %d\n",
            MAL_POSE_MAGIC, MAL_POSE_VERSION, kind, count);
}

static inline void pf_record_begin(FILE* f, const char* name) {
    fputs(name, f);
}

static inline void pf_field_int(FILE* f, const char* key, int v) {
    fprintf(f, " %s=%d", key, v);
}

static inline void pf_field_float(FILE* f, const char* key, float v) {
    /* 9 significant digits roundtrip a binary32 exactly. */
    fprintf(f, " %s=%.9g", key, (double)v);
}

static inline void pf_record_end(FILE* f) {
    fputc('\n', f);
}

/* --------------------------------------------------------------- reading */

/* Consume the rest of the current line (normally just '\n'). */
static inline void pf_skip_eol(FILE* f) {
    int c;
    do { c = fgetc(f); } while (c != '\n' && c != EOF);
}

/*
 * Read and validate "mal-pose <version>" / "kind <want_kind>" / "count <n>".
 * Returns 0 and stores n in *count on success, -1 otherwise.
 */
static inline int pf_read_header(FILE* f, const char* want_kind, int* count) {
    char a[MAL_POSE_KEY_MAX], b[MAL_POSE_KEY_MAX];
    int ver = 0, n = 0;

    if (fscanf(f, "%31s %d", a, &ver) != 2) return -1;
    if (strcmp(a, MAL_POSE_MAGIC) != 0 || ver != MAL_POSE_VERSION) return -1;

    if (fscanf(f, "%31s %31s", a, b) != 2) return -1;
    if (strcmp(a, "kind") != 0 || strcmp(b, want_kind) != 0) return -1;

    if (fscanf(f, "%31s %d", a, &n) != 2) return -1;
    if (strcmp(a, "count") != 0 || n < 0) return -1;

    pf_skip_eol(f);
    *count = n;
    return 0;
}

/*
 * Read one record line into r.  Blank lines and '#' comments are skipped.
 * Returns 0 on a record, -1 on EOF or malformed input.
 */
static inline int pf_read_record(FILE* f, PfRecord* r) {
    char line[MAL_POSE_LINE];

    while (fgets(line, sizeof line, f)) {
        char* p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (*p == '\0' || *p == '#') continue;

        char* tok = strtok(p, " \t\r\n");
        if (!tok) continue;
        snprintf(r->name, sizeof r->name, "%s", tok);

        r->n = 0;
        while ((tok = strtok(NULL, " \t\r\n")) != NULL) {
            if (r->n >= MAL_POSE_MAX_FIELDS) return -1;
            char* eq = strchr(tok, '=');
            if (!eq || eq == tok || eq[1] == '\0') return -1;
            *eq = '\0';
            snprintf(r->field[r->n].key, sizeof r->field[r->n].key, "%s", tok);

            char* end = NULL;
            double v = strtod(eq + 1, &end);
            if (end == eq + 1 || *end != '\0') return -1;
            r->field[r->n].val = v;
            r->n++;
        }
        return 0;
    }
    return -1;
}

/* Look up a field by key; returns 0 and stores the value, -1 if absent. */
static inline int pf_get(const PfRecord* r, const char* key, double* out) {
    for (int i = 0; i < r->n; i++) {
        if (strcmp(r->field[i].key, key) == 0) {
            *out = r->field[i].val;
            return 0;
        }
    }
    return -1;
}

/* Same, but for float-typed fields. */
static inline int pf_get_float(const PfRecord* r, const char* key, float* out) {
    double d;
    if (pf_get(r, key, &d) != 0) return -1;
    *out = (float)d;
    return 0;
}

/* Same, but for int-typed fields. */
static inline int pf_get_int(const PfRecord* r, const char* key, int* out) {
    double d;
    if (pf_get(r, key, &d) != 0) return -1;
    *out = (int)d;
    return 0;
}

#endif /* MAL_POSE_FORMAT_H */
