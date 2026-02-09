/*
 * main_gui.c — HalalSeq desktop GUI application.
 *
 * SDL2 + Nuklear immediate-mode GUI.  Single window (960x640), left
 * panel for file selection / progress, right panel for results or
 * database viewer.
 *
 * Build:  see Makefile "gui" target or CMakeLists.txt.
 */

/* ================================================================== */
/* Nuklear implementation (must come first, once per translation unit) */
/* ================================================================== */
#include "nuklear_setup.h"

/* ================================================================== */
/* Project headers                                                     */
/* ================================================================== */
#include "gui_analysis.h"
#include "tinyfiledialogs.h"
#include "refdb.h"
#include "index.h"
#include "kmer.h"

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <stdlib.h>

/* ================================================================== */
/* Constants                                                           */
/* ================================================================== */
#define WINDOW_W  960
#define WINDOW_H  640
#define LEFT_W    340
#define VERSION   "0.2.0"

/* Row heights (logical points) */
#define ROW_TITLE  30
#define ROW_LABEL  26
#define ROW_BTN    38
#define ROW_SMALL  22
#define ROW_TABLE  24
#define ROW_SPACE  8
#define ROW_BAR    20

/* ================================================================== */
/* GUI state bundle                                                    */
/* ================================================================== */
typedef struct {
    analysis_context_t analysis;
    char index_path[1024];
    int  right_panel_mode;    /* 0 = Results, 1 = Database */
    halal_refdb_t *db_info;   /* loaded at startup for database viewer */
    float min_display_pct;    /* minimum % to show in chart, default 0.1 */
    /* First-launch wizard */
    int  show_wizard;         /* 1 = wizard visible, 0 = normal UI */
    int  wizard_step;         /* 0=welcome, 1=index, 2=formats, 3=ready */
    int  index_found;         /* auto-detected on launch */
    volatile int wizard_building; /* 1 = index build in progress */
} gui_state_t;

/* ================================================================== */
/* Colour helpers                                                      */
/* ================================================================== */
static struct nk_color col_pass       = {34, 139, 34, 255};   /* forest green */
static struct nk_color col_fail       = {220, 20, 60, 255};   /* crimson      */
static struct nk_color col_inconc     = {218, 165, 32, 255};  /* goldenrod    */
static struct nk_color col_halal      = {34, 139, 34, 255};
static struct nk_color col_haram      = {220, 20, 60, 255};
static struct nk_color col_mashbooh   = {218, 165, 32, 255};
static struct nk_color col_unknown    = {160, 160, 160, 255};

static struct nk_color status_color(halal_status_t s) {
    switch (s) {
        case HALAL:            return col_halal;
        case HARAM:            return col_haram;
        case MASHBOOH:         return col_mashbooh;
        default:               return col_unknown;
    }
}

static struct nk_color verdict_color(verdict_t v) {
    switch (v) {
        case PASS:         return col_pass;
        case FAIL:         return col_fail;
        default:           return col_inconc;
    }
}

/* ================================================================== */
/* Friendly species name lookup                                        */
/* ================================================================== */
static const char *friendly_species_name(const char *species_id) {
    /* Map Latin names to common food names */
    struct { const char *latin; const char *common; } names[] = {
        { "Bos_taurus",       "Beef (Cow)"           },
        { "Sus_scrofa",       "Pork (Pig)"           },
        { "Ovis_aries",       "Lamb (Sheep)"         },
        { "Gallus_gallus",    "Chicken"              },
        { "Capra_hircus",     "Goat"                 },
        { "Equus_caballus",   "Horse"                },
        { "Bubalus_bubalis",  "Buffalo"              },
        { "Anas_platyrhynchos","Duck"                },
        { "Cervus_elaphus",   "Deer (Venison)"       },
        { "Meleagris_gallopavo","Turkey"             },
        { "Oryctolagus_cuniculus","Rabbit"            },
        { "Camelus_dromedarius","Camel"              },
        { "Canis_lupus",      "Dog"                  },
        { "Equus_asinus",     "Donkey"               },
        { NULL, NULL }
    };
    for (int i = 0; names[i].latin; i++) {
        if (strcmp(species_id, names[i].latin) == 0)
            return names[i].common;
    }
    return species_id;  /* fallback to Latin name */
}

/* Friendly verdict text */
static const char *friendly_verdict(verdict_t v) {
    switch (v) {
        case PASS: return "HALAL - No haram content detected";
        case FAIL: return "NOT HALAL - Haram content detected";
        default:   return "INCONCLUSIVE - Unable to determine";
    }
}

/* Friendly status text */
static const char *friendly_status(halal_status_t s) {
    switch (s) {
        case HALAL:    return "Halal";
        case HARAM:    return "Haram";
        case MASHBOOH: return "Doubtful";
        default:       return "Unknown";
    }
}

/* Confidence level from cross-marker agreement */
static const char *confidence_label(double agreement) {
    if (agreement >= 0.95) return "Very High";
    if (agreement >= 0.85) return "High";
    if (agreement >= 0.70) return "Moderate";
    return "Low";
}

/* ================================================================== */
/* File dialog — multi-select                                          */
/* ================================================================== */
static void open_file_dialog(gui_state_t *st) {
    const char *filters[] = { "*.fq", "*.fastq", "*.fq.gz", "*.fastq.gz",
                               "*.fa", "*.fasta", "*.fa.gz", "*.fasta.gz" };
    const char *result = tinyfd_openFileDialog(
        "Select DNA sample file(s)",           /* title        */
        "",                                    /* default path */
        8, filters,                            /* filter       */
        "DNA sample files (*.fq *.fa *.gz)",   /* description  */
        1                                      /* multi-select */
    );
    if (!result) return;

    /* Parse pipe-separated paths (tinyfiledialogs multi-select format) */
    analysis_context_t *a = &st->analysis;
    a->n_fastq_files = 0;
    const char *p = result;
    while (*p && a->n_fastq_files < 32) {
        const char *sep = strchr(p, '|');
        int len = sep ? (int)(sep - p) : (int)strlen(p);
        if (len > 0 && len < 1024) {
            memcpy(a->fastq_paths[a->n_fastq_files], p, (size_t)len);
            a->fastq_paths[a->n_fastq_files][len] = '\0';
            a->n_fastq_files++;
        }
        if (!sep) break;
        p = sep + 1;
    }
    detect_samples(a);
}

/* ================================================================== */
/* Locate default index file                                           */
/* ================================================================== */
static void find_default_index(char *path, size_t sz) {
#ifdef __APPLE__
    {
        char buf[1024];
        const char *base = SDL_GetBasePath();
        if (base) {
            snprintf(buf, sizeof(buf), "%s../Resources/default.idx", base);
            FILE *f = fopen(buf, "rb");
            if (f) { fclose(f); snprintf(path, sz, "%s", buf); return; }
        }
    }
#endif
    {
        const char *base = SDL_GetBasePath();
        if (base) {
            char buf[1024];
            snprintf(buf, sizeof(buf), "%sdefault.idx", base);
            FILE *f = fopen(buf, "rb");
            if (f) { fclose(f); snprintf(path, sz, "%s", buf); return; }
        }
    }
    {
        FILE *f = fopen("halal.idx", "rb");
        if (f) { fclose(f); snprintf(path, sz, "halal.idx"); return; }
    }
    path[0] = '\0';
}

/* ================================================================== */
/* Load database info from index for database viewer                   */
/* ================================================================== */
static halal_refdb_t *load_db_info(const char *index_path) {
    if (!index_path[0]) return NULL;
    halal_index_t *idx = index_load(index_path);
    if (!idx) return NULL;
    /* Steal the db pointer, then free index parts manually
       (index_destroy dereferences idx->db, so we can't NULL it) */
    halal_refdb_t *db = idx->db;
    int S = db->n_species;
    int M = db->n_markers;
    for (int s = 0; s < S; s++) fmh_destroy(idx->coarse[s]);
    free(idx->coarse);
    for (int m = 0; m < M; m++) {
        for (int s = 0; s < S; s++) kmer_set_destroy(idx->fine[m][s]);
        free(idx->fine[m]);
    }
    free(idx->fine);
    for (int m = 0; m < M; m++) kmer_set_destroy(idx->primer_index[m]);
    free(idx->primer_index);
    free(idx);  /* don't call refdb_destroy — we're keeping db */
    return db;
}

/* ================================================================== */
/* Format file size                                                    */
/* ================================================================== */
static void format_bytes(long bytes, char *out, size_t sz) {
    if (bytes < 1024)
        snprintf(out, sz, "%ld B", bytes);
    else if (bytes < 1024 * 1024)
        snprintf(out, sz, "%.1f KB", bytes / 1024.0);
    else if (bytes < 1024L * 1024 * 1024)
        snprintf(out, sz, "%.1f MB", bytes / (1024.0 * 1024.0));
    else
        snprintf(out, sz, "%.1f GB", bytes / (1024.0 * 1024.0 * 1024.0));
}

/* ================================================================== */
/* Draw horizontal bar chart with CI whiskers                          */
/* ================================================================== */
static void draw_horizontal_bars(struct nk_context *ctx,
                                  const halal_report_t *r,
                                  float chart_width,
                                  float min_display_pct)
{
    struct nk_command_buffer *canvas = nk_window_get_canvas(ctx);
    float name_col = 140.0f;
    float pct_col = 64.0f;
    float bar_area = chart_width - name_col - pct_col - 20.0f;
    if (bar_area < 40.0f) bar_area = 40.0f;

    /* Count hidden species (below threshold but non-zero) */
    int n_hidden = 0;
    for (int i = 0; i < r->n_species; i++) {
        const species_report_t *sp = &r->species[i];
        if (sp->weight_pct < min_display_pct && sp->read_pct < min_display_pct
            && (sp->weight_pct > 0.0001 || sp->read_pct > 0.0001))
            n_hidden++;
    }

    /* Find max value for scaling (among visible species) */
    double max_val = 0.0;
    for (int i = 0; i < r->n_species; i++) {
        const species_report_t *sp = &r->species[i];
        if (sp->weight_pct < min_display_pct && sp->read_pct < min_display_pct)
            continue;
        double hi = sp->ci_hi;
        if (hi > max_val) max_val = hi;
        if (sp->weight_pct > max_val) max_val = sp->weight_pct;
    }
    if (max_val < 1.0) max_val = 1.0;
    double scale_max = max_val * 1.1;
    if (scale_max > 100.0) scale_max = 100.0;

    int row_idx = 0;
    for (int i = 0; i < r->n_species; i++) {
        const species_report_t *sp = &r->species[i];
        if (sp->weight_pct < min_display_pct && sp->read_pct < min_display_pct)
            continue;

        nk_layout_row_dynamic(ctx, 28, 1);
        struct nk_rect row_bounds;
        nk_widget(&row_bounds, ctx);

        float x0 = row_bounds.x;
        float y_mid = row_bounds.y + row_bounds.h * 0.5f;

        /* Alternating row background */
        if (row_idx % 2 == 1) {
            nk_fill_rect(canvas, row_bounds, 0,
                         nk_rgba(255, 255, 255, 8));
        }
        row_idx++;

        /* Species name */
        struct nk_rect name_rect = nk_rect(x0, row_bounds.y,
                                            name_col, row_bounds.h);
        nk_draw_text(canvas, name_rect,
                     friendly_species_name(sp->species_id),
                     (int)strlen(friendly_species_name(sp->species_id)),
                     ctx->style.font, nk_rgba(0,0,0,0),
                     nk_rgb(210, 210, 210));

        /* Bar track (dim background) */
        float bar_x = x0 + name_col;
        float bar_h = row_bounds.h - 8.0f;
        float bar_y = row_bounds.y + 4.0f;
        nk_fill_rect(canvas, nk_rect(bar_x, bar_y, bar_area, bar_h),
                      4, nk_rgba(255, 255, 255, 15));

        /* Colored bar */
        float bar_w = (float)(sp->weight_pct / scale_max) * bar_area;
        if (bar_w < 3.0f && sp->weight_pct > 0.0001) bar_w = 3.0f;
        struct nk_color bar_color = status_color(sp->halal_status);
        if (bar_w > 0.5f)
            nk_fill_rect(canvas, nk_rect(bar_x, bar_y, bar_w, bar_h),
                          4, bar_color);

        /* CI whiskers (bar color at 40% opacity) */
        if (sp->ci_lo >= 0 && sp->ci_hi > 0) {
            float ci_lo_x = bar_x + (float)(sp->ci_lo / scale_max) * bar_area;
            float ci_hi_x = bar_x + (float)(sp->ci_hi / scale_max) * bar_area;
            float whisker_h = bar_h * 0.6f;
            struct nk_color wc = nk_rgba(bar_color.r, bar_color.g,
                                          bar_color.b, 100);
            nk_stroke_line(canvas, ci_lo_x, y_mid, ci_hi_x, y_mid, 2.0f, wc);
            nk_stroke_line(canvas, ci_lo_x, y_mid - whisker_h/2,
                           ci_lo_x, y_mid + whisker_h/2, 1.5f, wc);
            nk_stroke_line(canvas, ci_hi_x, y_mid - whisker_h/2,
                           ci_hi_x, y_mid + whisker_h/2, 1.5f, wc);
        }

        /* Percentage text (right-aligned) */
        char pct_buf[32];
        if (sp->weight_pct > 0.0001 && sp->weight_pct < 0.1)
            snprintf(pct_buf, sizeof(pct_buf), "< 0.1%%");
        else
            snprintf(pct_buf, sizeof(pct_buf), "%5.1f%%", sp->weight_pct);
        float pct_x = x0 + name_col + bar_area + 4.0f;
        struct nk_rect pct_rect = nk_rect(pct_x, row_bounds.y,
                                           pct_col, row_bounds.h);
        nk_draw_text(canvas, pct_rect, pct_buf, (int)strlen(pct_buf),
                     ctx->style.font, nk_rgba(0,0,0,0),
                     nk_rgb(210, 210, 210));
    }

    /* "N more below threshold" */
    if (n_hidden > 0) {
        nk_layout_row_dynamic(ctx, 22, 1);
        char hidden_buf[64];
        snprintf(hidden_buf, sizeof(hidden_buf),
                 "%d more below %.1f%% threshold", n_hidden, min_display_pct);
        nk_label_colored(ctx, hidden_buf, NK_TEXT_LEFT,
                         nk_rgb(120, 120, 120));
    }
}

/* ================================================================== */
/* Draw stacked species bar (bottom of left panel)                     */
/* ================================================================== */
static void draw_species_bar(struct nk_context *ctx, const halal_report_t *r) {
    struct nk_command_buffer *canvas = nk_window_get_canvas(ctx);
    struct nk_rect bounds;
    nk_layout_row_dynamic(ctx, ROW_BAR, 1);
    nk_widget(&bounds, ctx);

    /* Rounded background track */
    nk_fill_rect(canvas, bounds, 4, nk_rgba(255, 255, 255, 15));

    /* Count visible segments */
    int n_vis = 0;
    for (int i = 0; i < r->n_species; i++) {
        float w = (float)(r->species[i].weight_pct / 100.0) * bounds.w;
        if (w >= 1.0f) n_vis++;
    }

    /* Draw segments with 1px dark gap between them */
    float x = bounds.x;
    int seg = 0;
    for (int i = 0; i < r->n_species; i++) {
        float w = (float)(r->species[i].weight_pct / 100.0) * bounds.w;
        if (w < 1.0f) continue;

        /* 1px gap between segments (not before first or after last) */
        float gap = (seg > 0 && seg < n_vis) ? 1.0f : 0.0f;
        float draw_x = x + gap;
        float draw_w = w - gap;
        if (draw_w < 1.0f) draw_w = 1.0f;

        struct nk_color c = status_color(r->species[i].halal_status);

        /* Round left end of first segment, right end of last */
        float rounding = 0;
        if (seg == 0 && seg == n_vis - 1)
            rounding = 4;  /* only segment: round both ends */
        else if (seg == 0)
            rounding = 4;  /* first segment */
        else if (seg == n_vis - 1)
            rounding = 4;  /* last segment */

        nk_fill_rect(canvas, nk_rect(draw_x, bounds.y, draw_w, bounds.h),
                      rounding, c);
        x += w;
        seg++;
    }
}

/* ================================================================== */
/* Draw database viewer panel                                          */
/* ================================================================== */
static void draw_database_panel(struct nk_context *ctx,
                                 const halal_refdb_t *db)
{
    if (!db) {
        nk_layout_row_dynamic(ctx, 50, 1);
        nk_label_wrap(ctx, "No database loaded. Ensure the index file "
                      "is available.");
        return;
    }

    /* Summary line */
    nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
    {
        char summary[128];
        snprintf(summary, sizeof(summary),
                 "%d species, %d markers, %d references",
                 db->n_species, db->n_markers, db->n_marker_refs);
        nk_label(ctx, summary, NK_TEXT_LEFT);
    }

    nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
    nk_spacing(ctx, 1);

    /* --- Species table --- */
    nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
    nk_label(ctx, "Species", NK_TEXT_LEFT);

    /* Header */
    nk_layout_row_dynamic(ctx, ROW_TABLE, 4);
    nk_label(ctx, "Species ID",   NK_TEXT_LEFT);
    nk_label(ctx, "Common Name",  NK_TEXT_LEFT);
    nk_label(ctx, "Status",       NK_TEXT_CENTERED);
    nk_label(ctx, "Mito CN",      NK_TEXT_RIGHT);

    for (int s = 0; s < db->n_species; s++) {
        nk_layout_row_dynamic(ctx, ROW_TABLE, 4);
        nk_label(ctx, db->species[s].species_id, NK_TEXT_LEFT);
        nk_label(ctx, db->species[s].common_name, NK_TEXT_LEFT);
        nk_label_colored(ctx, friendly_status(db->species[s].status),
                         NK_TEXT_CENTERED,
                         status_color(db->species[s].status));
        char cn_buf[32];
        snprintf(cn_buf, sizeof(cn_buf), "%.0f",
                 db->species[s].mito_copy_number);
        nk_label(ctx, cn_buf, NK_TEXT_RIGHT);
    }

    nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
    nk_spacing(ctx, 1);

    /* --- Markers table --- */
    nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
    nk_label(ctx, "Markers & Primers", NK_TEXT_LEFT);

    nk_layout_row_dynamic(ctx, ROW_TABLE, 3);
    nk_label(ctx, "Marker",       NK_TEXT_LEFT);
    nk_label(ctx, "Forward",      NK_TEXT_LEFT);
    nk_label(ctx, "Reverse",      NK_TEXT_LEFT);

    for (int m = 0; m < db->n_markers; m++) {
        nk_layout_row_dynamic(ctx, ROW_TABLE, 3);
        nk_label(ctx, db->marker_ids[m], NK_TEXT_LEFT);
        nk_label(ctx, db->primer_f[m][0] ? db->primer_f[m] : "-",
                 NK_TEXT_LEFT);
        nk_label(ctx, db->primer_r[m][0] ? db->primer_r[m] : "-",
                 NK_TEXT_LEFT);
    }

    nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
    nk_spacing(ctx, 1);

    /* --- Coverage matrix --- */
    nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
    nk_label(ctx, "Reference Coverage (amplicon bp)", NK_TEXT_LEFT);

    /* Header row: blank + marker names */
    int n_cols = db->n_markers + 1;
    nk_layout_row_dynamic(ctx, ROW_TABLE, n_cols);
    nk_label(ctx, "Species", NK_TEXT_LEFT);
    for (int m = 0; m < db->n_markers; m++)
        nk_label(ctx, db->marker_ids[m], NK_TEXT_CENTERED);

    /* Data rows */
    for (int s = 0; s < db->n_species; s++) {
        nk_layout_row_dynamic(ctx, ROW_TABLE, n_cols);
        nk_label(ctx, friendly_species_name(db->species[s].species_id),
                 NK_TEXT_LEFT);
        for (int m = 0; m < db->n_markers; m++) {
            marker_ref_t *mr = refdb_get_marker_ref(db, s, m);
            char cell[16];
            if (mr && mr->seq_len > 0)
                snprintf(cell, sizeof(cell), "%d", mr->seq_len);
            else
                snprintf(cell, sizeof(cell), "-");
            nk_label(ctx, cell, NK_TEXT_CENTERED);
        }
    }
}

/* ================================================================== */
/* First-launch wizard                                                 */
/* ================================================================== */

/* Check if ~/.halalseq/setup_done exists */
static int wizard_setup_done_exists(void) {
    const char *home = getenv("HOME");
    if (!home) return 0;
    char path[1024];
    snprintf(path, sizeof(path), "%s/.halalseq/setup_done", home);
    FILE *f = fopen(path, "r");
    if (f) { fclose(f); return 1; }
    return 0;
}

/* Create ~/.halalseq/setup_done marker */
static void wizard_write_setup_done(void) {
    const char *home = getenv("HOME");
    if (!home) return;
    char dir[1024];
    snprintf(dir, sizeof(dir), "%s/.halalseq", home);
    mkdir(dir, 0755);
    char path[1024];
    snprintf(path, sizeof(path), "%s/.halalseq/setup_done", home);
    FILE *f = fopen(path, "w");
    if (f) { fprintf(f, "1\n"); fclose(f); }
}

/* Background thread for index build from wizard */
static int wizard_build_worker(void *data) {
    gui_state_t *st = (gui_state_t *)data;

    /* Build refdb then index via system() calls to halalseq CLI */
    const char *base = SDL_GetBasePath();
    char cmd[2048];
    if (base && base[0]) {
        snprintf(cmd, sizeof(cmd),
                 "%shalalseq build-db -o /tmp/_hs_wizard.db && "
                 "%shalalseq index -d /tmp/_hs_wizard.db -o halal.idx && "
                 "rm -f /tmp/_hs_wizard.db",
                 base, base);
    } else {
        snprintf(cmd, sizeof(cmd),
                 "./halalseq build-db -o /tmp/_hs_wizard.db && "
                 "./halalseq index -d /tmp/_hs_wizard.db -o halal.idx && "
                 "rm -f /tmp/_hs_wizard.db");
    }
    int ret = system(cmd);
    if (ret == 0) {
        find_default_index(st->index_path, sizeof(st->index_path));
        if (st->index_path[0]) {
            st->index_found = 1;
            if (st->db_info) refdb_destroy(st->db_info);
            st->db_info = load_db_info(st->index_path);
        }
    }
    st->wizard_building = 0;
    return ret;
}

/* Draw the wizard overlay */
static void draw_wizard(struct nk_context *ctx, gui_state_t *st,
                         int win_w, int win_h)
{
    /* Centered panel ~500x420 */
    float pw = 500, ph = 420;
    if (pw > win_w - 40) pw = (float)win_w - 40;
    if (ph > win_h - 40) ph = (float)win_h - 40;
    float px = ((float)win_w - pw) * 0.5f;
    float py = ((float)win_h - ph) * 0.5f;

    if (!nk_begin(ctx, "Setup Wizard",
                  nk_rect(px, py, pw, ph),
                  NK_WINDOW_BORDER | NK_WINDOW_TITLE |
                  NK_WINDOW_NO_SCROLLBAR)) {
        nk_end(ctx);
        return;
    }

    /* Step indicator — ASCII-safe, e.g. "Step 1 of 4  [*] [ ] [ ] [ ]" */
    nk_layout_row_dynamic(ctx, 20, 1);
    {
        char dot_str[64] = "";
        char tmp[16];
        snprintf(tmp, sizeof(tmp), "Step %d/4  ", st->wizard_step + 1);
        strcat(dot_str, tmp);
        for (int i = 0; i < 4; i++) {
            strcat(dot_str, i == st->wizard_step ? "[*] " : "[ ] ");
        }
        nk_label(ctx, dot_str, NK_TEXT_CENTERED);
    }

    nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
    nk_spacing(ctx, 1);

    /* Content area per step */
    switch (st->wizard_step) {
    case 0: /* Welcome */
        nk_layout_row_dynamic(ctx, 36, 1);
        nk_label(ctx, "Welcome to HalalSeq", NK_TEXT_CENTERED);

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        nk_layout_row_dynamic(ctx, 60, 1);
        nk_label_wrap(ctx,
            "Halal food authentication via DNA metabarcoding. "
            "This wizard will check that everything is set up "
            "for accurate analysis.");

        nk_layout_row_dynamic(ctx, 44, 1);
        nk_label_wrap(ctx,
            "HalalSeq identifies animal species in food samples "
            "using mitochondrial DNA markers.");
        break;

    case 1: /* Reference Index */
        nk_layout_row_dynamic(ctx, 36, 1);
        nk_label(ctx, "Reference Index", NK_TEXT_CENTERED);

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        if (st->wizard_building) {
            nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
            nk_label_colored(ctx, "Building reference index...",
                             NK_TEXT_LEFT, nk_rgb(218, 165, 32));
            nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
            nk_size bpv = 50;
            nk_progress(ctx, &bpv, 100, NK_FIXED);
        } else if (st->index_found) {
            nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
            nk_label_colored(ctx, "Index found",
                             NK_TEXT_LEFT, nk_rgb(34, 139, 34));

            nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
            nk_label(ctx, st->index_path, NK_TEXT_LEFT);

            if (st->db_info) {
                nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                char info[128];
                snprintf(info, sizeof(info),
                         "%d species, %d markers",
                         st->db_info->n_species,
                         st->db_info->n_markers);
                nk_label(ctx, info, NK_TEXT_LEFT);
            }
        } else {
            nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
            nk_label_colored(ctx,
                "No reference index found!",
                NK_TEXT_LEFT, nk_rgb(220, 20, 60));

            nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
            nk_label_wrap(ctx,
                "An index is required for species identification. "
                "Click below to build one from built-in references.");

            nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
            nk_spacing(ctx, 1);

            nk_layout_row_dynamic(ctx, ROW_BTN, 1);
            if (nk_button_label(ctx, "Build Index")) {
                st->wizard_building = 1;
                SDL_CreateThread(wizard_build_worker, "wizard_build", st);
            }
        }
        break;

    case 2: /* Supported Formats */
        nk_layout_row_dynamic(ctx, 36, 1);
        nk_label(ctx, "Supported Formats", NK_TEXT_CENTERED);

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
        nk_label(ctx, "Supported file types:", NK_TEXT_LEFT);

        nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
        nk_label(ctx, "  .fq  .fastq  .fq.gz  .fastq.gz", NK_TEXT_LEFT);
        nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
        nk_label(ctx, "  .fa  .fasta  .fa.gz  .fasta.gz", NK_TEXT_LEFT);

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
        nk_label_wrap(ctx,
            "HalalSeq works with raw sequencing data - "
            "no preprocessing required.");

        nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
        nk_label_wrap(ctx,
            "For best results, use amplicon-targeted sequencing "
            "(PCR + Illumina).");

        nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
        nk_label_wrap(ctx,
            "R1/R2 paired-end files are automatically detected "
            "and merged.");
        break;

    case 3: /* Ready */
        nk_layout_row_dynamic(ctx, 36, 1);
        nk_label(ctx, "Setup Complete", NK_TEXT_CENTERED);

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        nk_layout_row_dynamic(ctx, 60, 1);
        nk_label_wrap(ctx,
            "Everything is ready! Click Start to begin using "
            "HalalSeq for halal food authentication.");

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
        nk_label_colored(ctx, "You can re-run this wizard by deleting:",
                         NK_TEXT_LEFT, nk_rgb(140, 140, 140));
        nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
        nk_label_colored(ctx, "  ~/.halalseq/setup_done",
                         NK_TEXT_LEFT, nk_rgb(140, 140, 140));
        break;
    }

    /* Spacer to push nav buttons to bottom */
    nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
    nk_spacing(ctx, 1);
    nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
    nk_spacing(ctx, 1);

    /* Navigation buttons */
    nk_layout_row_dynamic(ctx, ROW_BTN, 2);

    /* Back button */
    if (st->wizard_step > 0) {
        if (nk_button_label(ctx, "Back"))
            st->wizard_step--;
    } else {
        /* Invisible placeholder */
        struct nk_style_button invis = ctx->style.button;
        invis.normal = nk_style_item_color(nk_rgba(0,0,0,0));
        invis.hover  = nk_style_item_color(nk_rgba(0,0,0,0));
        invis.active = nk_style_item_color(nk_rgba(0,0,0,0));
        invis.border = 0;
        invis.text_normal = nk_rgba(0,0,0,0);
        invis.text_hover  = nk_rgba(0,0,0,0);
        invis.text_active = nk_rgba(0,0,0,0);
        nk_button_label_styled(ctx, &invis, "");
    }

    /* Next / Start button */
    if (st->wizard_step < 3) {
        /* Disable Next on step 1 if building or no index */
        int can_next = 1;
        if (st->wizard_step == 1 && !st->index_found)
            can_next = 0;
        if (st->wizard_building)
            can_next = 0;

        if (can_next) {
            struct nk_style_button green = ctx->style.button;
            green.normal = nk_style_item_color(nk_rgb(34, 139, 34));
            green.hover  = nk_style_item_color(nk_rgb(0, 180, 0));
            green.text_normal = nk_rgb(255, 255, 255);
            green.text_hover  = nk_rgb(255, 255, 255);
            if (nk_button_label_styled(ctx, &green, "Next"))
                st->wizard_step++;
        } else {
            struct nk_style_button grey = ctx->style.button;
            grey.normal = nk_style_item_color(nk_rgb(80, 80, 80));
            grey.hover  = nk_style_item_color(nk_rgb(80, 80, 80));
            grey.text_normal = nk_rgb(140, 140, 140);
            grey.text_hover  = nk_rgb(140, 140, 140);
            nk_button_label_styled(ctx, &grey, "Next");
        }
    } else {
        /* Start button on final step */
        struct nk_style_button green = ctx->style.button;
        green.normal = nk_style_item_color(nk_rgb(34, 139, 34));
        green.hover  = nk_style_item_color(nk_rgb(0, 180, 0));
        green.text_normal = nk_rgb(255, 255, 255);
        green.text_hover  = nk_rgb(255, 255, 255);
        if (nk_button_label_styled(ctx, &green, "Start")) {
            st->show_wizard = 0;
            wizard_write_setup_done();
        }
    }

    nk_end(ctx);
}

/* ================================================================== */
/* Draw the full GUI layout                                            */
/* ================================================================== */
static void draw_gui(struct nk_context *ctx, gui_state_t *st,
                     int win_w, int win_h)
{
    /* Show wizard instead of normal UI if active */
    if (st->show_wizard) {
        draw_wizard(ctx, st, win_w, win_h);
        return;
    }

    analysis_context_t *analysis = &st->analysis;

    if (!nk_begin(ctx, "HalalSeq",
                  nk_rect(0, 0, (float)win_w, (float)win_h),
                  NK_WINDOW_NO_SCROLLBAR | NK_WINDOW_BACKGROUND)) {
        nk_end(ctx);
        return;
    }

    /* Title bar */
    nk_layout_row_dynamic(ctx, ROW_TITLE, 1);
    nk_label(ctx, "HalalSeq - Halal Food DNA Authentication", NK_TEXT_CENTERED);

    nk_layout_row_dynamic(ctx, 2, 1);
    nk_spacing(ctx, 1);

    /* Two-column layout */
    float col_widths[] = { (float)LEFT_W, (float)(win_w - LEFT_W - 20) };
    nk_layout_row(ctx, NK_STATIC, (float)(win_h - ROW_TITLE - 24), 2, col_widths);

    /* ============================================================== */
    /* LEFT PANEL                                                      */
    /* ============================================================== */
    if (nk_group_begin(ctx, "input_panel", NK_WINDOW_BORDER)) {

        nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
        nk_label(ctx, "DNA Sample Files:", NK_TEXT_LEFT);

        /* Sample list display */
        if (analysis->n_samples > 0) {
            int show = analysis->n_samples > 6 ? 6 : analysis->n_samples;
            for (int i = 0; i < show; i++) {
                nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                char label[300];
                if (analysis->samples[i].n_files == 2)
                    snprintf(label, sizeof(label), "%s (R1+R2)",
                             analysis->samples[i].sample_name);
                else
                    snprintf(label, sizeof(label), "%s",
                             analysis->samples[i].sample_name);
                nk_label(ctx, label, NK_TEXT_LEFT);
            }
            if (analysis->n_samples > 6) {
                nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                char more[32];
                snprintf(more, sizeof(more), "... +%d more",
                         analysis->n_samples - 6);
                nk_label(ctx, more, NK_TEXT_LEFT);
            }

            /* Total file size & sample count */
            {
                mem_estimate_t est = analysis_estimate_memory(analysis);
                char sz_buf[64];
                format_bytes(est.total_file_bytes, sz_buf, sizeof(sz_buf));
                nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                char info[128];
                snprintf(info, sizeof(info),
                         "%d sample%s (%d file%s, %s)",
                         analysis->n_samples,
                         analysis->n_samples > 1 ? "s" : "",
                         analysis->n_fastq_files,
                         analysis->n_fastq_files > 1 ? "s" : "",
                         sz_buf);
                nk_label(ctx, info, NK_TEXT_LEFT);

                /* RAM estimate */
                nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                char ram_info[128];
                snprintf(ram_info, sizeof(ram_info),
                         "Est. reads: %dK, RAM: %d MB",
                         est.estimated_reads / 1000,
                         est.estimated_ram_mb);
                nk_label(ctx, ram_info, NK_TEXT_LEFT);

                /* RAM warning */
                if (est.estimated_ram_mb > 1024) {
                    nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                    nk_label_colored(ctx,
                        "Warning: >1 GB RAM estimated",
                        NK_TEXT_LEFT, nk_rgb(255, 200, 50));

                    nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                    nk_checkbox_label(ctx, "Subsample to 500K reads/sample",
                                     &analysis->subsample_enabled);
                }
            }
        } else {
            nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
            nk_label(ctx, "(drop files or click Choose Files)",
                     NK_TEXT_LEFT);
        }

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        /* Browse + Clear buttons */
        nk_layout_row_dynamic(ctx, ROW_BTN, 2);
        if (nk_button_label(ctx, "Choose Files...")) {
            open_file_dialog(st);
        }
        {
            int has_files = analysis->n_fastq_files > 0;
            if (has_files && nk_button_label(ctx, "Clear Files")) {
                analysis->n_fastq_files = 0;
                analysis->n_samples = 0;
            } else if (!has_files) {
                /* Draw disabled clear button */
                struct nk_style_button grey = ctx->style.button;
                grey.normal = nk_style_item_color(nk_rgb(50, 50, 50));
                grey.hover  = nk_style_item_color(nk_rgb(50, 50, 50));
                grey.text_normal = nk_rgb(100, 100, 100);
                grey.text_hover  = nk_rgb(100, 100, 100);
                nk_button_label_styled(ctx, &grey, "Clear Files");
            }
        }

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        /* Analyze button */
        nk_layout_row_dynamic(ctx, ROW_BTN, 1);
        {
            int can_run = (analysis->n_fastq_files > 0 &&
                          (analysis->state == ANALYSIS_IDLE ||
                           analysis->state == ANALYSIS_DONE ||
                           analysis->state == ANALYSIS_ERROR) &&
                          st->index_path[0]);
            if (can_run) {
                struct nk_style_button green = ctx->style.button;
                green.normal  = nk_style_item_color(nk_rgb(34, 139, 34));
                green.hover   = nk_style_item_color(nk_rgb(0, 180, 0));
                green.active  = nk_style_item_color(nk_rgb(0, 140, 0));
                green.text_normal  = nk_rgb(255, 255, 255);
                green.text_hover   = nk_rgb(255, 255, 255);
                green.text_active  = nk_rgb(255, 255, 255);
                if (nk_button_label_styled(ctx, &green, "Run Analysis")) {
                    snprintf(analysis->index_path,
                             sizeof(analysis->index_path), "%s",
                             st->index_path);
                    analysis_start(analysis);
                }
            } else {
                struct nk_style_button grey = ctx->style.button;
                grey.normal = nk_style_item_color(nk_rgb(80, 80, 80));
                grey.hover  = nk_style_item_color(nk_rgb(80, 80, 80));
                grey.active = nk_style_item_color(nk_rgb(80, 80, 80));
                grey.text_normal = nk_rgb(140, 140, 140);
                grey.text_hover  = nk_rgb(140, 140, 140);
                grey.text_active = nk_rgb(140, 140, 140);
                nk_button_label_styled(ctx, &grey, "Run Analysis");
            }
        }

        /* Status line — friendly labels with per-sample progress */
        nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
        {
            static char status_buf[128];
            const char *status_text;
            int si = analysis->progress_sample_idx;
            int ns = analysis->n_samples;
            const char *sname = (ns > 0 && si < ns)
                ? analysis->samples[si].sample_name : "";

            switch (analysis->state) {
                case ANALYSIS_IDLE:              status_text = "Ready"; break;
                case ANALYSIS_LOADING_INDEX:     status_text = "Preparing..."; break;
                case ANALYSIS_READING_FASTQ:
                    if (ns > 1)
                        snprintf(status_buf, sizeof(status_buf),
                                 "Sample %d/%d: %s (%d reads)...",
                                 si + 1, ns, sname,
                                 analysis->progress_reads);
                    else
                        snprintf(status_buf, sizeof(status_buf),
                                 "Reading sample (%d reads)...",
                                 analysis->progress_reads);
                    status_text = status_buf;
                    break;
                case ANALYSIS_CLASSIFYING:
                    if (ns > 1)
                        snprintf(status_buf, sizeof(status_buf),
                                 "Sample %d/%d: Identifying species...",
                                 si + 1, ns);
                    else
                        snprintf(status_buf, sizeof(status_buf),
                                 "Identifying species...");
                    status_text = status_buf;
                    break;
                case ANALYSIS_RUNNING_EM:
                    if (ns > 1)
                        snprintf(status_buf, sizeof(status_buf),
                                 "Sample %d/%d: Calculating amounts...",
                                 si + 1, ns);
                    else
                        snprintf(status_buf, sizeof(status_buf),
                                 "Calculating amounts...");
                    status_text = status_buf;
                    break;
                case ANALYSIS_GENERATING_REPORT:
                    if (ns > 1)
                        snprintf(status_buf, sizeof(status_buf),
                                 "Sample %d/%d: Generating report...",
                                 si + 1, ns);
                    else
                        snprintf(status_buf, sizeof(status_buf),
                                 "Generating report...");
                    status_text = status_buf;
                    break;
                case ANALYSIS_DONE:              status_text = "Analysis complete"; break;
                case ANALYSIS_ERROR:             status_text = "Error occurred"; break;
                default:                         status_text = ""; break;
            }
            nk_label(ctx, status_text, NK_TEXT_LEFT);
        }

        /* Progress bar */
        if (analysis->state > ANALYSIS_IDLE &&
            analysis->state < ANALYSIS_DONE) {
            nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
            nk_size pv = 0;
            switch (analysis->state) {
                case ANALYSIS_LOADING_INDEX:     pv = 10; break;
                case ANALYSIS_READING_FASTQ:     pv = 25; break;
                case ANALYSIS_CLASSIFYING:       pv = 50; break;
                case ANALYSIS_RUNNING_EM:        pv = 75; break;
                case ANALYSIS_GENERATING_REPORT: pv = 90; break;
                default: break;
            }
            nk_progress(ctx, &pv, 100, NK_FIXED);
        }

        /* Stacked species bar + legend (after results) */
        {
            halal_report_t *sel_rpt = NULL;
            if (analysis->state == ANALYSIS_DONE &&
                analysis->reports &&
                analysis->selected_sample < analysis->n_reports)
                sel_rpt = analysis->reports[analysis->selected_sample];
            if (sel_rpt) {
                nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
                nk_spacing(ctx, 1);
                nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
                nk_label(ctx, "Sample Composition:", NK_TEXT_LEFT);
                draw_species_bar(ctx, sel_rpt);

                /* Legend: show species names next to color blocks */
                nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
                nk_spacing(ctx, 1);
                for (int i = 0; i < sel_rpt->n_species; i++) {
                    species_report_t *sp = &sel_rpt->species[i];
                    if (sp->weight_pct < 0.5) continue;
                    nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                    char legend[128];
                    snprintf(legend, sizeof(legend), "  %s: %.1f%%",
                             friendly_species_name(sp->species_id),
                             sp->weight_pct);
                    nk_label_colored(ctx, legend, NK_TEXT_LEFT,
                                     status_color(sp->halal_status));
                }
            }
        }

        nk_group_end(ctx);
    }

    /* ============================================================== */
    /* RIGHT PANEL                                                     */
    /* ============================================================== */
    if (nk_group_begin(ctx, "results_panel", NK_WINDOW_BORDER)) {

        /* Tab toggle: Results / Database */
        nk_layout_row_dynamic(ctx, ROW_LABEL, 2);
        if (nk_option_label(ctx, "Results", st->right_panel_mode == 0))
            st->right_panel_mode = 0;
        if (nk_option_label(ctx, "Database", st->right_panel_mode == 1))
            st->right_panel_mode = 1;

        nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
        nk_spacing(ctx, 1);

        if (st->right_panel_mode == 1) {
            /* --- Database viewer --- */
            draw_database_panel(ctx, st->db_info);

        } else if (analysis->reports && analysis->state == ANALYSIS_DONE) {
            /* --- Sample selector tabs (if multiple samples) --- */
            if (analysis->n_reports > 1) {
                /* Show up to 8 sample tabs per row */
                int tabs_per_row = analysis->n_reports < 8
                    ? analysis->n_reports : 8;
                nk_layout_row_dynamic(ctx, ROW_BTN, tabs_per_row);
                for (int i = 0; i < analysis->n_reports; i++) {
                    halal_report_t *sr = analysis->reports[i];
                    const char *tab_label = sr ? sr->sample_id : "?";
                    if (analysis->selected_sample == i) {
                        /* Active tab — highlight */
                        struct nk_style_button active = ctx->style.button;
                        active.normal = nk_style_item_color(nk_rgb(60, 120, 60));
                        active.hover  = nk_style_item_color(nk_rgb(70, 140, 70));
                        active.text_normal = nk_rgb(255, 255, 255);
                        active.text_hover  = nk_rgb(255, 255, 255);
                        if (nk_button_label_styled(ctx, &active, tab_label))
                            analysis->selected_sample = i;
                    } else {
                        if (nk_button_label(ctx, tab_label))
                            analysis->selected_sample = i;
                    }
                }
                nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
                nk_spacing(ctx, 1);
            }

            /* --- Results for selected sample --- */
            halal_report_t *rpt = NULL;
            if (analysis->selected_sample < analysis->n_reports)
                rpt = analysis->reports[analysis->selected_sample];

            if (rpt) {
                /* Verdict — large, descriptive */
                nk_layout_row_dynamic(ctx, 44, 1);
                {
                    struct nk_color vc = verdict_color(rpt->verdict);
                    nk_label_colored(ctx, friendly_verdict(rpt->verdict),
                                     NK_TEXT_CENTERED, vc);
                }

                nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
                nk_spacing(ctx, 1);

                /* Threshold slider + horizontal bar chart */
                nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
                {
                    char slider_label[64];
                    snprintf(slider_label, sizeof(slider_label),
                             "Show species above: %.1f%%",
                             st->min_display_pct);
                    nk_label(ctx, slider_label, NK_TEXT_LEFT);
                }
                nk_layout_row_dynamic(ctx, ROW_SMALL, 1);
                nk_slider_float(ctx, 0.0f, &st->min_display_pct, 5.0f, 0.1f);

                nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
                nk_spacing(ctx, 1);

                nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
                nk_label(ctx, "Species Proportions", NK_TEXT_LEFT);

                float right_w = col_widths[1] - 16.0f;
                draw_horizontal_bars(ctx, rpt, right_w,
                                      st->min_display_pct);

                nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
                nk_spacing(ctx, 1);

                /* Column headers — friendly names, 4 columns */
                nk_layout_row_dynamic(ctx, ROW_TABLE, 4);
                nk_label(ctx, "Animal",  NK_TEXT_LEFT);
                nk_label(ctx, "Status",  NK_TEXT_CENTERED);
                nk_label(ctx, "Amount",  NK_TEXT_RIGHT);
                nk_label(ctx, "Range",   NK_TEXT_RIGHT);

                /* Species rows — show all detected */
                for (int s = 0; s < rpt->n_species; s++) {
                    species_report_t *sp = &rpt->species[s];
                    if (sp->weight_pct < 0.001 && sp->read_pct < 0.001)
                        continue;

                    nk_layout_row_dynamic(ctx, ROW_TABLE, 4);

                    /* Common name */
                    nk_label(ctx, friendly_species_name(sp->species_id),
                             NK_TEXT_LEFT);

                    /* Status (colour-coded) */
                    nk_label_colored(ctx, friendly_status(sp->halal_status),
                                     NK_TEXT_CENTERED,
                                     status_color(sp->halal_status));

                    /* Amount */
                    char buf[64];
                    snprintf(buf, sizeof(buf), "%.1f%%", sp->weight_pct);
                    nk_label(ctx, buf, NK_TEXT_RIGHT);

                    /* Range */
                    snprintf(buf, sizeof(buf), "%.1f-%.1f%%",
                             sp->ci_lo, sp->ci_hi);
                    nk_label(ctx, buf, NK_TEXT_RIGHT);
                }

                /* Summary — plain English */
                nk_layout_row_dynamic(ctx, ROW_SPACE, 1);
                nk_spacing(ctx, 1);

                nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
                {
                    char meta[256];
                    snprintf(meta, sizeof(meta),
                             "%d DNA fragments analyzed",
                             rpt->total_reads);
                    nk_label(ctx, meta, NK_TEXT_LEFT);
                }

                if (rpt->cross_marker_agreement > 0) {
                    nk_layout_row_dynamic(ctx, ROW_LABEL, 1);
                    char conf[128];
                    snprintf(conf, sizeof(conf), "Confidence: %s",
                             confidence_label(rpt->cross_marker_agreement));
                    nk_label(ctx, conf, NK_TEXT_LEFT);
                }
            }

        } else if (analysis->state == ANALYSIS_IDLE) {
            nk_layout_row_dynamic(ctx, 50, 1);
            nk_label_wrap(ctx,
                "Choose DNA sample file(s) and click "
                "Run Analysis to check halal status.");
        } else if (analysis->state == ANALYSIS_ERROR) {
            nk_layout_row_dynamic(ctx, 50, 1);
            nk_label_colored_wrap(ctx, analysis->error_msg,
                                  nk_rgb(220, 60, 60));
        } else {
            nk_layout_row_dynamic(ctx, 50, 1);
            nk_label_wrap(ctx, "Analyzing your sample, please wait...");
        }

        nk_group_end(ctx);
    }

    nk_end(ctx);
}

/* ================================================================== */
/* Main                                                                */
/* ================================================================== */
int main(int argc, char *argv[]) {
    (void)argc; (void)argv;

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        return 1;
    }

    SDL_SetHint(SDL_HINT_VIDEO_HIGHDPI_DISABLED, "0");

    SDL_Window *window = SDL_CreateWindow(
        "HalalSeq - Halal Food DNA Authentication",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WINDOW_W, WINDOW_H,
        SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE |
        SDL_WINDOW_ALLOW_HIGHDPI);
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer(
        window, -1,
        SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    /* --- Compute DPI scale and set renderer scale -------------------- */
    float dpi_scale = 1.0f;
    {
        int render_w, window_w;
        SDL_GetRendererOutputSize(renderer, &render_w, NULL);
        SDL_GetWindowSize(window, &window_w, NULL);
        if (window_w > 0)
            dpi_scale = (float)render_w / (float)window_w;
        if (dpi_scale < 1.0f) dpi_scale = 1.0f;
    }
    SDL_RenderSetScale(renderer, dpi_scale, dpi_scale);

    /* --- Nuklear init ---------------------------------------------- */
    struct nk_context *ctx = nk_sdl_init(window, renderer);
    {
        struct nk_font_atlas *atlas;
        struct nk_font_config cfg = nk_font_config(0);
        struct nk_font *font = NULL;

        cfg.oversample_h = 3;
        cfg.oversample_v = 2;

        float font_size = 18.0f * dpi_scale;

        nk_sdl_font_stash_begin(&atlas);

        static const char *font_paths[] = {
            "/System/Library/Fonts/Supplemental/Arial.ttf",
            "/System/Library/Fonts/Helvetica.ttc",
            "C:\\Windows\\Fonts\\arial.ttf",
            "C:\\Windows\\Fonts\\segoeui.ttf",
            NULL
        };
        for (const char **p = font_paths; *p; p++) {
            font = nk_font_atlas_add_from_file(atlas, *p, font_size, &cfg);
            if (font) break;
        }
        if (!font)
            font = nk_font_atlas_add_default(atlas, font_size, &cfg);

        nk_sdl_font_stash_end();

        font->handle.height /= dpi_scale;
        nk_style_set_font(ctx, &font->handle);
    }

    /* Dark theme */
    {
        struct nk_color table[NK_COLOR_COUNT];
        table[NK_COLOR_TEXT]                    = nk_rgba(210, 210, 210, 255);
        table[NK_COLOR_WINDOW]                  = nk_rgba(35, 35, 38, 255);
        table[NK_COLOR_HEADER]                  = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_BORDER]                  = nk_rgba(65, 65, 70, 255);
        table[NK_COLOR_BUTTON]                  = nk_rgba(60, 60, 65, 255);
        table[NK_COLOR_BUTTON_HOVER]            = nk_rgba(75, 75, 80, 255);
        table[NK_COLOR_BUTTON_ACTIVE]           = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_TOGGLE]                  = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_TOGGLE_HOVER]            = nk_rgba(55, 55, 60, 255);
        table[NK_COLOR_TOGGLE_CURSOR]           = nk_rgba(44, 160, 44, 255);
        table[NK_COLOR_SELECT]                  = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_SELECT_ACTIVE]           = nk_rgba(44, 160, 44, 255);
        table[NK_COLOR_SLIDER]                  = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_SLIDER_CURSOR]           = nk_rgba(44, 160, 44, 255);
        table[NK_COLOR_SLIDER_CURSOR_HOVER]     = nk_rgba(60, 180, 60, 255);
        table[NK_COLOR_SLIDER_CURSOR_ACTIVE]    = nk_rgba(34, 140, 34, 255);
        table[NK_COLOR_PROPERTY]                = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_EDIT]                    = nk_rgba(45, 45, 50, 255);
        table[NK_COLOR_EDIT_CURSOR]             = nk_rgba(210, 210, 210, 255);
        table[NK_COLOR_COMBO]                   = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_CHART]                   = nk_rgba(50, 50, 55, 255);
        table[NK_COLOR_CHART_COLOR]             = nk_rgba(44, 160, 44, 255);
        table[NK_COLOR_CHART_COLOR_HIGHLIGHT]   = nk_rgba(255, 0, 0, 255);
        table[NK_COLOR_SCROLLBAR]               = nk_rgba(40, 40, 45, 255);
        table[NK_COLOR_SCROLLBAR_CURSOR]        = nk_rgba(60, 60, 65, 255);
        table[NK_COLOR_SCROLLBAR_CURSOR_HOVER]  = nk_rgba(75, 75, 80, 255);
        table[NK_COLOR_SCROLLBAR_CURSOR_ACTIVE] = nk_rgba(55, 55, 60, 255);
        table[NK_COLOR_TAB_HEADER]              = nk_rgba(50, 50, 55, 255);
        nk_style_from_table(ctx, table);
    }

    /* --- App state ------------------------------------------------- */
    gui_state_t state;
    memset(&state, 0, sizeof(state));
    analysis_init(&state.analysis);
    state.min_display_pct = 0.1f;

    find_default_index(state.index_path, sizeof(state.index_path));
    state.index_found = (state.index_path[0] != '\0') ? 1 : 0;

    /* Load database info for database viewer */
    state.db_info = load_db_info(state.index_path);

    /* First-launch wizard: show if no setup_done marker or no index */
    if (!wizard_setup_done_exists() || !state.index_found) {
        state.show_wizard = 1;
        state.wizard_step = 0;
    }

    /* --- Main loop ------------------------------------------------- */
    int running = 1;
    while (running) {
        SDL_Event evt;
        nk_input_begin(ctx);
        while (SDL_PollEvent(&evt)) {
            if (evt.type == SDL_QUIT) {
                running = 0;
            }
            if (evt.type == SDL_DROPFILE) {
                /* Accumulate dropped files */
                if (state.analysis.n_fastq_files < 32) {
                    snprintf(state.analysis.fastq_paths[state.analysis.n_fastq_files],
                             1024, "%s", evt.drop.file);
                    state.analysis.n_fastq_files++;
                    detect_samples(&state.analysis);
                }
                SDL_free(evt.drop.file);
            }
            nk_sdl_handle_event(&evt);
        }
        nk_input_end(ctx);

        int win_w, win_h;
        SDL_GetWindowSize(window, &win_w, &win_h);

        draw_gui(ctx, &state, win_w, win_h);

        SDL_SetRenderDrawColor(renderer, 30, 30, 30, 255);
        SDL_RenderClear(renderer);
        nk_sdl_render(NK_ANTI_ALIASING_ON);
        SDL_RenderPresent(renderer);
    }

    analysis_cleanup(&state.analysis);
    if (state.db_info) refdb_destroy(state.db_info);
    nk_sdl_shutdown();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
