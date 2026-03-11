# ============================================================
# Fenêtre Analyse Double Brin — Schéma + Séquence colorée
# ============================================================

import tkinter as tk
from tkinter import ttk

try:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                                    NavigationToolbar2Tk)
    MPL = True
except ImportError:
    MPL = False

# ── Palette identique à visualize_frame ──────────────────────
BG_MAIN    = '#F0F9F4'
BG_SURFACE = '#FFFFFF'
BG_PANEL   = '#E8F5EE'
BORDER     = '#B8DFC9'
FG_MAIN    = '#1A3A2A'
FG_DIM     = '#6B9E80'
FG_LABEL   = '#2D6A4F'
ACCENT     = '#1B7A4A'
ACCENT2    = '#2196A6'

MPL_BG   = '#FAFFFE'
MPL_AXES = '#F0F9F4'

C_GREEN  = '#52B788'    # ORF brin +
C_TEAL   = '#48CAE4'    # ORF brin -
C_ORANGE = '#F4A261'    # terminateur / Box-10
C_PURPLE = '#B19CD8'    # Shine-Dalgarno
C_BLUE   = '#2196A6'    # Box-35
C_GOLD   = '#F39C12'    # meilleur ORF / bannière

# Fond séquence colorée
SEQ_ORF_BEST = '#A9DFBF'
SEQ_ORF      = '#D5F5E3'
SEQ_BOX35    = '#AED6F1'
SEQ_BOX10    = '#F5B7B1'
SEQ_TERM     = '#FDEBD0'
SEQ_SD       = '#E8DAEF'


# ─────────────────────────────────────────────────────────────
class DualStrandWindow(tk.Toplevel):

    def __init__(self, parent, seq, results):
        super().__init__(parent)
        self.title('🧬  Analyse Double Brin — Résultats')
        self.geometry('1200x860')
        self.minsize(900, 650)
        self.configure(bg=BG_MAIN)
        self.resizable(True, True)
        self.seq = seq.upper()
        self._prepare()
        self._build_ui()

    # ══════════════════════════════════════════════════════════
    # PRÉPARATION DES DONNÉES
    # ══════════════════════════════════════════════════════════
    def _prepare(self):
        from analysis.orf_finder   import (find_reading_frames,
                                            find_best_coding_orf,
                                            reverse_complement)
        from analysis.motif_finder import (find_promoters,
                                            find_shine_dalgarno,
                                            find_terminators)

        seq = self.seq
        rc  = reverse_complement(seq)
        n   = len(seq)

        def gc_stem_ratio(t):
            s = t['arm1'] + t['arm2']
            return (s.count('G') + s.count('C')) / max(1, len(s))

        def analyse(s, sign):
            frames = [f for f in find_reading_frames(s, 20) if f['strand'] == '+']
            best, scored = find_best_coding_orf(frames, n) if frames else (None, [])
            promo  = find_promoters(s)
            sd     = find_shine_dalgarno(s)
            all_t  = find_terminators(s)
            # Terminateurs riches en GC (tige ≥ 50 %) — filtre biologique
            gc_t   = [t for t in all_t if gc_stem_ratio(t) >= 0.50]
            terms  = gc_t if gc_t else all_t[:4]
            return {
                'sign':        sign,
                'label':       "Brin sens 5'→3' (+)" if sign == '+' else "Brin antisens 3'→5' (-)",
                'seq':         s,
                'frames':      frames,
                'best':        best,
                'scored':      scored,
                'score':       scored[0][0] if scored else 0,
                'promoters':   promo['promoters'],
                'sd':          sd,
                'terminators': terms,
                'all_terms':   all_t,
                'orf_color':   C_GREEN if sign == '+' else C_TEAL,
            }

        self.sp  = analyse(seq, '+')
        self.sm  = analyse(rc,  '-')
        self.n   = n

        if self.sm['score'] > self.sp['score']:
            self.sp['is_best'] = False
            self.sm['is_best'] = True
            self.best = self.sm
        else:
            self.sp['is_best'] = True
            self.sm['is_best'] = False
            self.best = self.sp

    # ══════════════════════════════════════════════════════════
    # INTERFACE
    # ══════════════════════════════════════════════════════════
    def _build_ui(self):
        # Barre titre
        hbar = tk.Frame(self, bg=ACCENT, height=50)
        hbar.pack(fill='x')
        hbar.pack_propagate(False)
        tk.Label(hbar, text='🧬  Analyse Double Brin — Résultats',
                 bg=ACCENT, fg='#FFFFFF',
                 font=('Georgia', 13, 'bold')).pack(side='left', padx=14, pady=8)
        tk.Label(hbar,
                 text=(f"★  Meilleur candidat : {self.best['label']}"
                       f"  •  Score {self.best['score']}/100"),
                 bg=ACCENT, fg='#A8D5B8',
                 font=('Georgia', 9, 'italic')).pack(side='right', padx=14)

        # Notebook onglets
        nb = ttk.Notebook(self)
        nb.pack(fill='both', expand=True, padx=6, pady=6)
        style = ttk.Style()
        style.configure('TNotebook.Tab',
                         background=BG_PANEL, foreground=FG_DIM,
                         font=('Georgia', 9), padding=[10, 5])
        style.map('TNotebook.Tab',
                  background=[('selected', BG_SURFACE)],
                  foreground=[('selected', ACCENT)])

        # Onglet 1 — Carte double brin
        t1 = tk.Frame(nb, bg=BG_MAIN)
        nb.add(t1, text='  📊  Carte Double Brin  ')
        self._build_map_tab(t1)

        # Onglet 2 — Séquence colorée brin +
        t2 = tk.Frame(nb, bg=BG_MAIN)
        nb.add(t2, text="  🔤  Séquence Brin +  ")
        self._build_seq_tab(t2, self.sp)

        # Onglet 3 — Séquence colorée brin -
        t3 = tk.Frame(nb, bg=BG_MAIN)
        nb.add(t3, text="  🔤  Séquence Brin -  ")
        self._build_seq_tab(t3, self.sm)

        # Onglet 4 — Résultat optimal
        t4 = tk.Frame(nb, bg=BG_MAIN)
        nb.add(t4, text='  ★  Résultat Optimal  ')
        self._build_optimal_tab(t4)

    # ══════════════════════════════════════════════════════════
    # ONGLET 1 — CARTE MATPLOTLIB (style identique à visualize_frame)
    # ══════════════════════════════════════════════════════════
    def _build_map_tab(self, parent):
        if not MPL:
            tk.Label(parent, text='matplotlib requis.',
                     bg=BG_MAIN, fg=FG_DIM).pack(expand=True)
            return

        info = tk.Frame(parent, bg=BG_PANEL,
                        highlightthickness=1, highlightbackground=BORDER)
        info.pack(fill='x', padx=6, pady=(6, 0))
        tk.Label(info,
                 text='  Carte linéaire des deux brins — ORFs, promoteurs, terminateurs GC-riche, Shine-Dalgarno',
                 bg=BG_PANEL, fg=FG_LABEL,
                 font=('Georgia', 9)).pack(side='left', pady=6)

        fig_frame = tk.Frame(parent, bg=BG_SURFACE,
                             highlightthickness=1, highlightbackground=BORDER)
        fig_frame.pack(fill='both', expand=True, padx=6, pady=6)

        fig, (ax_p, ax_m) = plt.subplots(
            2, 1, figsize=(12, 6.5),
            facecolor=MPL_BG,
            gridspec_kw={'hspace': 0.55}
        )
        fig.suptitle(f"Carte Double Brin — {self.n} pb",
                     color=FG_MAIN, fontsize=11,
                     fontfamily='serif', y=0.98)

        self._draw_strand_ax(ax_p, self.sp)
        self._draw_strand_ax(ax_m, self.sm)

        # Légende commune en bas
        legend_patches = [
            mpatches.Patch(facecolor=C_GREEN,  label='ORF sens (+)'),
            mpatches.Patch(facecolor=C_TEAL,   label='ORF antisens (-)'),
            mpatches.Patch(facecolor=C_BLUE,   label='Promoteur Box-35', alpha=0.75),
            mpatches.Patch(facecolor=C_ORANGE, label='Promoteur Box-10', alpha=0.75),
            mpatches.Patch(facecolor=C_PURPLE, label='Shine-Dalgarno'),
            mpatches.Patch(facecolor=C_ORANGE, label='Terminateur GC-riche', alpha=0.3),
        ]
        fig.legend(handles=legend_patches,
                   loc='lower center', ncol=6,
                   facecolor=BG_SURFACE, edgecolor=BORDER,
                   labelcolor=FG_MAIN, fontsize=8,
                   prop={'family': 'serif'},
                   bbox_to_anchor=(0.5, 0.01))

        plt.tight_layout(rect=[0, 0.07, 1, 0.97])
        self._embed_fig(fig, fig_frame)

    def _draw_strand_ax(self, ax, strand):
        n        = self.n
        is_best  = strand.get('is_best', False)
        oc       = strand['orf_color']

        ax.set_facecolor(MPL_AXES)
        ax.set_xlim(0, n)
        ax.set_ylim(-0.5, 1.9)
        ax.set_xlabel('Position (pb)', color=FG_DIM, fontsize=8, fontfamily='serif')
        ax.tick_params(colors=FG_DIM, labelsize=7)
        for sp in ax.spines.values():
            sp.set_edgecolor(BORDER)
            sp.set_linewidth(0.7)
        ax.set_yticks([])

        mark = '  ★  MEILLEUR CANDIDAT' if is_best else ''
        title_color = C_GOLD if is_best else FG_MAIN
        ax.set_title(
            f"{strand['label']}{mark}  —  score {strand['score']}/100",
            color=title_color, fontsize=9, fontfamily='serif', pad=4,
            fontweight='bold' if is_best else 'normal'
        )

        mid = 0.5
        # Backbone
        ax.axhline(y=mid, color=BORDER, linewidth=2, zorder=1)
        ax.text(-3,  mid, "5'", va='center', ha='right',
                color=FG_DIM, fontsize=8, fontfamily='serif')
        ax.text(n+3, mid, "3'", va='center', ha='left',
                color=FG_DIM, fontsize=8, fontfamily='serif')

        # ── Terminateurs GC-riche (fond orange, en arrière-plan) ──
        for t in strand['terminators']:
            ax.axvspan(t['pos'], t['end'],
                       ymin=0.05, ymax=0.95,
                       facecolor=C_ORANGE, alpha=0.12, zorder=1)

        # ── Promoteurs (lignes verticales colorées) ────────────
        for p in strand['promoters'][:8]:
            ax.axvline(x=p['pos35'], color=C_BLUE,   linewidth=1.3, alpha=0.75, zorder=2)
            ax.axvline(x=p['pos10'], color=C_ORANGE, linewidth=1.3, alpha=0.75, zorder=2)

        # ── Shine-Dalgarno (triangles) ─────────────────────────
        for sd in strand['sd'][:8]:
            ax.plot(sd['pos'], 1.35, 'v',
                    color=C_PURPLE, markersize=7, alpha=0.85, zorder=3)

        # ── ORFs (rectangles FancyBbox comme visualize_frame) ──
        best_orf = strand.get('best')
        for i, seg in enumerate(strand['frames'][:12]):
            y      = mid + 0.18 if i % 2 == 0 else mid + 0.20
            is_bst = (best_orf and
                      seg['start'] == best_orf['start'] and
                      seg['end']   == best_orf['end'])
            color  = C_GOLD if is_bst else oc
            height = 0.22 if is_bst else 0.18
            alpha  = 0.95 if is_bst else 0.7

            rect = mpatches.FancyBboxPatch(
                (seg['start'], y - height/2),
                max(2, seg['end'] - seg['start']),
                height,
                boxstyle='round,pad=0.005',
                facecolor=color,
                edgecolor='#FFFFFF', linewidth=0.8,
                alpha=alpha, zorder=4
            )
            ax.add_patch(rect)

            if seg['length'] > n * 0.04:
                label = f"★ ORF{i+1}" if is_bst else f"ORF{i+1}"
                ax.text(
                    seg['start'] + seg['length']/2, y,
                    label,
                    ha='center', va='center',
                    color='#FFFFFF', fontsize=7,
                    fontweight='bold' if is_bst else 'normal',
                    fontfamily='serif', zorder=5
                )

            # Marque ATG
            if seg.get('has_atg') and seg.get('atg_pos') is not None:
                ax.plot(seg['atg_pos'], y + height/2 + 0.06,
                        'v', color='#2ECC71', markersize=5, zorder=5)

        # Graduations verticales légères
        step = max(50, (n // 5 // 50) * 50)
        for pos in range(0, n+1, step):
            ax.axvline(x=pos, color='#BDC3C7',
                       linewidth=0.4, linestyle=':', alpha=0.5, zorder=0)

    # ══════════════════════════════════════════════════════════
    # ONGLETS 2/3 — SÉQUENCE ADN COLORÉE
    # ══════════════════════════════════════════════════════════
    def _build_seq_tab(self, parent, strand):
        seq = strand['seq']
        n   = len(seq)

        # ── Construire la carte pos → couleur ─────────────────
        cmap = {}   # pos -> {'bg', 'fg', 'bold'}

        # 1) ORFs — fond vert
        best_orf = strand.get('best')
        for seg in strand['frames']:
            is_bst = (best_orf and
                      seg['start'] == best_orf['start'] and
                      seg['end']   == best_orf['end'])
            bg = SEQ_ORF_BEST if is_bst else SEQ_ORF
            for pos in range(seg['start'], min(seg['end'], n)):
                cmap.setdefault(pos, {})['bg'] = bg

        # 2) Terminateurs GC-riche — fond orange (priorité sur ORF)
        for t in strand['terminators']:
            for pos in range(t['pos'], min(t['end'], n)):
                cmap.setdefault(pos, {})['bg'] = SEQ_TERM

        # 3) Promoteurs — box35 bleu / box10 rouge
        for p in strand['promoters'][:8]:
            for pos in range(p['pos35'], min(p['pos35']+6, n)):
                cmap.setdefault(pos, {})['bg'] = SEQ_BOX35
            for pos in range(p['pos10'], min(p['pos10']+6, n)):
                cmap.setdefault(pos, {})['bg'] = SEQ_BOX10

        # 4) Shine-Dalgarno — fond violet
        for sd in strand['sd']:
            for pos in range(sd['pos'], min(sd['pos']+len(sd.get('seq','AGGAGG')), n)):
                cmap.setdefault(pos, {})['bg'] = SEQ_SD

        # 5) Codon ATG en gras vert foncé
        if best_orf and best_orf.get('has_atg') and best_orf.get('atg_pos') is not None:
            for dp in range(3):
                pos = best_orf['atg_pos'] + dp
                if 0 <= pos < n:
                    cmap.setdefault(pos, {})
                    cmap[pos]['fg']   = '#0B5345'
                    cmap[pos]['bold'] = True

        # 6) Codons stop des ORFs — texte rouge gras
        for seg in strand['frames']:
            for dp in range(3):
                pos = seg['end'] - 3 + dp
                if 0 <= pos < n:
                    cmap.setdefault(pos, {})
                    cmap[pos]['fg']   = '#922B21'
                    cmap[pos]['bold'] = True

        # ── Interface ─────────────────────────────────────────
        is_best = strand.get('is_best', False)
        hdr_bg  = C_GOLD if is_best else strand['orf_color']
        mark    = '  ★  MEILLEUR CANDIDAT' if is_best else ''

        hdr = tk.Frame(parent, bg=hdr_bg, height=40)
        hdr.pack(fill='x', padx=6, pady=(6, 0))
        hdr.pack_propagate(False)
        tk.Label(hdr,
                 text=f"  {strand['label']}{mark}  —  {n} pb",
                 bg=hdr_bg, fg='#FFFFFF',
                 font=('Georgia', 10, 'bold')).pack(side='left', pady=8)
        score_bg = '#27AE60' if strand['score'] >= 70 else '#E67E22'
        tk.Label(hdr,
                 text=f"Score : {strand['score']}/100",
                 bg=score_bg, fg='#FFFFFF',
                 font=('Georgia', 9, 'bold'),
                 padx=10).pack(side='right', pady=8, padx=10)

        # Légende
        leg = tk.Frame(parent, bg=BG_PANEL,
                        highlightthickness=1, highlightbackground=BORDER)
        leg.pack(fill='x', padx=6, pady=4)
        tk.Label(leg, text='  Légende : ', bg=BG_PANEL, fg=FG_MAIN,
                 font=('Georgia', 8, 'bold')).pack(side='left', pady=4)
        for bg, lbl in [
            (SEQ_ORF_BEST, 'ORF ★ meilleur'),
            (SEQ_ORF,      'ORF'),
            (SEQ_BOX35,    'Box -35'),
            (SEQ_BOX10,    'Box -10'),
            (SEQ_TERM,     'Term. GC-riche'),
            (SEQ_SD,       'Shine-Dalgarno'),
        ]:
            f = tk.Frame(leg, bg=bg, padx=5, pady=1,
                         highlightthickness=1, highlightbackground='#BDC3C7')
            f.pack(side='left', padx=3, pady=4)
            tk.Label(f, text=lbl, bg=bg, fg='#1A3A2A',
                     font=('Georgia', 7)).pack()

        # Zone texte scrollable
        outer = tk.Frame(parent, bg=BORDER, padx=1, pady=1)
        outer.pack(fill='both', expand=True, padx=6, pady=(0, 6))
        cv = tk.Canvas(outer, bg=BG_SURFACE, bd=0, highlightthickness=0)
        vsb = ttk.Scrollbar(outer, orient='vertical', command=cv.yview)
        cv.configure(yscrollcommand=vsb.set)
        vsb.pack(side='right', fill='y')
        cv.pack(side='left', fill='both', expand=True)
        sf = tk.Frame(cv, bg=BG_SURFACE)
        win = cv.create_window((0, 0), window=sf, anchor='nw')
        sf.bind('<Configure>', lambda e: cv.configure(scrollregion=cv.bbox('all')))
        cv.bind('<Configure>', lambda e: cv.itemconfig(win, width=e.width))
        cv.bind_all('<MouseWheel>',
                    lambda e: cv.yview_scroll(int(-1*(e.delta/120)), 'units'))

        self._render_seq(sf, seq, cmap, strand)

    def _render_seq(self, parent, seq, cmap, strand):
        N_PER_LINE = 60

        # Info-ligne ORF
        info = tk.Frame(parent, bg=BG_SURFACE, padx=10, pady=6)
        info.pack(fill='x')
        best = strand.get('best')
        if best:
            tk.Label(info,
                     text=(f"ORF : Frame {best['frame']:+d}  |  "
                           f"pos. {best['start']+1}–{best['end']}  |  "
                           f"{best['length']} pb  |  {best['num_aa']} aa  |  "
                           f"ATG {'✓' if best.get('has_atg') else '✗'}  |  "
                           f"Stop : {best.get('stop_codon','?')}"),
                     bg=BG_SURFACE, fg=ACCENT,
                     font=('Courier New', 9, 'bold')).pack(anchor='w')
        tk.Label(info,
                 text=(f"{len(strand['promoters'])} promoteur(s)  •  "
                       f"{len(strand['terminators'])} terminateur(s) GC-riche  •  "
                       f"{len(strand['sd'])} site(s) SD"),
                 bg=BG_SURFACE, fg=FG_DIM,
                 font=('Georgia', 8, 'italic')).pack(anchor='w')

        tk.Frame(parent, bg=BORDER, height=1).pack(fill='x', padx=8, pady=4)

        txt = tk.Frame(parent, bg=BG_SURFACE, padx=8)
        txt.pack(fill='both', expand=True)

        n = len(seq)
        for line_s in range(0, n, N_PER_LINE):
            line_e = min(line_s + N_PER_LINE, n)
            row = tk.Frame(txt, bg=BG_SURFACE)
            row.pack(anchor='w', pady=1)

            # Numéro de position
            tk.Label(row,
                     text=f"{line_s+1:>5} ",
                     bg=BG_SURFACE, fg='#95A5A6',
                     font=('Courier New', 9)).pack(side='left')

            # Bases groupées par couleur identique
            pos = line_s
            while pos < line_e:
                cm   = cmap.get(pos, {})
                bg   = cm.get('bg', BG_SURFACE)
                fg   = cm.get('fg', '#1B3A2E')
                bold = cm.get('bold', False)

                # Grouper jusqu'à 3 bases consécutives de même couleur
                end = pos + 1
                while (end < line_e and end - pos < 3 and
                       not cmap.get(end, {}).get('bold', False) and
                       cmap.get(end, {}).get('bg', BG_SURFACE) == bg):
                    end += 1

                bases = seq[pos:end]
                font  = ('Courier New', 10, 'bold') if bold else ('Courier New', 10)
                tk.Label(row, text=bases, bg=bg, fg=fg,
                         font=font, padx=0, pady=1, bd=0
                         ).pack(side='left')

                # Espace inter-codon
                if end % 3 == 0 and end < line_e:
                    tk.Label(row, text=' ', bg=BG_SURFACE,
                             font=('Courier New', 10)).pack(side='left')
                pos = end

            # Annotations à droite
            annots = []
            for p in strand['promoters'][:8]:
                if line_s <= p['pos35'] < line_e:
                    annots.append(f"─35@{p['pos35']+1}")
                if line_s <= p['pos10'] < line_e:
                    annots.append(f"─10@{p['pos10']+1}")
            for t in strand['terminators']:
                if line_s <= t['pos'] < line_e:
                    s2 = t['arm1'] + t['arm2']
                    gc = (s2.count('G')+s2.count('C'))/max(1,len(s2))*100
                    annots.append(f"TERM@{t['pos']+1}(GC:{gc:.0f}%)")
            if annots:
                tk.Label(row,
                         text='  ← ' + '  '.join(annots),
                         bg=BG_SURFACE, fg=FG_DIM,
                         font=('Courier New', 7, 'italic')).pack(side='left')

    # ══════════════════════════════════════════════════════════
    # ONGLET 4 — RÉSULTAT OPTIMAL
    # ══════════════════════════════════════════════════════════
    def _build_optimal_tab(self, parent):
        # Scroll wrapper
        container = tk.Frame(parent, bg=BG_MAIN)
        container.pack(fill='both', expand=True)
        cv  = tk.Canvas(container, bg=BG_MAIN, bd=0, highlightthickness=0)
        vsb = ttk.Scrollbar(container, orient='vertical', command=cv.yview)
        cv.configure(yscrollcommand=vsb.set)
        vsb.pack(side='right', fill='y')
        cv.pack(side='left', fill='both', expand=True)
        inner = tk.Frame(cv, bg=BG_MAIN)
        wid   = cv.create_window((0, 0), window=inner, anchor='nw')
        inner.bind('<Configure>', lambda e: cv.configure(scrollregion=cv.bbox('all')))
        cv.bind('<Configure>', lambda e: cv.itemconfig(wid, width=e.width))
        cv.bind_all('<MouseWheel>',
                    lambda e: cv.yview_scroll(int(-1*(e.delta/120)), 'units'))

        best = self.best
        orf  = best.get('best')

        # ── Bannière or ────────────────────────────────────────
        banner = tk.Frame(inner, bg='#FEF9E7',
                           highlightthickness=2,
                           highlightbackground=C_GOLD)
        banner.pack(fill='x', padx=8, pady=(8, 6))
        tk.Label(banner,
                 text='★  RÉSULTAT OPTIMAL — ANALYSE DOUBLE BRIN',
                 bg='#FEF9E7', fg='#7D6608',
                 font=('Georgia', 13, 'bold'),
                 padx=16, pady=10).pack(anchor='w')

        cards = tk.Frame(banner, bg='#FEF9E7')
        cards.pack(fill='x', padx=12, pady=(0, 8))
        infos = [
            ('Brin optimal', best['label'],
             C_GREEN if best['sign']=='+' else C_TEAL),
            ('Score codant', f"{best['score']}/100",
             '#27AE60' if best['score']>=70 else '#E67E22'),
            ('Cadre lecture', f"Frame {orf['frame']:+d}" if orf else '—', ACCENT),
            ('Position ORF',
             f"{orf['start']+1} → {orf['end']}" if orf else '—', ACCENT),
            ('Longueur ADN',
             f"{orf['length']} pb" if orf else '—', ACCENT),
            ('Taille protéine',
             f"{orf['num_aa']} aa" if orf else '—', ACCENT),
            ('Codon start',
             '✓ ATG' if orf and orf.get('has_atg') else '✗ Absent',
             '#27AE60' if orf and orf.get('has_atg') else '#C0392B'),
            ('Codon stop',
             orf.get('stop_codon','—') if orf else '—', ACCENT),
            ('Promoteurs',     str(len(best['promoters'])), C_BLUE),
            ('Term. GC-riche', str(len(best['terminators'])), C_ORANGE),
            ('Sites SD',       str(len(best['sd'])), C_PURPLE),
            ('Masse estimée',
             f"≈ {round(orf['num_aa']*110/1000,1)} kDa" if orf else '—', ACCENT),
        ]
        for j, (lbl, val, color) in enumerate(infos):
            col = j % 4
            row = j // 4
            f = tk.Frame(cards, bg='#FFFDE7',
                         highlightthickness=1,
                         highlightbackground='#F9E79F',
                         padx=10, pady=8)
            f.grid(row=row, column=col, padx=4, pady=4, sticky='ew')
            tk.Label(f, text=lbl, bg='#FFFDE7', fg='#7D6608',
                     font=('Georgia', 8)).pack(anchor='w')
            tk.Label(f, text=val, bg='#FFFDE7', fg=color,
                     font=('Courier New', 10, 'bold'),
                     wraplength=160).pack(anchor='w')
            cards.columnconfigure(col, weight=1)

        # Justification
        if best['scored']:
            jf = tk.Frame(inner, bg=BG_SURFACE,
                          highlightthickness=1, highlightbackground=BORDER)
            jf.pack(fill='x', padx=8, pady=4)
            tk.Label(jf, text='  Critères de sélection du brin optimal :',
                     bg=BG_SURFACE, fg=FG_MAIN,
                     font=('Georgia', 9, 'bold'), pady=6).pack(anchor='w')
            for note in best['scored'][0][1]:
                tk.Label(jf, text=f'    ✓  {note}',
                         bg=BG_SURFACE, fg='#145C38',
                         font=('Courier New', 9)).pack(anchor='w', pady=1)
            tk.Label(jf, text='', bg=BG_SURFACE).pack(pady=2)

        # ── Mini-carte matplotlib du brin optimal ──────────────
        if MPL:
            map_hdr = tk.Frame(inner, bg=ACCENT, height=32)
            map_hdr.pack(fill='x', padx=8, pady=(8, 0))
            map_hdr.pack_propagate(False)
            tk.Label(map_hdr,
                     text=f'  Carte du brin optimal : {best["label"]}',
                     bg=ACCENT, fg='#FFFFFF',
                     font=('Georgia', 9, 'bold')).pack(side='left', pady=4)

            mf = tk.Frame(inner, bg=BG_SURFACE,
                           highlightthickness=1, highlightbackground=BORDER)
            mf.pack(fill='x', padx=8, pady=(0, 6))

            fig2, ax2 = plt.subplots(figsize=(11, 2.8), facecolor=MPL_BG)
            self._draw_strand_ax(ax2, best)
            # Ajouter mini-légende dans cet axe
            legend_patches = [
                mpatches.Patch(facecolor=best['orf_color'], label='ORF'),
                mpatches.Patch(facecolor=C_GOLD,   label='ORF ★ meilleur'),
                mpatches.Patch(facecolor=C_BLUE,   label='Box-35'),
                mpatches.Patch(facecolor=C_ORANGE, label='Box-10 / Terminateur'),
                mpatches.Patch(facecolor=C_PURPLE, label='Shine-Dalgarno'),
            ]
            ax2.legend(handles=legend_patches, loc='upper right',
                       facecolor=BG_SURFACE, edgecolor=BORDER,
                       labelcolor=FG_MAIN, fontsize=7,
                       prop={'family': 'serif'}, framealpha=0.95)
            plt.tight_layout()
            self._embed_fig(fig2, mf, toolbar=False)

        # ── Protéine colorée ───────────────────────────────────
        if orf:
            ph = tk.Frame(inner, bg=ACCENT, height=32)
            ph.pack(fill='x', padx=8, pady=(8, 0))
            ph.pack_propagate(False)
            tk.Label(ph,
                     text=f"  Protéine traduite — {orf['num_aa']} acides aminés",
                     bg=ACCENT, fg='#FFFFFF',
                     font=('Georgia', 9, 'bold')).pack(side='left', pady=4)

            pf = tk.Frame(inner, bg=BORDER, padx=1, pady=1)
            pf.pack(fill='x', padx=8, pady=(0, 6))
            tw = tk.Text(pf, height=4, bg=BG_SURFACE, fg=FG_MAIN,
                         font=('Courier New', 11), wrap='word',
                         bd=0, padx=10, pady=8, state='normal')
            tw.pack(fill='x')
            tw.tag_configure('hydro',  foreground='#E67E22',
                             font=('Courier New', 11, 'bold'))
            tw.tag_configure('polar',  foreground='#2196A6')
            tw.tag_configure('charge', foreground='#8E44AD')
            tw.tag_configure('other',  foreground='#27AE60')
            from data.codon_table import AA_HYDROPHOBIC, AA_POLAR, AA_CHARGED
            for aa in orf.get('protein', ''):
                tag = ('hydro'  if aa in AA_HYDROPHOBIC else
                       'polar'  if aa in AA_POLAR       else
                       'charge' if aa in AA_CHARGED      else 'other')
                tw.insert('end', aa, tag)
            tw.config(state='disabled')

            leg_aa = tk.Frame(inner, bg=BG_SURFACE)
            leg_aa.pack(fill='x', padx=12, pady=(0, 8))
            for col, lbl in [('#E67E22','Hydrophobe'), ('#2196A6','Polaire'),
                              ('#8E44AD','Chargé'),    ('#27AE60','Autre')]:
                tk.Label(leg_aa, text=f'■ {lbl}', bg=BG_SURFACE, fg=col,
                         font=('Georgia', 8)).pack(side='left', padx=8)

        # ── Tableau comparaison des deux brins ─────────────────
        ch = tk.Frame(inner, bg='#5D6D7E', height=32)
        ch.pack(fill='x', padx=8, pady=(8, 0))
        ch.pack_propagate(False)
        tk.Label(ch, text='  Comparaison brin (+) vs brin (−)',
                 bg='#5D6D7E', fg='#FFFFFF',
                 font=('Georgia', 9, 'bold')).pack(side='left', pady=4)

        comp = tk.Frame(inner, bg=BG_SURFACE,
                        highlightthickness=1, highlightbackground=BORDER)
        comp.pack(fill='x', padx=8, pady=(0, 12))
        comp.columnconfigure(1, weight=1)
        comp.columnconfigure(2, weight=1)

        headers = ['Critère', 'Brin sens (+)', 'Brin antisens (-)']
        for c, (h, fg_c) in enumerate(zip(headers,
                                           [FG_MAIN, C_GREEN, C_TEAL])):
            tk.Label(comp, text=h, bg=BG_PANEL, fg=fg_c,
                     font=('Georgia', 9, 'bold'), padx=12, pady=6
                     ).grid(row=0, column=c, sticky='ew', padx=2, pady=2)

        rows_data = [
            ('Score codant',
             f"{self.sp['score']}/100", f"{self.sm['score']}/100"),
            ('Nb ORFs détectés',
             str(len(self.sp['frames'])), str(len(self.sm['frames']))),
            ('Meilleur ORF — longueur',
             f"{self.sp['best']['length']} pb" if self.sp['best'] else '—',
             f"{self.sm['best']['length']} pb" if self.sm['best'] else '—'),
            ('Codon ATG',
             '✓' if self.sp['best'] and self.sp['best'].get('has_atg') else '✗',
             '✓' if self.sm['best'] and self.sm['best'].get('has_atg') else '✗'),
            ('Promoteurs détectés',
             str(len(self.sp['promoters'])), str(len(self.sm['promoters']))),
            ('Terminateurs GC-riche',
             str(len(self.sp['terminators'])), str(len(self.sm['terminators']))),
            ('Sites Shine-Dalgarno',
             str(len(self.sp['sd'])), str(len(self.sm['sd']))),
        ]
        for r, (label, vp, vm) in enumerate(rows_data, 1):
            bg_r  = BG_SURFACE if r % 2 == 0 else BG_PANEL
            win_p = (str(vp) > str(vm)) or (vp == '✓' and vm != '✓')
            win_m = (str(vm) > str(vp)) or (vm == '✓' and vp != '✓')
            fp    = '#C8E6C9' if win_p else bg_r
            fm    = '#C8E6C9' if win_m else bg_r

            tk.Label(comp, text=label, bg=bg_r, fg=FG_MAIN,
                     font=('Georgia', 9), padx=12, pady=5
                     ).grid(row=r, column=0, sticky='ew', padx=2, pady=1)
            tk.Label(comp, text=vp, bg=fp,
                     fg=C_GREEN if win_p else FG_MAIN,
                     font=('Courier New', 9, 'bold' if win_p else 'normal'),
                     padx=12, pady=5
                     ).grid(row=r, column=1, sticky='ew', padx=2, pady=1)
            tk.Label(comp, text=vm, bg=fm,
                     fg=C_TEAL if win_m else FG_MAIN,
                     font=('Courier New', 9, 'bold' if win_m else 'normal'),
                     padx=12, pady=5
                     ).grid(row=r, column=2, sticky='ew', padx=2, pady=1)

    # ══════════════════════════════════════════════════════════
    # UTILITAIRE — embed matplotlib
    # ══════════════════════════════════════════════════════════
    def _embed_fig(self, fig, frame, toolbar=True):
        canvas = FigureCanvasTkAgg(fig, master=frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)
        if toolbar:
            tb_frame = tk.Frame(frame, bg=BG_PANEL,
                                highlightthickness=1,
                                highlightbackground=BORDER)
            tb_frame.pack(fill='x')
            tb = NavigationToolbar2Tk(canvas, tb_frame)
            tb.config(bg=BG_PANEL)
            for child in tb.winfo_children():
                try:
                    child.config(bg=BG_PANEL)
                except Exception:
                    pass
            tb.update()