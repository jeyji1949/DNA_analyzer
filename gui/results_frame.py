# ============================================================
# Panneau de résultats — Thème biologique clair
# ============================================================

import tkinter as tk
from tkinter import ttk
import math

from analysis.orf_finder import reverse_complement

BG_MAIN    = '#F0F9F4'
BG_SURFACE = '#FFFFFF'
BG_PANEL   = '#E8F5EE'
BG_CARD    = '#F7FBF8'
BORDER     = '#B8DFC9'
FG_MAIN    = '#1A3A2A'
FG_DIM     = '#6B9E80'
FG_LABEL   = '#2D6A4F'
ACCENT     = '#1B7A4A'
ACCENT2    = '#2196A6'
ACCENT3    = '#E67E22'
ACCENT4    = '#8E44AD'
GREEN_S    = '#27AE60'
RED_S      = '#C0392B'
YELLOW_S   = '#F39C12'

# Couleurs de tag (pastel)
TAG_GREEN_BG  = '#D4EDDA'
TAG_BLUE_BG   = '#D1ECF1'
TAG_ORANGE_BG = '#FDEBD0'
TAG_PURPLE_BG = '#E8DAEF'
TAG_RED_BG    = '#FADBD8'


def _badge(parent, text, bg, fg=None):
    fg = fg or FG_MAIN
    lbl = tk.Label(parent, text=f' {text} ', bg=bg, fg=fg,
                   font=('Georgia', 8, 'bold'),
                   relief='flat', padx=4, pady=1)
    return lbl


class ResultsFrame(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent, bg=BG_MAIN)
        self._build()

    def _build(self):
        # Canvas de fond avec motifs ADN subtils
        self.bg_canvas = tk.Canvas(self, bg=BG_MAIN, bd=0, highlightthickness=0)
        self.bg_canvas.place(x=0, y=0, relwidth=1.0, relheight=1.0)
        self.bg_canvas.bind('<Configure>', lambda e: self._draw_bg())

        self.inner_nb = ttk.Notebook(self)
        self.inner_nb.pack(fill='both', expand=True, padx=6, pady=6)

        self.tab_orfs        = self._tab('🔍  ORFs')
        self.tab_protein     = self._tab('⚙️  Protéine')
        self.tab_promoters   = self._tab('📍  Promoteurs')
        self.tab_sd          = self._tab('🎯  Shine-Dalgarno')
        self.tab_terminators = self._tab('🔚  Terminateurs')
        self.tab_restriction = self._tab('✂️  Restriction')
        self.tab_stats       = self._tab('📊  Statistiques')

    def _tab(self, title):
        f = tk.Frame(self.inner_nb, bg=BG_SURFACE)
        self.inner_nb.add(f, text=f'  {title}  ')
        return f

    def _draw_bg(self):
        c = self.bg_canvas
        c.delete('bgpat')
        w = c.winfo_width()
        h = c.winfo_height()
        if w < 10:
            return
        # Hélice décorative en arrière-plan (très pâle)
        step = 20
        amp  = 14
        freq = 0.035
        for strand_offset in [0, math.pi]:
            pts = []
            for x in range(0, w+step, step):
                y = h*0.5 + amp * math.sin(freq*x + strand_offset)
                pts.append((x, y))
            for i in range(len(pts)-1):
                c.create_line(pts[i][0], pts[i][1], pts[i+1][0], pts[i+1][1],
                               fill='#D8F3DC', width=3, smooth=True, tags='bgpat')
        # Petits cercles (bases)
        for i, x in enumerate(range(0, w, 40)):
            y1 = h*0.5 + amp * math.sin(freq*x)
            y2 = h*0.5 + amp * math.sin(freq*x + math.pi)
            c.create_line(x, y1, x, y2, fill='#C8E6C9', width=1.5, tags='bgpat')

    def display(self, results):
        self._draw_bg()
        self._show_orfs(results)
        self._show_protein(results)
        self._show_promoters(results)
        self._show_sd(results)
        self._show_terminators(results)
        self._show_restriction(results)
        self._show_stats(results)

    def _clear(self, tab):
        for w in tab.winfo_children():
            w.destroy()

    def _section_header(self, parent, text, count=None, ok=True):
        f = tk.Frame(parent, bg=BG_PANEL,
                     highlightthickness=1, highlightbackground=BORDER)
        f.pack(fill='x', padx=8, pady=(8,4))
        color = ACCENT if ok else FG_DIM
        tk.Label(f, text=text, bg=BG_PANEL, fg=color,
                 font=('Georgia', 11, 'bold'),
                 padx=12, pady=6).pack(side='left')
        if count is not None:
            badge_bg = TAG_GREEN_BG if ok and count > 0 else TAG_RED_BG if count == 0 else TAG_BLUE_BG
            badge_fg = GREEN_S if count > 0 else RED_S
            _badge(f, str(count), badge_bg, badge_fg).pack(side='left', padx=4)
        return f

    def _info_hint(self, parent, text):
        tk.Label(parent, text=text, bg=BG_SURFACE, fg=FG_DIM,
                 font=('Georgia', 8, 'italic'),
                 padx=10, pady=2).pack(anchor='w', padx=8)

    def _make_tree(self, parent, columns, widths=None):
        frame = tk.Frame(parent, bg=BG_SURFACE,
                         highlightthickness=1, highlightbackground=BORDER)
        frame.pack(fill='both', expand=True, padx=8, pady=4)
        tree = ttk.Treeview(frame, columns=columns, show='headings', height=14)
        for i, col in enumerate(columns):
            tree.heading(col, text=col)
            w = widths[i] if widths else 100
            tree.column(col, width=w, anchor='center', minwidth=40)
        vsb = ttk.Scrollbar(frame, orient='vertical',   command=tree.yview)
        hsb = ttk.Scrollbar(frame, orient='horizontal', command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.grid(row=0, column=0, sticky='nsew')
        vsb.grid(row=0, column=1, sticky='ns')
        hsb.grid(row=1, column=0, sticky='ew')
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        return tree

    def _empty(self, parent, msg):
        tk.Label(parent, text=msg, bg=BG_SURFACE, fg=FG_DIM,
                 font=('Georgia', 10)).pack(expand=True, pady=40)

    # ── ORFs ──────────────────────────────────────────────────
    def _show_orfs(self, results):
        tab = self.tab_orfs
        self._clear(tab)
        orfs = results.get('orfs', [])
        self._section_header(tab, '  Open Reading Frames (ORFs)', len(orfs), bool(orfs))
        self._info_hint(tab, 'Détection sur les 6 cadres de lecture • longueur min. réglable dans le panneau gauche')
        if not orfs:
            msg = ('Aucun ORF détecté.\n\n'
                   '→ Réduisez la longueur minimale dans le panneau gauche\n'
                   '→ La séquence est peut-être trop courte\n'
                   '→ Vérifiez que la séquence contient bien ATGC')
            self._empty(tab, msg)
            return
        cols   = ('#', 'Frame', 'Début', 'Fin', 'Longueur (pb)', 'Nb AA', 'Brin', 'Codon Stop', 'Statut')
        widths = (40,  70,      80,      80,    120,            70,      60,     100,           110)
        tree = self._make_tree(tab, cols, widths)
        for i, orf in enumerate(orfs):
            status = '★ Meilleur candidat' if i == 0 and orf['strand'] == '+' else 'Candidat'
            tag    = 'best' if i == 0 and orf['strand'] == '+' else \
                     'plus' if orf['strand'] == '+' else 'minus'
            tree.insert('', 'end', values=(
                i+1, f"{'+' if orf['strand']=='+' else ''}{orf['frame']}",
                orf['start']+1, orf['end'], orf['length'],
                orf['num_aa'], orf['strand'], orf['stop_codon'], status,
            ), tags=(tag,))
        tree.tag_configure('best',  background='#C8E6C9', foreground='#145C38', font=('Courier New', 9, 'bold'))
        tree.tag_configure('plus',  foreground=ACCENT)
        tree.tag_configure('minus', foreground=ACCENT2)

    # ── Protéine ──────────────────────────────────────────────
    def _show_protein(self, results):
        tab = self.tab_protein
        self._clear(tab)
        orfs = [o for o in results.get('orfs', []) if o['strand'] == '+']
        if not orfs:
            self._section_header(tab, '  Traduction Protéique', 0, False)
            msg = ('Aucun ORF détecté sur le brin sens (+).\n\n'
                   '→ Des ORFs peuvent exister sur le brin antisens (-)\n'
                   '→ Consultez l\'onglet ORFs pour voir tous les cadres\n'
                   '→ Réduisez la longueur minimale si nécessaire')
            self._empty(tab, msg)
            return
        best = orfs[0]
        protein = best['protein']
        self._section_header(tab, f"  ORF codant — Frame +{best['frame']} — {len(protein)} aa", None, True)

        # Cartes d'info
        info_frame = tk.Frame(tab, bg=BG_SURFACE)
        info_frame.pack(fill='x', padx=8, pady=4)
        infos = [
            ('Position', f"{best['start']+1} → {best['end']}"),
            ('Longueur ADN', f"{best['length']} pb"),
            ('Protéine', f"{len(protein)} aa"),
            ('Poids mol.', f"≈ {len(protein)*110/1000:.1f} kDa"),
            ('Codon Start', 'ATG'),
            ('Codon Stop', best['stop_codon']),
        ]
        for j, (label, val) in enumerate(infos):
            col = j % 3
            row = j // 3
            card = tk.Frame(info_frame, bg=BG_CARD,
                            highlightthickness=1, highlightbackground=BORDER,
                            padx=10, pady=8)
            card.grid(row=row, column=col, padx=4, pady=4, sticky='ew')
            tk.Label(card, text=label, bg=BG_CARD, fg=FG_DIM,
                     font=('Georgia', 8)).pack(anchor='w')
            tk.Label(card, text=val, bg=BG_CARD, fg=ACCENT,
                     font=('Courier New', 11, 'bold')).pack(anchor='w')
        for col in range(3):
            info_frame.columnconfigure(col, weight=1)

        # Séquence protéique
        tk.Label(tab, text='Séquence protéique traduite (colorée par propriété) :',
                 bg=BG_SURFACE, fg=FG_LABEL,
                 font=('Georgia', 9)).pack(anchor='w', padx=10, pady=(6, 2))

        pf = tk.Frame(tab, bg=BORDER, padx=1, pady=1)
        pf.pack(fill='both', expand=True, padx=8, pady=(0, 8))
        tw = tk.Text(pf, bg=BG_SURFACE, fg=FG_MAIN,
                     font=('Courier New', 10), wrap='word',
                     state='normal', bd=0, padx=8, pady=8)
        tw.pack(fill='both', expand=True)
        tw.tag_configure('hydro',  foreground='#E67E22', font=('Courier New', 10, 'bold'))
        tw.tag_configure('polar',  foreground='#2196A6')
        tw.tag_configure('charge', foreground='#8E44AD')
        tw.tag_configure('other',  foreground='#27AE60')
        tw.tag_configure('stop',   foreground='#C0392B', font=('Courier New', 10, 'bold'))
        from data.codon_table import AA_HYDROPHOBIC, AA_POLAR, AA_CHARGED
        for aa in protein:
            tag = ('hydro'  if aa in AA_HYDROPHOBIC else
                   'polar'  if aa in AA_POLAR       else
                   'charge' if aa in AA_CHARGED      else
                   'stop'   if aa == '*'              else
                   'other')
            tw.insert('end', aa, tag)
        tw.config(state='disabled')

        # Légende
        leg = tk.Frame(tab, bg=BG_SURFACE)
        leg.pack(fill='x', padx=10, pady=4)
        for col, lbl in [('#E67E22','Hydrophobe'), ('#2196A6','Polaire'),
                          ('#8E44AD','Chargé'),    ('#27AE60','Autre/spécial')]:
            tk.Label(leg, text=f'■ {lbl}', bg=BG_SURFACE, fg=col,
                     font=('Georgia', 8)).pack(side='left', padx=8)

    # ── Promoteurs ────────────────────────────────────────────
    def _show_promoters(self, results):
        tab = self.tab_promoters
        self._clear(tab)
        promoters = results.get('promoters', [])
        self._section_header(tab, '  Promoteurs Bactériens (-10 / -35)', len(promoters), bool(promoters))
        self._info_hint(tab, 'Box -35: TTGACA  •  Box -10: TATAAT  •  Espacement optimal: 17 pb (14–22 accepté)')
        if not promoters:
            self._empty(tab, 'Aucun promoteur complet détecté.')
            return
        cols   = ('#', 'Pos -35', 'Seq -35', 'MM', 'Pos -10', 'Seq -10', 'MM', 'Espacement', 'Qualité')
        widths = (40,  90,        90,         50,   90,        90,         50,   100,           130)
        tree = self._make_tree(tab, cols, widths)
        for i, p in enumerate(promoters):
            tag = {'Consensus parfait':'best','Fort':'good','Modéré':'mod','Faible':'weak'}.get(p['quality'],'weak')
            tree.insert('', 'end', values=(
                i+1, p['pos35']+1, p['seq35'], p['dist35'],
                p['pos10']+1, p['seq10'], p['dist10'],
                p['spacing'], p['quality'],
            ), tags=(tag,))
        tree.tag_configure('best', background='#C8E6C9', foreground='#145C38', font=('Courier New', 9,'bold'))
        tree.tag_configure('good', foreground=GREEN_S)
        tree.tag_configure('mod',  foreground=YELLOW_S)
        tree.tag_configure('weak', foreground=FG_DIM)

    # ── Shine-Dalgarno ────────────────────────────────────────
    def _show_sd(self, results):
        tab = self.tab_sd
        self._clear(tab)
        sd_sites = results.get('sd_sites', [])
        self._section_header(tab, '  Sites Shine-Dalgarno', len(sd_sites), bool(sd_sites))
        self._info_hint(tab, 'Consensus: AGGAGG  •  Localisé 7–9 pb en amont de l\'ATG')
        if not sd_sites:
            self._empty(tab, 'Aucun site Shine-Dalgarno détecté.')
            return
        cols   = ('#', 'Position SD', 'Séquence SD', 'Position ATG', 'Espacement (pb)', 'Mismatches', 'Qualité')
        widths = (40,  110,           110,            110,             140,               100,           110)
        tree = self._make_tree(tab, cols, widths)
        for i, sd in enumerate(sd_sites):
            tag = 'best' if sd['dist']==0 else 'good' if sd['dist']==1 else 'mod'
            tree.insert('', 'end', values=(
                i+1, sd['pos']+1, sd['seq'], sd['atg_pos']+1,
                sd['spacing'], sd['dist'], sd['quality'],
            ), tags=(tag,))
        tree.tag_configure('best', background=TAG_PURPLE_BG, foreground='#5B2C6F', font=('Courier New',9,'bold'))
        tree.tag_configure('good', foreground=ACCENT4)
        tree.tag_configure('mod',  foreground=FG_DIM)

    # ── Terminateurs ─────────────────────────────────────────
    def _show_terminators(self, results):
        tab = self.tab_terminators
        self._clear(tab)
        terminators = results.get('terminators', [])
        self._section_header(tab, '  Terminateurs Rho-Indépendants', len(terminators), bool(terminators))
        self._info_hint(tab, 'Structure : tige-boucle palindromique + séquence polyT (≥ 3 thymine)')
        if not terminators:
            self._empty(tab, 'Aucun terminateur Rho-indépendant détecté.')
            return
        cols   = ('#', 'Position', 'Fin', 'Bras 1', 'Boucle', 'Bras 2', 'PolyT', 'Long. tige')
        widths = (40,  90,         80,    100,      80,       100,      90,       110)
        tree = self._make_tree(tab, cols, widths)
        for i, t in enumerate(terminators):
            tree.insert('', 'end', values=(
                i+1, t['pos']+1, t['end'],
                t['arm1'], t['loop'], t['arm2'],
                t['poly_t'], t['stem_len'],
            ))

    # ── Restriction ───────────────────────────────────────────
    def _show_restriction(self, results):
        tab = self.tab_restriction
        self._clear(tab)
        restriction = results.get('restriction', {})
        self._section_header(tab, '  Sites de Restriction', len(restriction), bool(restriction))
        self._info_hint(tab, 'Enzymes de restriction courantes — positions 1-based')
        if not restriction:
            self._empty(tab, 'Aucun site de restriction détecté pour les enzymes standard.')
            return
        cols   = ('Enzyme', 'Site de reconnaissance', 'Nb coupures', 'Positions (max 20)')
        widths = (120,       180,                       120,           500)
        tree = self._make_tree(tab, cols, widths)
        for enzyme, info in sorted(restriction.items(), key=lambda x: -x[1]['count']):
            pos_str = ', '.join(map(str, info['positions'][:20]))
            if len(info['positions']) > 20:
                pos_str += f'  … (+{len(info["positions"])-20})'
            tree.insert('', 'end', values=(enzyme, info['site'], info['count'], pos_str))

    # ── Statistiques ─────────────────────────────────────────
    def _show_stats(self, results):
        tab = self.tab_stats
        self._clear(tab)
        stats = results.get('stats', {})
        if not stats:
            self._empty(tab, 'Aucune statistique disponible.')
            return
        self._section_header(tab, '  Statistiques Physico-Chimiques', None, True)

        info_frame = tk.Frame(tab, bg=BG_SURFACE)
        info_frame.pack(fill='x', padx=8, pady=4)
        n    = stats.get('length', 0)
        pcts = stats.get('pcts', {})
        infos = [
            ('Longueur',        f"{n} pb"),
            ('Teneur GC',       f"{stats.get('gc', 0):.2f}%"),
            ('Teneur AT',       f"{stats.get('at', 0):.2f}%"),
            ('Ratio G/C',       f"{stats.get('ratio_gc', 0):.3f}"),
            ('Tm (Wallace)',     f"{stats.get('tm_wallace', 0)} °C"),
            ('Tm (Nearest-Nb)', f"{stats.get('tm_nn', 0)} °C"),
            ('Masse ADN ss',    f"{stats.get('mw_ss', 0)} kDa"),
            ('Masse ADN ds',    f"{stats.get('mw_ds', 0)} kDa"),
            ('A', f"{stats.get('counts',{}).get('A',0)}  ({pcts.get('A',0):.1f}%)"),
            ('T', f"{stats.get('counts',{}).get('T',0)}  ({pcts.get('T',0):.1f}%)"),
            ('G', f"{stats.get('counts',{}).get('G',0)}  ({pcts.get('G',0):.1f}%)"),
            ('C', f"{stats.get('counts',{}).get('C',0)}  ({pcts.get('C',0):.1f}%)"),
        ]
        for j, (label, val) in enumerate(infos):
            col = j % 4
            row = j // 4
            card = tk.Frame(info_frame, bg=BG_CARD,
                            highlightthickness=1, highlightbackground=BORDER,
                            padx=10, pady=8)
            card.grid(row=row, column=col, padx=4, pady=4, sticky='ew')
            tk.Label(card, text=label, bg=BG_CARD, fg=FG_DIM, font=('Georgia', 8)).pack(anchor='w')
            tk.Label(card, text=val,   bg=BG_CARD, fg=ACCENT, font=('Courier New', 11, 'bold')).pack(anchor='w')
        for col in range(4):
            info_frame.columnconfigure(col, weight=1)

        # Barre de composition
        bar_lbl = tk.Label(tab, text='Composition nucléotidique :',
                            bg=BG_SURFACE, fg=FG_LABEL, font=('Georgia', 9))
        bar_lbl.pack(anchor='w', padx=10, pady=(8, 2))

        bar_outer = tk.Frame(tab, bg=BORDER, padx=1, pady=1, height=28)
        bar_outer.pack(fill='x', padx=10)
        bar_outer.pack_propagate(False)
        bar_canvas = tk.Canvas(bar_outer, bg=BG_SURFACE, height=26, bd=0, highlightthickness=0)
        bar_canvas.pack(fill='x')

        bar_colors = {'G': '#52B788', 'C': '#48CAE4', 'A': '#F4A261', 'T': '#A29BFE'}

        def draw_bars(event=None, canvas=bar_canvas):
            canvas.delete('all')
            w = canvas.winfo_width()
            if w < 10:
                return
            x = 0
            for base, color in bar_colors.items():
                pct = pcts.get(base, 0)
                bw = int(pct / 100 * w)
                canvas.create_rectangle(x, 0, x+bw, 26, fill=color, outline='')
                if bw > 35:
                    canvas.create_text(x+bw//2, 13, text=f'{base} {pct:.1f}%',
                                        fill='#FFFFFF', font=('Georgia', 8, 'bold'))
                x += bw

        bar_canvas.bind('<Configure>', draw_bars)

        # Légende couleurs bases
        leg = tk.Frame(tab, bg=BG_SURFACE)
        leg.pack(fill='x', padx=10, pady=6)
        for base, color in bar_colors.items():
            tk.Label(leg, text=f'■ {base}', bg=BG_SURFACE, fg=color,
                     font=('Georgia', 9, 'bold')).pack(side='left', padx=10)

def display_reading_frames(seq):
    seq = seq.upper().replace(" ", "").replace("\n", "")
    rc = reverse_complement(seq)

    print("\n===== BRIN 5' → 3' =====\n")

    for frame in range(3):
        codons = [seq[i:i+3] for i in range(frame, len(seq)-2, 3)]
        print(f"Frame +{frame+1}:")
        print(" | ".join(codons))
        print()

    print("\n===== BRIN 3' → 5' (reverse complement) =====\n")

    for frame in range(3):
        codons = [rc[i:i+3] for i in range(frame, len(rc)-2, 3)]
        print(f"Frame -{frame+1}:")
        print(" | ".join(codons))
        print()