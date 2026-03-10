# ============================================================
# Panneau de visualisation — Thème biologique clair
# ============================================================

import tkinter as tk
from tkinter import ttk
import math

try:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    MATPLOTLIB_OK = True
except ImportError:
    MATPLOTLIB_OK = False

BG_MAIN    = '#F0F9F4'
BG_SURFACE = '#FFFFFF'
BG_PANEL   = '#E8F5EE'
BORDER     = '#B8DFC9'
FG_MAIN    = '#1A3A2A'
FG_DIM     = '#6B9E80'
FG_LABEL   = '#2D6A4F'
ACCENT     = '#1B7A4A'
ACCENT2    = '#2196A6'
ACCENT3    = '#E67E22'
ACCENT4    = '#8E44AD'

# Palette matplotlib
C_GREEN  = '#52B788'
C_TEAL   = '#48CAE4'
C_ORANGE = '#F4A261'
C_PURPLE = '#B19CD8'
C_LIME   = '#90EE90'
C_PINK   = '#FFB3BA'

MPL_BG   = '#FAFFFE'   # fond graphique (quasi blanc)
MPL_AXES = '#F0F9F4'   # fond axes


class VisualizeFrame(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent, bg=BG_MAIN)
        self.seq     = ''
        self.results = {}
        self.fig     = None
        self.canvas  = None

        if not MATPLOTLIB_OK:
            tk.Label(self, text='Matplotlib non installé.\npip install matplotlib',
                     bg=BG_MAIN, fg=FG_DIM, font=('Georgia', 12)).pack(expand=True)
            return
        self._build()

    def _build(self):
        # Barre de sélection
        ctrl = tk.Frame(self, bg=BG_PANEL,
                        highlightthickness=1, highlightbackground=BORDER)
        ctrl.pack(fill='x', padx=6, pady=(6, 0))

        tk.Label(ctrl, text='Affichage :',
                 bg=BG_PANEL, fg=FG_LABEL, font=('Georgia', 9, 'bold')
                 ).pack(side='left', padx=12, pady=8)

        self.graph_var = tk.StringVar(value='Carte de séquence')
        choices = [
            'Carte de séquence',
            'GC% (fenêtre glissante)',
            'Composition nucléotidique',
            'Longueurs des ORFs',
        ]
        cb = ttk.Combobox(ctrl, textvariable=self.graph_var,
                          values=choices, state='readonly', width=30,
                          font=('Georgia', 9))
        cb.pack(side='left', padx=4, pady=6)
        cb.bind('<<ComboboxSelected>>', lambda e: self._refresh())

        # Zone matplotlib
        self.fig_frame = tk.Frame(self, bg=BG_SURFACE,
                                   highlightthickness=1, highlightbackground=BORDER)
        self.fig_frame.pack(fill='both', expand=True, padx=6, pady=6)

    def draw(self, seq, results):
        self.seq     = seq
        self.results = results
        self._refresh()

    def _refresh(self):
        if not MATPLOTLIB_OK:
            return
        choice = self.graph_var.get()
        for w in self.fig_frame.winfo_children():
            w.destroy()
        if self.fig:
            plt.close(self.fig)

        if choice == 'Carte de séquence':
            self._draw_sequence_map()
        elif choice == 'GC% (fenêtre glissante)':
            self._draw_gc_sliding()
        elif choice == 'Composition nucléotidique':
            self._draw_composition()
        elif choice == 'Longueurs des ORFs':
            self._draw_orf_lengths()

    # ── Carte de séquence ─────────────────────────────────────
    def _draw_sequence_map(self):
        seq         = self.seq
        orfs        = self.results.get('orfs', [])
        promoters   = self.results.get('promoters', [])
        sd_sites    = self.results.get('sd_sites', [])
        terminators = self.results.get('terminators', [])

        self.fig, ax = plt.subplots(figsize=(11, 4), facecolor=MPL_BG)
        ax.set_facecolor(MPL_AXES)
        ax.set_xlim(0, len(seq))
        ax.set_ylim(-0.5, 1.8)
        ax.set_xlabel('Position (pb)', color=FG_DIM, fontsize=9, fontfamily='serif')
        ax.tick_params(colors=FG_DIM, labelsize=8)
        for spine in ax.spines.values():
            spine.set_edgecolor(BORDER)
            spine.set_linewidth(0.8)

        # Ligne de base
        ax.axhline(y=0.5, color=BORDER, linewidth=2, zorder=1)
        ax.text(0, 0.5, "5'", va='center', ha='right', color=FG_DIM, fontsize=8, fontfamily='serif')
        ax.text(len(seq), 0.5, "3'", va='center', ha='left', color=FG_DIM, fontsize=8, fontfamily='serif')

        # ORFs
        for i, orf in enumerate(orfs[:8]):
            y     = 0.68 if orf['strand'] == '+' else 0.32
            color = C_GREEN if orf['strand'] == '+' else C_TEAL
            alpha = 0.9 if i == 0 else 0.6
            rect  = mpatches.FancyBboxPatch(
                (orf['start'], y-0.1), orf['length'], 0.2,
                boxstyle='round,pad=0.005',
                facecolor=color, edgecolor='#FFFFFF',
                linewidth=0.8, alpha=alpha, zorder=2
            )
            ax.add_patch(rect)
            if orf['length'] > len(seq)*0.04:
                ax.text(orf['start']+orf['length']/2, y, f"ORF{i+1}",
                        ha='center', va='center',
                        color='#FFFFFF', fontsize=7, fontweight='bold',
                        fontfamily='serif', zorder=3)

        # Promoteurs
        for j, p in enumerate(promoters[:6]):
            ax.axvline(x=p['pos35'], color=ACCENT2,  linewidth=1.2, alpha=0.7, zorder=2)
            ax.axvline(x=p['pos10'], color=C_ORANGE,  linewidth=1.2, alpha=0.7, zorder=2)

        # SD
        for sd in sd_sites[:10]:
            ax.plot(sd['pos'], 1.35, 'v', color=C_PURPLE, markersize=7,
                    alpha=0.85, zorder=3)

        # Terminateurs
        for t in terminators[:5]:
            ax.axvspan(t['pos'], t['end'], ymin=0.05, ymax=0.95,
                       facecolor=C_ORANGE, alpha=0.08, zorder=1)

        # Légende
        legend_items = [
            mpatches.Patch(facecolor=C_GREEN,  label="ORF sens (+)"),
            mpatches.Patch(facecolor=C_TEAL,   label="ORF antisens (-)"),
            mpatches.Patch(facecolor=ACCENT2,  label="Promoteur Box-35", alpha=0.7),
            mpatches.Patch(facecolor=C_ORANGE, label="Promoteur Box-10", alpha=0.7),
            mpatches.Patch(facecolor=C_PURPLE, label="Shine-Dalgarno"),
            mpatches.Patch(facecolor=C_ORANGE, label="Terminateur", alpha=0.3),
        ]
        legend = ax.legend(handles=legend_items, loc='upper right',
                           facecolor=BG_SURFACE, edgecolor=BORDER,
                           labelcolor=FG_MAIN, fontsize=7.5,
                           framealpha=0.95, prop={'family': 'serif'})

        ax.set_title(f'Carte de la séquence — {len(seq)} pb',
                     color=FG_MAIN, fontsize=10, fontfamily='serif', pad=8)
        ax.set_yticks([])
        self._embed(self.fig)

    # ── GC sliding window ─────────────────────────────────────
    def _draw_gc_sliding(self):
        from analysis.statistics import sliding_gc
        positions, gc_values = sliding_gc(self.seq, window=100, step=10)

        self.fig, ax = plt.subplots(figsize=(11, 4), facecolor=MPL_BG)
        ax.set_facecolor(MPL_AXES)

        ax.plot(positions, gc_values, color=C_GREEN, linewidth=1.8, zorder=3)
        ax.fill_between(positions, gc_values, alpha=0.18, color=C_GREEN, zorder=2)
        ax.axhline(y=50, color=BORDER, linewidth=1, linestyle='--', alpha=0.8, zorder=1)

        gc_avg = self.results.get('stats', {}).get('gc', 0)
        ax.axhline(y=gc_avg, color=C_ORANGE, linewidth=1.2, linestyle='--',
                   alpha=0.9, label=f'GC% moyen : {gc_avg:.1f}%', zorder=2)

        ax.set_xlabel('Position (pb)', color=FG_DIM, fontsize=9, fontfamily='serif')
        ax.set_ylabel('GC (%)', color=FG_DIM, fontsize=9, fontfamily='serif')
        ax.set_title('Teneur en GC — fenêtre glissante de 100 pb',
                     color=FG_MAIN, fontsize=10, fontfamily='serif')
        ax.tick_params(colors=FG_DIM, labelsize=8)
        for spine in ax.spines.values():
            spine.set_edgecolor(BORDER)
        ax.legend(facecolor=BG_SURFACE, edgecolor=BORDER,
                  labelcolor=FG_MAIN, fontsize=8, prop={'family':'serif'})
        self._embed(self.fig)

    # ── Composition ───────────────────────────────────────────
    def _draw_composition(self):
        stats  = self.results.get('stats', {})
        counts = stats.get('counts', {})
        pcts   = stats.get('pcts', {})
        bases  = ['A', 'T', 'G', 'C']
        colors = [C_ORANGE, C_PURPLE, C_GREEN, C_TEAL]

        self.fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), facecolor=MPL_BG)
        for ax in (ax1, ax2):
            ax.set_facecolor(MPL_AXES)
            ax.tick_params(colors=FG_DIM, labelsize=8)
            for spine in ax.spines.values():
                spine.set_edgecolor(BORDER)

        # Histogramme
        bars = ax1.bar(bases, [counts.get(b, 0) for b in bases],
                       color=colors, edgecolor='#FFFFFF', linewidth=1.2,
                       width=0.6)
        ax1.set_title('Comptage des bases', color=FG_MAIN, fontsize=10, fontfamily='serif')
        ax1.set_ylabel('Nombre', color=FG_DIM, fontsize=9, fontfamily='serif')
        for bar, base in zip(bars, bases):
            height = bar.get_height()
            ax1.text(bar.get_x()+bar.get_width()/2,
                     height + max(counts.values(), default=1)*0.015,
                     f'{pcts.get(base,0):.1f}%',
                     ha='center', va='bottom', color=FG_MAIN,
                     fontsize=8, fontfamily='serif')

        # Camembert
        wedges, _ = ax2.pie([pcts.get(b, 0) for b in bases],
                             colors=colors, startangle=90,
                             wedgeprops=dict(edgecolor='#FFFFFF', linewidth=1.5))
        ax2.legend(wedges, [f'{b}: {pcts.get(b,0):.1f}%' for b in bases],
                   loc='lower center', ncol=2,
                   facecolor=BG_SURFACE, edgecolor=BORDER,
                   labelcolor=FG_MAIN, fontsize=8, prop={'family':'serif'})
        ax2.set_title('Répartition (%)', color=FG_MAIN, fontsize=10, fontfamily='serif')

        self.fig.suptitle('Composition Nucléotidique', color=FG_MAIN,
                           fontsize=11, fontfamily='serif')
        self._embed(self.fig)

    # ── Longueurs ORFs ────────────────────────────────────────
    def _draw_orf_lengths(self):
        orfs = self.results.get('orfs', [])
        if not orfs:
            tk.Label(self.fig_frame, text='Aucun ORF à afficher.',
                     bg=BG_SURFACE, fg=FG_DIM, font=('Georgia', 11)).pack(expand=True)
            return

        self.fig, ax = plt.subplots(figsize=(11, 4), facecolor=MPL_BG)
        ax.set_facecolor(MPL_AXES)

        lengths = [o['length'] for o in orfs[:20]]
        frames  = [str(o['frame']) for o in orfs[:20]]
        colors  = [C_GREEN if o['strand']=='+' else C_TEAL for o in orfs[:20]]

        bars = ax.barh(range(len(lengths)), lengths, color=colors,
                       edgecolor='#FFFFFF', linewidth=0.8, height=0.6)
        ax.set_yticks(range(len(frames)))
        ax.set_yticklabels([f"Frame {f}" for f in frames],
                            color=FG_DIM, fontsize=8, fontfamily='serif')
        ax.set_xlabel('Longueur (pb)', color=FG_DIM, fontsize=9, fontfamily='serif')
        ax.set_title(f'Longueurs des {len(lengths)} premiers ORFs',
                     color=FG_MAIN, fontsize=10, fontfamily='serif')
        ax.tick_params(colors=FG_DIM, labelsize=8)
        for spine in ax.spines.values():
            spine.set_edgecolor(BORDER)

        for bar, length in zip(bars, lengths):
            ax.text(bar.get_width()+5, bar.get_y()+bar.get_height()/2,
                    f'{length} pb', va='center', color=FG_MAIN,
                    fontsize=7, fontfamily='serif')

        legend_items = [
            mpatches.Patch(facecolor=C_GREEN, label="Brin sens (+)"),
            mpatches.Patch(facecolor=C_TEAL,  label="Brin antisens (-)"),
        ]
        ax.legend(handles=legend_items, facecolor=BG_SURFACE, edgecolor=BORDER,
                  labelcolor=FG_MAIN, fontsize=8, prop={'family':'serif'})
        self._embed(self.fig)

    # ── Embed ─────────────────────────────────────────────────
    def _embed(self, fig):
        plt.tight_layout()

        # Appliquer le style à la figure
        fig.patch.set_linewidth(0)

        canvas = FigureCanvasTkAgg(fig, master=self.fig_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)

        toolbar_frame = tk.Frame(self.fig_frame, bg=BG_PANEL,
                                  highlightthickness=1,
                                  highlightbackground=BORDER)
        toolbar_frame.pack(fill='x')

        toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
        toolbar.config(bg=BG_PANEL)
        for child in toolbar.winfo_children():
            try:
                child.config(bg=BG_PANEL)
            except Exception:
                pass
        toolbar.update()
        self.canvas = canvas