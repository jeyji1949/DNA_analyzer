# ============================================================
# Panneau gauche — saisie & options — Thème biologique clair
# ============================================================

import tkinter as tk
from tkinter import filedialog
import math

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

DEMO_SEQ = (
    'CATTTCTCTTAAGATTTATTCTATCTTAACACAACAACTTTTAATAAAAGATATGTAGAT'
    'TACAATTTAAATAGATTGTAATATTTGTAACACTAACATTAATATAGTTGTTATTTTTGT'
    'TACATAAACCACTAATAACTCATAATCTTTTAAAACTTATATTTGAGATAACATCAACTT'
    'TACATTACAAGTTATAAAACAAAAGAAGTGGGACACAGAATTCGTCTTGAACACTGTGTC'
    'CCACCTCGTCCCCAAAACTTGCTCTGTCCGTAGAAAAATAAAAAGGGGCCCCCTTTGTTG'
)


class InputFrame(tk.Frame):
    def __init__(self, parent, on_analyze_callback):
        super().__init__(parent, bg=BG_PANEL)
        self.callback = on_analyze_callback
        self._build()

    def _build(self):
        # Scrollable content via Canvas + Frame
        outer = tk.Frame(self, bg=BG_PANEL)
        outer.pack(fill='both', expand=True)

        # ── Canvas décoratif en fond (ne capte pas les clics) ──
        self.bio_canvas = tk.Canvas(outer, bg=BG_PANEL, bd=0,
                                     highlightthickness=0, width=40)
        self.bio_canvas.pack(side='right', fill='y')
        self.bio_canvas.bind('<Configure>', lambda e: self._draw_bio_side())

        # ── Contenu principal ──────────────────────────────────
        content = tk.Frame(outer, bg=BG_PANEL)
        content.pack(side='left', fill='both', expand=True)

        # Titre
        tk.Label(content, text="Séquence ADN (5'→3')",
                 bg=BG_PANEL, fg=ACCENT,
                 font=('Georgia', 11, 'bold')).pack(anchor='w', padx=12, pady=(12, 4))

        # Zone de texte
        txt_outer = tk.Frame(content, bg=BORDER, padx=1, pady=1)
        txt_outer.pack(fill='x', padx=10, pady=(0, 4))
        txt_inner = tk.Frame(txt_outer, bg=BG_SURFACE)
        txt_inner.pack(fill='x')

        self.seq_text = tk.Text(
            txt_inner, height=8,
            bg=BG_SURFACE, fg='#1B5E3B',
            font=('Courier New', 9),
            insertbackground=ACCENT,
            selectbackground='#C8E6C9',
            selectforeground=FG_MAIN,
            wrap='word', bd=0, padx=8, pady=8,
            relief='flat'
        )
        self.seq_text.pack(fill='x')
        self.seq_text.insert('end', DEMO_SEQ)
        self.seq_text.bind('<KeyRelease>', self._update_stats)

        # Stats
        self.stats_label = tk.Label(content, text='',
                                    bg=BG_PANEL, fg=FG_DIM,
                                    font=('Georgia', 8, 'italic'))
        self.stats_label.pack(anchor='w', padx=12, pady=(2, 6))
        self._update_stats()

        # Boutons utilitaires
        btn_row = tk.Frame(content, bg=BG_PANEL)
        btn_row.pack(fill='x', padx=10, pady=(0, 8))

        self._mkbtn(btn_row, '📁 Importer FASTA', self._import_fasta,
                    ACCENT2, '#FFFFFF').pack(side='left', fill='x', expand=True, padx=(0, 3))
        self._mkbtn(btn_row, '↔ RC', self._reverse_complement,
                    '#4A7C59', '#FFFFFF').pack(side='left', padx=(0, 3))
        self._mkbtn(btn_row, '✕', self._clear,
                    '#C0392B', '#FFFFFF').pack(side='left')

        # Séparateur
        self._separator(content)

        # Options
        tk.Label(content, text='Analyses à effectuer :',
                 bg=BG_PANEL, fg=FG_LABEL,
                 font=('Georgia', 9, 'bold')).pack(anchor='w', padx=12, pady=(4, 4))

        self.options = {}
        opts = [
            ('ORFs',           '🔍  Détection ORFs (6 cadres)'),
            ('Promoteurs',     '📍  Promoteurs -10 / -35'),
            ('Shine-Dalgarno', '🎯  Shine-Dalgarno'),
            ('Terminateurs',   '🔚  Terminateurs Rho-indép.'),
            ('Restriction',    '✂️   Sites de restriction'),
        ]
        for key, label in opts:
            var = tk.BooleanVar(value=True)
            tk.Checkbutton(
                content, text=label, variable=var,
                bg=BG_PANEL, fg=FG_MAIN,
                activebackground=BG_PANEL, activeforeground=ACCENT,
                selectcolor=BG_SURFACE,
                font=('Georgia', 9), cursor='hand2'
            ).pack(anchor='w', padx=16, pady=1)
            self.options[key] = var

        # Longueur min ORF
        self._separator(content)
        row = tk.Frame(content, bg=BG_PANEL)
        row.pack(fill='x', padx=12, pady=(4, 8))
        tk.Label(row, text='Longueur min. ORF (pb) :',
                 bg=BG_PANEL, fg=FG_DIM,
                 font=('Georgia', 9)).pack(side='left')
        entry_frame = tk.Frame(row, bg=BORDER, padx=1, pady=1)
        entry_frame.pack(side='right')
        self.min_orf = tk.Entry(entry_frame, width=6,
                                bg=BG_SURFACE, fg=FG_MAIN,
                                insertbackground=ACCENT,
                                font=('Courier New', 9), bd=0)
        self.min_orf.pack()
        self.min_orf.insert(0, '30')

        # ── Bouton LANCER ──────────────────────────────────────
        run_btn = tk.Button(
            content,
            text='▶  LANCER L\'ANALYSE',
            bg=ACCENT, fg='#FFFFFF',
            font=('Georgia', 11, 'bold'),
            relief='flat', cursor='hand2',
            activebackground='#145C38',
            activeforeground='#FFFFFF',
            command=self._run,          # <-- connexion correcte
            bd=0
        )
        run_btn.pack(fill='x', padx=10, ipady=10, pady=(4, 12))
        run_btn.bind('<Enter>', lambda e: run_btn.config(bg='#145C38'))
        run_btn.bind('<Leave>', lambda e: run_btn.config(bg=ACCENT))

    # ── Bande décorative latérale avec motifs biologiques ─────
    def _draw_bio_side(self):
        c = self.bio_canvas
        c.delete('bio')
        w = c.winfo_width()
        h = c.winfo_height()
        if w < 4 or h < 10:
            return
        cx = w / 2
        # Hélice verticale
        step = 20
        amp  = 6
        freq = 0.14
        pts1, pts2 = [], []
        for y in range(0, h + step, step):
            pts1.append((cx + amp * math.sin(freq * y),       y))
            pts2.append((cx + amp * math.sin(freq * y + math.pi), y))
        for i in range(len(pts1) - 1):
            c.create_line(pts1[i][0], pts1[i][1], pts1[i+1][0], pts1[i+1][1],
                           fill='#B7E4C7', width=1.5, smooth=True, tags='bio')
        for i in range(len(pts2) - 1):
            c.create_line(pts2[i][0], pts2[i][1], pts2[i+1][0], pts2[i+1][1],
                           fill='#D8F3DC', width=1.5, smooth=True, tags='bio')
        # Barreaux
        for i, (p1, p2) in enumerate(zip(pts1[::2], pts2[::2])):
            c.create_line(p1[0], p1[1], p2[0], p2[1],
                           fill='#95D5B2', width=1, tags='bio')
            r = 2
            c.create_oval(p1[0]-r, p1[1]-r, p1[0]+r, p1[1]+r,
                           fill='#C8E6C9', outline='', tags='bio')

    # ── Séparateur décoratif ──────────────────────────────────
    def _separator(self, parent):
        sep_frame = tk.Frame(parent, bg=BG_PANEL, height=18)
        sep_frame.pack(fill='x', padx=10, pady=2)
        sep_canvas = tk.Canvas(sep_frame, bg=BG_PANEL, height=18,
                                bd=0, highlightthickness=0)
        sep_canvas.pack(fill='x')

        def draw(event=None, cv=sep_canvas):
            cv.delete('all')
            w2 = cv.winfo_width()
            if w2 < 10:
                return
            mid = 9
            cv.create_line(8, mid, w2 - 8, mid, fill=BORDER, width=1)
            cv.create_oval(w2//2-4, mid-4, w2//2+4, mid+4,
                            fill='#74C69D', outline=BORDER, width=1)

        sep_canvas.bind('<Configure>', draw)

    # ── Helpers ───────────────────────────────────────────────
    def _mkbtn(self, parent, text, cmd, bg, fg):
        return tk.Button(parent, text=text, command=cmd,
                         bg=bg, fg=fg, font=('Georgia', 8),
                         relief='flat', cursor='hand2',
                         activebackground='#145C38',
                         activeforeground='#FFFFFF',
                         padx=6, pady=5, bd=0)

    def _update_stats(self, event=None):
        seq = self.get_sequence()
        if not seq:
            self.stats_label.config(text='')
            return
        a, t, g, c_ = seq.count('A'), seq.count('T'), seq.count('G'), seq.count('C')
        gc = (g + c_) / len(seq) * 100
        self.stats_label.config(
            text=f'{len(seq)} pb  •  GC: {gc:.1f}%  •  A:{a}  T:{t}  G:{g}  C:{c_}'
        )

    def _import_fasta(self):
        path = filedialog.askopenfilename(
            filetypes=[('FASTA', '*.fasta *.fa *.txt *.seq'), ('Tous', '*.*')]
        )
        if path:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            lines = [l.strip() for l in content.splitlines() if not l.startswith('>')]
            self.set_sequence(''.join(lines).upper())

    def _reverse_complement(self):
        from analysis.orf_finder import reverse_complement
        seq = self.get_sequence()
        if seq:
            self.set_sequence(reverse_complement(seq))

    def _clear(self):
        self.seq_text.delete('1.0', 'end')
        self._update_stats()

    def _run(self):
        seq  = self.get_sequence()
        opts = self.get_options()
        self.callback(seq, opts)

    # ── Interface publique ─────────────────────────────────────
    def get_sequence(self):
        return (self.seq_text.get('1.0', 'end').strip()
                .upper().replace(' ', '').replace('\n', ''))

    def set_sequence(self, seq):
        self.seq_text.delete('1.0', 'end')
        self.seq_text.insert('end', seq)
        self._update_stats()

    def get_options(self):
        opts = {k: v.get() for k, v in self.options.items()}
        try:
            opts['min_orf'] = int(self.min_orf.get())
        except ValueError:
            opts['min_orf'] = 90
        return opts