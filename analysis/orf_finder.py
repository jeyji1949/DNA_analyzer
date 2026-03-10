# ============================================================
# Détection des ORFs — 6 cadres de lecture
# Mode 1 : segments entre codons stop (depuis le 1er nucléotide)
# Mode 2 : ORFs classiques ATG→Stop
# ============================================================

from data.codon_table import CODON_TABLE, STOP_CODONS, START_CODON


def reverse_complement(seq):
    comp = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(comp.get(c,'N') for c in reversed(seq.upper()))


def translate(seq):
    return ''.join(CODON_TABLE.get(seq[i:i+3], '?') for i in range(0, len(seq)-2, 3))


# ──────────────────────────────────────────────────────────────
# MODE 1 — SEGMENTS ENTRE CODONS STOP (depuis le 1er nucléotide)
# Répond à la question : "quels sont les ORFs possibles à partir
# du 1er nucléotide sur les 2 brins ?"
# ──────────────────────────────────────────────────────────────
def find_reading_frames(seq, min_length=30):
    """
    Détecte tous les segments de lecture sur les 6 cadres en
    partant du 1er nucléotide, délimités par des codons stop.
    Chaque segment peut ou non contenir un ATG interne.

    Retourne une liste triée par longueur décroissante.
    """
    seq = seq.upper().replace(' ','').replace('\n','')
    seq_len = len(seq)
    rc_seq  = reverse_complement(seq)
    results = []

    for strand, s in [('+', seq), ('-', rc_seq)]:
        for frame in range(3):
            codons      = [s[i:i+3] for i in range(frame, len(s)-2, 3)
                           if len(s[i:i+3]) == 3]
            seg_start   = 0   # index codon du début du segment courant

            for k, codon in enumerate(codons):
                is_last = (k == len(codons) - 1)
                if codon in STOP_CODONS or is_last:
                    end_k = k if codon in STOP_CODONS else k + 1
                    seg_codons = codons[seg_start:end_k]
                    seg_len_nt = len(seg_codons) * 3

                    if seg_len_nt >= min_length:
                        nt_start = frame + seg_start * 3
                        nt_end   = frame + end_k * 3
                        seg_seq  = s[nt_start:nt_end]

                        # Cherche le premier ATG dans le segment
                        first_atg = None
                        for ci, sc in enumerate(seg_codons):
                            if sc == 'ATG':
                                first_atg = seg_start + ci
                                break

                        prot = translate(seg_seq)
                        has_atg = first_atg is not None
                        atg_pos_nt = (frame + first_atg * 3) if has_atg else None

                        if strand == '+':
                            rs = nt_start
                            re = nt_end
                            atg_real = atg_pos_nt
                        else:
                            rs = seq_len - nt_end
                            re = seq_len - nt_start
                            atg_real = (seq_len - (frame + (first_atg+1)*3)) if has_atg else None

                        stop_codon = codon if codon in STOP_CODONS else '—'

                        results.append({
                            'start':      rs,
                            'end':        re,
                            'length':     seg_len_nt,
                            'frame':      frame + 1 if strand == '+' else -(frame + 1),
                            'strand':     strand,
                            'seq':        seg_seq,
                            'protein':    prot,
                            'num_aa':     len(prot),
                            'has_atg':    has_atg,
                            'atg_pos':    atg_real,
                            'stop_codon': stop_codon,
                            'mode':       'frame',
                        })
                    seg_start = k + 1

    return sorted(results, key=lambda x: -x['length'])


# ──────────────────────────────────────────────────────────────
# MODE 2 — ORFs classiques ATG → Stop
# ──────────────────────────────────────────────────────────────
def find_orfs(seq, min_length=30, start_codon='ATG'):
    """
    Détecte tous les ORFs ATG→Stop sur les 6 cadres.
    """
    seq = seq.upper().replace(' ','').replace('\n','')
    seq_len = len(seq)
    rc_seq  = reverse_complement(seq)
    orfs    = []

    for strand, s in [('+', seq), ('-', rc_seq)]:
        for frame in range(3):
            i = frame
            while i < len(s) - 2:
                if s[i:i+3] == start_codon:
                    start_pos = i
                    for j in range(i + 3, len(s) - 2, 3):
                        stop = s[j:j+3]
                        if stop in STOP_CODONS:
                            orf_seq = s[start_pos:j+3]
                            if len(orf_seq) >= min_length:
                                if strand == '+':
                                    rs, re = start_pos, j + 3
                                else:
                                    rs = seq_len - (j + 3)
                                    re = seq_len - start_pos
                                protein = translate(orf_seq)
                                orfs.append({
                                    'start':      rs,
                                    'end':        re,
                                    'length':     len(orf_seq),
                                    'frame':      frame + 1 if strand == '+' else -(frame + 1),
                                    'strand':     strand,
                                    'seq':        orf_seq,
                                    'protein':    protein.rstrip('*'),
                                    'stop_codon': stop,
                                    'num_aa':     len(orf_seq) // 3 - 1,
                                    'has_atg':    True,
                                    'mode':       'orf',
                                })
                            i = j + 3
                            break
                    else:
                        i += 3
                else:
                    i += 3

    return sorted(orfs, key=lambda x: -x['length'])


# ──────────────────────────────────────────────────────────────
# SCORE — ORF potentiellement codant
# ──────────────────────────────────────────────────────────────
def score_coding_potential(orf, seq_length):
    """
    Calcule un score de potentiel codant (0–100).
    Critères biologiques :
      - Longueur (plus c'est long, plus probable)
      - Présence d'un ATG de départ
      - Brin sens préféré
      - Longueur protéine ≥ 50 aa (seuil standard bactérien)
      - % de codons rares faible
    """
    score  = 0
    notes  = []

    # 1. Longueur relative à la séquence totale
    pct_seq = orf['length'] / seq_length * 100
    len_score = min(40, int(pct_seq * 2))
    score += len_score
    notes.append(f"Longueur {orf['length']} pb ({pct_seq:.1f}% de la séq.) → +{len_score}")

    # 2. Présence ATG
    if orf.get('has_atg'):
        score += 20
        notes.append('+20 (possède un ATG)')

    # 3. Longueur protéine ≥ 50 aa (seuil protéine bactérienne)
    if orf['num_aa'] >= 50:
        score += 20
        notes.append('+20 (≥50 acides aminés)')
    elif orf['num_aa'] >= 30:
        score += 10
        notes.append('+10 (≥30 acides aminés)')

    # 4. Brin sens légèrement favorisé
    if orf['strand'] == '+':
        score += 5
        notes.append('+5 (brin sens)')

    # 5. Codon stop présent (ORF complet)
    if orf.get('stop_codon') and orf['stop_codon'] != '—':
        score += 10
        notes.append('+10 (codon stop présent)')

    # 6. Pas de codon stop interne (protéine sans interruption)
    protein = orf.get('protein','')
    if '*' not in protein:
        score += 5
        notes.append('+5 (pas de stop interne)')

    return min(100, score), notes


def find_best_coding_orf(frames, seq_length):
    """
    Parmi tous les segments/ORFs, identifie le meilleur candidat codant.
    Retourne (best_orf, all_scored) où all_scored = [(score, notes, orf), ...]
    """
    scored = []
    for orf in frames:
        s, notes = score_coding_potential(orf, seq_length)
        scored.append((s, notes, orf))
    scored.sort(key=lambda x: -x[0])
    best = scored[0][2] if scored else None
    return best, scored