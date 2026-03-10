# ============================================================
# Statistiques et propriétés physico-chimiques de la séquence
# ============================================================

from collections import Counter


def nucleotide_composition(seq):
    """Retourne le comptage et le pourcentage de chaque base."""
    seq = seq.upper()
    n = len(seq)
    counts = {b: seq.count(b) for b in 'ATGC'}
    pcts   = {b: round(counts[b] / n * 100, 2) if n else 0 for b in 'ATGC'}
    gc     = round((counts['G'] + counts['C']) / n * 100, 2) if n else 0
    at     = round((counts['A'] + counts['T']) / n * 100, 2) if n else 0
    return {
        'counts':  counts,
        'pcts':    pcts,
        'gc':      gc,
        'at':      at,
        'length':  n,
        'ratio_gc': round(counts['G'] / counts['C'], 3) if counts['C'] else 0,
        'ratio_at': round(counts['A'] / counts['T'], 3) if counts['T'] else 0,
    }


def calc_tm_wallace(seq):
    """
    Température de fusion selon la règle de Wallace.
    Fiable pour séquences courtes (< 14 pb).
    Tm = 2(A+T) + 4(G+C)
    """
    seq = seq.upper()
    a = seq.count('A')
    t = seq.count('T')
    g = seq.count('G')
    c = seq.count('C')
    return 2 * (a + t) + 4 * (g + c)


def calc_tm_nearest_neighbor(seq):
    """
    Température de fusion par la méthode nearest-neighbor.
    Plus précise que Wallace, surtout pour séquences > 14 pb.
    """
    seq = seq.upper()

    # Paramètres thermodynamiques (dH en cal/mol, dS en cal/mol/K)
    params = {
        'AA': (-7900,  -22.2), 'AT': (-7200,  -20.4),
        'AC': (-8400,  -22.4), 'AG': (-7800,  -21.0),
        'TA': (-7200,  -21.3), 'TT': (-7900,  -22.2),
        'TC': (-8200,  -22.2), 'TG': (-8500,  -22.7),
        'CA': (-8500,  -22.7), 'CT': (-7800,  -21.0),
        'CC': (-8000,  -19.9), 'CG': (-10600, -27.2),
        'GA': (-8200,  -22.2), 'GT': (-8400,  -22.4),
        'GC': (-9800,  -24.4), 'GG': (-8000,  -19.9),
    }

    R = 1.987       # constante des gaz (cal/mol/K)
    C = 250e-9      # concentration oligonucléotide (250 nM)

    dH = sum(params.get(seq[i:i+2], (0, 0))[0] for i in range(len(seq) - 1))
    dS = sum(params.get(seq[i:i+2], (0, 0))[1] for i in range(len(seq) - 1))

    # Corrections d'initiation
    dH += -3200
    dS += -9.0

    if dS == 0:
        return 0.0

    Tm = dH / (dS + R * (2.303 * (1.0 / C))) - 273.15
    return round(Tm, 1)


def molecular_weight(seq, double_stranded=False):
    """
    Masse moléculaire approximative de la séquence ADN.
    Poids moyens des nucléotides monophosphates (Da) :
        dAMP: 313.2, dTMP: 304.2, dGMP: 329.2, dCMP: 289.2
    """
    seq = seq.upper()
    weights = {'A': 313.2, 'T': 304.2, 'G': 329.2, 'C': 289.2}
    mw = sum(weights.get(b, 0) for b in seq) - (len(seq) - 1) * 18.0
    if double_stranded:
        mw *= 2
    return round(mw / 1000, 2)   # retourne en kDa


def sliding_gc(seq, window=100, step=10):
    """
    Calcule le %GC dans une fenêtre glissante.

    Retourne :
        positions : liste des positions centrales
        gc_values : liste des %GC correspondants
    """
    seq = seq.upper()
    positions, gc_values = [], []
    for i in range(0, len(seq) - window, step):
        w = seq[i:i+window]
        gc = (w.count('G') + w.count('C')) / window * 100
        positions.append(i + window // 2)
        gc_values.append(round(gc, 2))
    return positions, gc_values


def codon_usage(orf_seq):
    """
    Calcule la fréquence d'utilisation de chaque codon dans un ORF.

    Retourne :
        dict {codon: pourcentage}  trié par fréquence décroissante
    """
    orf_seq = orf_seq.upper()
    codons = [orf_seq[i:i+3] for i in range(0, len(orf_seq) - 2, 3)]
    counts = Counter(codons)
    total = sum(counts.values())
    return {
        codon: round(count / total * 100, 2)
        for codon, count in counts.most_common()
    }


def protein_stats(protein_seq):
    """
    Statistiques basiques sur une séquence protéique.
    Retourne composition en AA et poids moléculaire estimé.
    """
    from data.codon_table import AA_HYDROPHOBIC, AA_POLAR, AA_CHARGED
    protein_seq = protein_seq.replace('*', '')
    n = len(protein_seq)
    if n == 0:
        return {}
    counts = Counter(protein_seq)
    hydro = sum(counts[aa] for aa in AA_HYDROPHOBIC)
    polar = sum(counts[aa] for aa in AA_POLAR)
    charg = sum(counts[aa] for aa in AA_CHARGED)
    return {
        'length':       n,
        'mw_kda':       round(n * 110 / 1000, 1),   # estimation grossière
        'hydrophobic':  round(hydro / n * 100, 1),
        'polar':        round(polar / n * 100, 1),
        'charged':      round(charg / n * 100, 1),
        'composition':  dict(counts.most_common()),
    }
