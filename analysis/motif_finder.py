# ============================================================
# Recherche de motifs régulateurs
# Promoteurs, Shine-Dalgarno, Terminateurs Rho-indépendants
# ============================================================

from analysis.orf_finder import reverse_complement


def hamming(s1, s2):
    """Nombre de différences (mismatches) entre deux séquences de même longueur."""
    if len(s1) != len(s2):
        return max(len(s1), len(s2))
    return sum(a != b for a, b in zip(s1, s2))


# ------------------------------------------------------------------
# PROMOTEURS BACTÉRIENS (-10 / -35)
# ------------------------------------------------------------------

def find_promoters(seq, max_mm_box10=2, max_mm_box35=2):
    """
    Détecte les promoteurs bactériens (paires box -10 / box -35).

    Consensus :
        Box -10 : TATAAT
        Box -35 : TTGACA
    Espacement optimal : 17 pb (accepté : 14 à 22 pb)

    Retourne :
        dict avec clés 'promoters', 'box10', 'box35'
    """
    box10_consensus = 'TATAAT'
    box35_consensus = 'TTGACA'

    box10_hits = []
    box35_hits = []

    for i in range(len(seq) - 5):
        sub = seq[i:i+6]
        d10 = hamming(sub, box10_consensus)
        d35 = hamming(sub, box35_consensus)
        if d10 <= max_mm_box10:
            box10_hits.append({'pos': i, 'seq': sub, 'dist': d10})
        if d35 <= max_mm_box35:
            box35_hits.append({'pos': i, 'seq': sub, 'dist': d35})

    promoters = []
    for b35 in box35_hits:
        for b10 in box10_hits:
            spacing = b10['pos'] - (b35['pos'] + 6)
            if 14 <= spacing <= 22:
                total_mm = b35['dist'] + b10['dist']
                quality = (
                    'Consensus parfait' if total_mm == 0 else
                    'Fort'              if total_mm <= 2 else
                    'Modéré'            if total_mm <= 4 else
                    'Faible'
                )
                promoters.append({
                    'pos35':   b35['pos'],
                    'seq35':   b35['seq'],
                    'dist35':  b35['dist'],
                    'pos10':   b10['pos'],
                    'seq10':   b10['seq'],
                    'dist10':  b10['dist'],
                    'spacing': spacing,
                    'quality': quality,
                })

    # Tri par qualité
    promoters.sort(key=lambda x: x['dist35'] + x['dist10'])

    return {
        'promoters': promoters,
        'box10':     box10_hits,
        'box35':     box35_hits,
    }


# ------------------------------------------------------------------
# SHINE-DALGARNO
# ------------------------------------------------------------------

def find_shine_dalgarno(seq, max_mm=2):
    """
    Détecte les sites Shine-Dalgarno (SD).

    Consensus : AGGAGG
    Position   : 7 à 9 pb en amont du codon ATG de départ.

    Retourne :
        Liste de dicts avec position SD, séquence, position ATG, espacement
    """
    sd_consensus = 'AGGAGG'
    results = []
    seen_pos = set()

    for i in range(len(seq) - 2):
        if seq[i:i+3] == 'ATG':
            for dist in range(7, 10):           # espacement 7, 8 ou 9
                sd_start = i - dist - 6
                if sd_start < 0:
                    continue
                candidate = seq[sd_start:sd_start+6]
                mm = hamming(candidate, sd_consensus)
                if mm <= max_mm and sd_start not in seen_pos:
                    quality = (
                        'Consensus' if mm == 0 else
                        'Fort'      if mm == 1 else
                        'Modéré'
                    )
                    results.append({
                        'pos':     sd_start,
                        'seq':     candidate,
                        'atg_pos': i,
                        'spacing': dist,
                        'dist':    mm,
                        'quality': quality,
                    })
                    seen_pos.add(sd_start)

    return sorted(results, key=lambda x: x['dist'])


# ------------------------------------------------------------------
# TERMINATEURS RHO-INDÉPENDANTS
# ------------------------------------------------------------------

def find_terminators(seq, min_stem=4, max_stem=12, min_poly_t=3):
    """
    Détecte les terminateurs de transcription Rho-indépendants.

    Structure : tige-boucle palindromique (stem-loop) suivie d'une
    séquence polyT (≥ min_poly_t thymine consécutives).

    Retourne :
        Liste de dicts décrivant chaque terminateur trouvé
    """
    terminators = []
    i = 0
    while i < len(seq) - 20:
        found = False
        for stem_len in range(min_stem, max_stem + 1):
            if found:
                break
            arm1 = seq[i:i+stem_len]
            rc_arm1 = reverse_complement(arm1)

            for loop_len in range(3, 9):
                arm2_start = i + stem_len + loop_len
                if arm2_start + stem_len > len(seq):
                    break
                arm2 = seq[arm2_start:arm2_start+stem_len]

                if hamming(arm2, rc_arm1) <= 1:
                    # Vérifie le polyT après la tige-boucle
                    after = arm2_start + stem_len
                    poly_t_seq = seq[after:after+8]
                    t_count = sum(1 for c in poly_t_seq[:6] if c == 'T')

                    if t_count >= min_poly_t:
                        terminators.append({
                            'pos':      i,
                            'end':      after + 6,
                            'arm1':     arm1,
                            'loop':     seq[i+stem_len:arm2_start],
                            'arm2':     arm2,
                            'stem_len': stem_len,
                            'loop_len': loop_len,
                            'poly_t':   poly_t_seq[:6],
                            'seq':      seq[i:after+6],
                        })
                        i = after        # évite les chevauchements
                        found = True
                        break
        if not found:
            i += 1

    return terminators


# ------------------------------------------------------------------
# SITES DE RESTRICTION
# ------------------------------------------------------------------

def find_restriction_sites(seq, sites_dict):
    """
    Localise tous les sites de restriction dans la séquence.

    Paramètres :
        seq        : séquence ADN
        sites_dict : dict {nom_enzyme: (séquence_site, pos_coupure)}

    Retourne :
        dict {nom_enzyme: {'site': str, 'positions': [int], 'count': int}}
    """
    results = {}
    for enzyme, (site, cut_pos) in sites_dict.items():
        positions = []
        idx = seq.find(site)
        while idx != -1:
            positions.append(idx + 1)       # position 1-based
            idx = seq.find(site, idx + 1)
        if positions:
            results[enzyme] = {
                'site':      site,
                'positions': positions,
                'count':     len(positions),
                'cut_pos':   cut_pos,
            }
    return results
