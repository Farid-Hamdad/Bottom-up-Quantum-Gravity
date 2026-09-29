from pathlib import Path
import argparse, csv, hashlib, json, sys
import numpy as np
from scipy.io import loadmat
from scipy.stats import pearsonr, spearmanr, rankdata

HERE = Path(__file__).resolve().parent
SCRIPTS = HERE.parent / 'scripts'
sys.path.insert(0, str(SCRIPTS))
from shadow_mi_test import mi_matrix
from raw_pairwise_mi_test import load_beta_profiles

DOI = '10.5281/zenodo.8279583'

def sha256(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def corr(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float('nan'), float('nan')
    return float(pearsonr(x, y).statistic), float(spearmanr(x, y).statistic)

def residualize(y, x):
    y = np.asarray(y, float)
    x = np.asarray(x, float)
    X = np.column_stack([np.ones(len(x)), x])
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    return y - X @ b

def partial_rank(x, y, z):
    rx = rankdata(x, method='average')
    ry = rankdata(y, method='average')
    rz = rankdata(z, method='average')
    ex, ey = residualize(rx, rz), residualize(ry, rz)
    if np.std(ex) == 0 or np.std(ey) == 0:
        return float('nan')
    return float(pearsonr(ex, ey).statistic)

def num(x):
    x = float(x)
    return x if np.isfinite(x) else None

def axes_from_history(h):
    h = np.asarray(h, int)
    a = np.full(h.shape, -1, int)
    a[np.isin(h, [1, 2])] = 0
    a[np.isin(h, [3, 4])] = 1
    a[np.isin(h, [5, 6])] = 2
    if np.any(a < 0):
        raise RuntimeError('unknown tomography code')
    return a

def audit_design(gs_raw, es_raw):
    g = loadmat(gs_raw)
    e = loadmat(es_raw)
    hg = np.asarray(g['combtomohistory'], int)
    he = np.asarray(e['combtomohistory'], int)
    xg = np.asarray(g['combdatasorted'])
    a = axes_from_history(hg)
    n = a.shape[1]
    single = np.array([[np.sum(a[:, j] == q) for q in range(3)] for j in range(n)])
    complete = np.eye(n, dtype=bool)
    exact27 = 0
    anyzero = 0
    mincell = 10**9
    maxcell = -1
    noncomplete = []

    for i in range(n):
        for j in range(i + 1, n):
            c = np.array([
                [np.sum((a[:, i] == q) & (a[:, j] == r)) for r in range(3)]
                for q in range(3)
            ])
            exact27 += int(np.all(c == 27))
            anyzero += int(np.any(c == 0))
            mincell = min(mincell, int(c.min()))
            maxcell = max(maxcell, int(c.max()))
            ok = bool(np.all(c > 0))
            complete[i, j] = complete[j, i] = ok
            if not ok:
                noncomplete.append((i + 1, j + 1))

    mismatches = []
    for i in range(n):
        for j in range(i + 1, n):
            if (not complete[i, j]) != (((i - j) % 5) == 0):
                mismatches.append((i + 1, j + 1))

    return {
        'gs_es_history_identical': bool(np.array_equal(hg, he)),
        'history_shape': list(hg.shape),
        'history_codes': np.unique(hg).astype(int).tolist(),
        'shot_columns': int(xg.shape[1] - 1),
        'first_column_is_0_to_242': bool(np.array_equal(xg[:, 0], np.arange(243))),
        'all_sites_81_81_81': bool(np.all(single == 81)),
        'n_pairs': int(n * (n - 1) // 2),
        'pairs_exact_27_each': int(exact27),
        'pairs_with_any_zero': int(anyzero),
        'min_joint_axis_count': int(mincell),
        'max_joint_axis_count': int(maxcell),
        'n_complete_pairs': int(np.sum(np.triu(complete, 1))),
        'n_noncomplete_pairs': int(len(noncomplete)),
        'mod5_rule_exact': bool(len(mismatches) == 0),
        'mod5_rule_mismatches': int(len(mismatches)),
    }

def channels(mi, sites):
    n = mi.shape[0]
    in_a = np.zeros(n, bool)
    in_a[sites] = True
    out = {}

    for name in ['all', 'in', 'cut']:
        sums, means, counts = [], [], []

        for j in sites:
            if name == 'all':
                mask = np.ones(n, bool)
                mask[j] = False
            elif name == 'in':
                mask = in_a.copy()
                mask[j] = False
            else:
                mask = ~in_a

            vals = mi[j, mask]
            good = np.isfinite(vals)
            s = float(np.nansum(vals))
            c = int(good.sum())

            sums.append(s)
            counts.append(c)
            means.append(s / c if c else float('nan'))

        out[name] = {
            'sum': np.asarray(sums),
            'mean': np.asarray(means),
            'count': np.asarray(counts),
        }

    return out

def edge_control(cut, invb):
    L = len(cut)
    p = np.arange(1, L + 1, dtype=float)
    g = 1.0 / (p * (L + 1.0 - p))
    g /= g.sum()

    rg, sg = corr(invb, g)
    rcg, scg = corr(cut, g)
    rci, sci = corr(cut, invb)

    e1, e2 = residualize(cut, g), residualize(invb, g)
    rp = float('nan') if np.std(e1) == 0 or np.std(e2) == 0 else float(pearsonr(e1, e2).statistic)

    return {
        'invbeta_vs_edge_pearson': num(rg),
        'invbeta_vs_edge_spearman': num(sg),
        'cut_vs_edge_pearson': num(rcg),
        'cut_vs_edge_spearman': num(scg),
        'cut_vs_invbeta_pearson': num(rci),
        'cut_vs_invbeta_spearman': num(sci),
        'partial_edge_pearson': num(rp),
        'partial_edge_spearman': num(partial_rank(cut, invb, g)),
    }

def distance_control(mi, sites, invb):
    n = mi.shape[0]
    in_a = np.zeros(n, bool)
    in_a[sites] = True
    outside = ~in_a

    kernel = {}
    for d in range(1, n):
        vals = [
            mi[i, i + d]
            for i in range(n - d)
            if not in_a[i]
            and not in_a[i + d]
            and np.isfinite(mi[i, i + d])
        ]
        if vals:
            kernel[d] = float(np.median(vals))

    obs, null = [], []

    for j in sites:
        o, z = [], []

        for k in np.where(outside)[0]:
            if not np.isfinite(mi[j, k]):
                continue

            d = abs(int(j) - int(k))

            if d in kernel:
                o.append(mi[j, k])
                z.append(kernel[d])

        obs.append(float(np.mean(o)))
        null.append(float(np.mean(z)))

    obs, null = np.asarray(obs), np.asarray(null)

    rci, sci = corr(obs, invb)
    rni, sni = corr(null, invb)
    rcn, scn = corr(obs, null)

    e1, e2 = residualize(obs, null), residualize(invb, null)
    rp = float('nan') if np.std(e1) == 0 or np.std(e2) == 0 else float(pearsonr(e1, e2).statistic)

    return {
        'cut_vs_invbeta_pearson': num(rci),
        'cut_vs_invbeta_spearman': num(sci),
        'null_vs_invbeta_pearson': num(rni),
        'null_vs_invbeta_spearman': num(sni),
        'cut_vs_null_pearson': num(rcn),
        'cut_vs_null_spearman': num(scn),
        'partial_distance_pearson': num(rp),
        'partial_distance_spearman': num(partial_rank(obs, invb, null)),
    }

def distance_permutation(mi, profiles, n_perm, seed):
    rng = np.random.default_rng(seed)
    n = mi.shape[0]
    by_d = {}

    for i in range(n):
        for j in range(i + 1, n):
            if np.isfinite(mi[i, j]):
                by_d.setdefault(j - i, []).append((i, j, float(mi[i, j])))

    rows = []

    for p in profiles:
        sites = (p['sites'] - 1).astype(int)
        invb = np.asarray(p['inv_beta_norm'], float)

        in_a = np.zeros(n, bool)
        in_a[sites] = True
        outside = ~in_a

        def cut_profile(mat):
            z = []

            for j in sites:
                vals = mat[j, outside]
                vals = vals[np.isfinite(vals)]
                z.append(float(np.mean(vals)))

            return np.asarray(z)

        obs = cut_profile(mi)
        ro, so = corr(obs, invb)

        rn = np.empty(n_perm)
        sn = np.empty(n_perm)

        for b in range(n_perm):
            m = np.full_like(mi, np.nan, dtype=float)
            np.fill_diagonal(m, 0.0)

            for pairs in by_d.values():
                vals = rng.permutation(np.asarray([x[2] for x in pairs], float))

                for (i, j, _), v in zip(pairs, vals):
                    m[i, j] = m[j, i] = v

            rn[b], sn[b] = corr(cut_profile(m), invb)

        rows.append({
            'cell_index': int(p['cell_index']),
            'n_sites': int(len(sites)),
            'obs_pearson': num(ro),
            'null_pearson_mean': num(rn.mean()),
            'null_pearson_sd': num(rn.std(ddof=1)),
            'p_pearson': num((1 + np.sum(rn >= ro)) / (n_perm + 1)),
            'obs_spearman': num(so),
            'null_spearman_mean': num(sn.mean()),
            'null_spearman_sd': num(sn.std(ddof=1)),
            'p_spearman': num((1 + np.sum(sn >= so)) / (n_perm + 1)),
        })

    return rows

def leave_a_out(mi, profiles, n_perm, seed):
    rng = np.random.default_rng(seed)
    n = mi.shape[0]
    rows = []

    for p in profiles:
        sites = (p['sites'] - 1).astype(int)
        invb = np.asarray(p['inv_beta_norm'], float)

        in_a = np.zeros(n, bool)
        in_a[sites] = True
        outside = ~in_a

        donor = {}

        for i in range(n):
            if in_a[i]:
                continue

            for j in range(i + 1, n):
                if not in_a[j] and np.isfinite(mi[i, j]):
                    donor.setdefault(j - i, []).append(float(mi[i, j]))

        specs, obs = {}, []

        for j in sites:
            vals, ds = [], []

            for k in np.where(outside)[0]:
                if not np.isfinite(mi[j, k]):
                    continue

                d = abs(int(j) - int(k))

                if d in donor and donor[d]:
                    vals.append(float(mi[j, k]))
                    ds.append(d)

            obs.append(float(np.mean(vals)))
            specs[int(j)] = ds

        obs = np.asarray(obs)
        ro, so = corr(obs, invb)

        rn = np.empty(n_perm)
        sn = np.empty(n_perm)

        for b in range(n_perm):
            prof = []

            for j in sites:
                draws = []

                for d in specs[int(j)]:
                    pool = donor[d]
                    draws.append(pool[rng.integers(0, len(pool))])

                prof.append(float(np.mean(draws)))

            rn[b], sn[b] = corr(np.asarray(prof), invb)

        rows.append({
            'cell_index': int(p['cell_index']),
            'n_sites': int(len(sites)),
            'obs_pearson': num(ro),
            'null_pearson_mean': num(rn.mean()),
            'null_pearson_sd': num(rn.std(ddof=1)),
            'p_pearson': num((1 + np.sum(rn >= ro)) / (n_perm + 1)),
            'obs_spearman': num(so),
            'null_spearman_mean': num(sn.mean()),
            'null_spearman_sd': num(sn.std(ddof=1)),
            'p_spearman': num((1 + np.sum(sn >= so)) / (n_perm + 1)),
        })

    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-root', type=Path, required=True)
    ap.add_argument('--n-perm', type=int, default=5000)
    args = ap.parse_args()

    root = args.data_root.resolve()
    outdir = HERE / 'results'
    outdir.mkdir(parents=True, exist_ok=True)

    gs_raw = root / 'RawData' / 'dataRawParam_GS_Delta_1.mat'
    es_raw = root / 'RawData' / 'dataRawParam_ES_Delta_1.mat'
    gs_beta = root / 'AnalyzedData' / 'Figure2' / 'GroundStateBetas.mat'
    es_beta = root / 'AnalyzedData' / 'Figure2' / 'ExcitedStateBetas.mat'
    readme = root / 'Readme.docx'

    inputs = [gs_raw, es_raw, gs_beta, es_beta]

    for p in inputs:
        if not p.exists():
            raise FileNotFoundError(p)

    provenance = {
        'doi': DOI,
        'estimator': 'conditional',
        'physical_projection': True,
        'n_permutations': int(args.n_perm),
        'seeds': {
            'gs_perm': 23023,
            'es_perm': 23024,
            'gs_leave_a': 23025,
            'es_leave_a': 23026,
        },
        'source_sha256': {
            p.relative_to(root).as_posix(): sha256(p)
            for p in inputs
        },
    }

    if readme.exists():
        provenance['source_sha256'][readme.relative_to(root).as_posix()] = sha256(readme)

    design = audit_design(gs_raw, es_raw)

    print('RECONSTRUCTING GS CONDITIONAL MI...')
    mi_gs, _ = mi_matrix(gs_raw, estimator='conditional', physical=True)

    print('RECONSTRUCTING ES CONDITIONAL MI...')
    mi_es, _ = mi_matrix(es_raw, estimator='conditional', physical=True)

    gs = [
        p for p in load_beta_profiles(gs_beta, group='databulk')
        if p['kind'] == 'experimental'
    ]

    es = [
        p for p in load_beta_profiles(es_beta, group='databulk')
        if p['kind'] == 'experimental'
    ]

    profile_rows = []
    edge_rows = []
    dist_rows = []

    for state, mi, profiles in [('GS', mi_gs, gs), ('ES', mi_es, es)]:
        for p in profiles:
            sites = (p['sites'] - 1).astype(int)
            invb = np.asarray(p['inv_beta_norm'], float)
            ch = channels(mi, sites)

            row = {
                'state': state,
                'cell_index': int(p['cell_index']),
                'n_sites': int(len(sites)),
            }

            for name in ['all', 'in', 'cut']:
                for mode in ['sum', 'mean']:
                    rp, rs = corr(ch[name][mode], invb)
                    row[name + '_' + mode + '_pearson'] = num(rp)
                    row[name + '_' + mode + '_spearman'] = num(rs)

                rp, rs = corr(ch[name]['count'], invb)
                row[name + '_coverage_pearson'] = num(rp)
                row[name + '_coverage_spearman'] = num(rs)

            profile_rows.append(row)

            if state == 'GS':
                e = edge_control(ch['cut']['mean'], invb)
                e.update({
                    'cell_index': int(p['cell_index']),
                    'n_sites': int(len(sites)),
                })
                edge_rows.append(e)

                d = distance_control(mi, sites, invb)
                d.update({
                    'cell_index': int(p['cell_index']),
                    'n_sites': int(len(sites)),
                })
                dist_rows.append(d)

    print('RUNNING DISTANCE-PRESERVING NULLS...')

    perm = {
        'GS': distance_permutation(mi_gs, gs, args.n_perm, 23023),
        'ES': distance_permutation(mi_es, es, args.n_perm, 23024),
    }

    print('RUNNING LEAVE-A-OUT DISTANCE-MATCHED NULLS...')

    leave = {
        'GS': leave_a_out(mi_gs, gs, args.n_perm, 23025),
        'ES': leave_a_out(mi_es, es, args.n_perm, 23026),
    }

    gsr = [r for r in profile_rows if r['state'] == 'GS']

    means = {}

    for key in [
        'all_sum_pearson',
        'all_sum_spearman',
        'all_mean_pearson',
        'all_mean_spearman',
        'in_sum_pearson',
        'in_sum_spearman',
        'in_mean_pearson',
        'in_mean_spearman',
        'cut_sum_pearson',
        'cut_sum_spearman',
        'cut_mean_pearson',
        'cut_mean_spearman',
    ]:
        vals = [
            np.nan if r[key] is None else r[key]
            for r in gsr
        ]
        means[key] = num(np.nanmean(vals))

    edge_valid = [r for r in edge_rows if r['n_sites'] >= 5]
    dist_valid = [r for r in dist_rows if r['n_sites'] >= 5]

    summary = {
        'gs_five_profile_means': means,
        'gs_edge_partial_L_ge_5': {
            'mean_pearson': num(np.mean([
                r['partial_edge_pearson']
                for r in edge_valid
            ])),
            'mean_spearman': num(np.mean([
                r['partial_edge_spearman']
                for r in edge_valid
            ])),
        },
        'gs_distance_partial_L_ge_5': {
            'mean_pearson': num(np.mean([
                r['partial_distance_pearson']
                for r in dist_valid
            ])),
            'mean_spearman': num(np.mean([
                r['partial_distance_spearman']
                for r in dist_valid
            ])),
        },
    }

    audit = {
        'provenance': provenance,
        'measurement_design': design,
        'profile_correlations': profile_rows,
        'gs_edge_geometry_control': edge_rows,
        'gs_distance_decay_control': dist_rows,
        'distance_preserving_permutation_null': perm,
        'leave_a_out_distance_matched_null': leave,
        'summary': summary,
    }

    jp = outdir / 'audit_summary.json'
    cp = outdir / 'profile_correlations.csv'

    jp.write_text(
        json.dumps(
            audit,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        ) + '\n',
        encoding='utf-8',
    )

    with cp.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(
            f,
            fieldnames=list(profile_rows[0].keys()),
        )
        w.writeheader()
        w.writerows(profile_rows)

    print('AUDIT_COMPLETED')
    print('SUMMARY_JSON=', jp)
    print('PROFILE_CSV=', cp)
    print('SUMMARY_SHA256=', sha256(jp))
    print('PROFILE_CSV_SHA256=', sha256(cp))
    print('DESIGN_COMPLETE_PAIRS=', design['n_complete_pairs'])
    print('DESIGN_NONCOMPLETE_PAIRS=', design['n_noncomplete_pairs'])
    print('MOD5_RULE_EXACT=', design['mod5_rule_exact'])
    print('GS_CUT_SUM_MEAN_PEARSON=', means['cut_sum_pearson'])
    print('GS_CUT_MEAN_MEAN_PEARSON=', means['cut_mean_pearson'])
    print('GS_CUT_MEAN_MEAN_SPEARMAN=', means['cut_mean_spearman'])
    print(
        'GS_EDGE_PARTIAL_MEAN_PEARSON=',
        summary['gs_edge_partial_L_ge_5']['mean_pearson'],
    )
    print(
        'GS_DISTANCE_PARTIAL_MEAN_PEARSON=',
        summary['gs_distance_partial_L_ge_5']['mean_pearson'],
    )

if __name__ == '__main__':
    main()
