"""Before any paired test: can their reported per-fold C-index be recomputed from the pickle?

If not, my reading of `risk` (sign) or `censorship` (which value means event) is wrong, and every
number downstream is wrong with it. The published per-fold values from this run's own log are
0.5875, 0.7834, 0.5328, 0.5222, 0.6475.
"""
import pickle, glob, numpy as np, os, sys

def cindex(risk, time, event):
    ii, jj = [], []
    for i in np.flatnonzero(event == 1):
        later = np.flatnonzero(time > time[i])
        if later.size:
            ii.append(np.full(later.size, i)); jj.append(later)
    if not ii: return float("nan")
    ii, jj = np.concatenate(ii), np.concatenate(jj)
    a, b = risk[ii], risk[jj]
    return float((np.sum(a > b) + 0.5 * np.sum(a == b)) / ii.size)

D = sys.argv[1]
logged = [0.5875, 0.7834, 0.5328, 0.5222, 0.6475]
print("%-6s %8s %10s %10s %10s" % ("fold", "n", "events", "from pkl", "logged"))
allc = {}
for k, f in enumerate(sorted(glob.glob(D + "/split_*_results.pkl"))):
    d = pickle.load(open(f, "rb"))
    cid = list(d)
    t = np.array([d[c]["time"] for c in cid], float)
    r = np.array([float(d[c]["risk"]) for c in cid])
    # SurvPath convention, stated in its own code: censorship 0 = event observed
    e = np.array([1.0 if float(d[c]["censorship"]) == 0.0 else 0.0 for c in cid])
    got = cindex(r, t, e)
    print("%-6d %8d %10d %10.4f %10.4f" % (k, len(cid), int(e.sum()), got, logged[k]))
    for c in cid:
        allc[c] = (float(d[c]["risk"]), float(d[c]["time"]),
                   1.0 if float(d[c]["censorship"]) == 0.0 else 0.0, k)
print()
print("cases pooled: %d unique" % len(allc))
np.save(os.path.join(D, "survpath_percase.npy"), allc, allow_pickle=True)
print("saved per-case risks ->", os.path.join(D, "survpath_percase.npy"))
