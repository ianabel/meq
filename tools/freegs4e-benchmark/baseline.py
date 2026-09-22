"""The reference ladder as a COST AND ACCURACY baseline.

Two things a future run needs and neither is in the .npz: what the reference
COST to produce at each rung, and how converged it actually is.  The second is
the one that decides whether a comparison against it is measuring MEQ or
measuring the reference.
"""
import json, os, re, sys, numpy as np

LAD = sys.argv[1]
rows = []
for n in (129, 257, 513):
    d = os.path.join(LAD, "n%d" % n)
    if not os.path.isdir(d):
        continue
    f = [x for x in os.listdir(d) if x.endswith(".npz")]
    if not f:
        continue
    z = np.load(os.path.join(d, f[0]), allow_pickle=True)
    j = json.load(open(os.path.join(d, "H_limited_circular.json")))
    t = open(os.path.join(LAD, "time-%d.log" % n)).read()
    wall = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", t)
    rss = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", t)
    secs = None
    if wall:
        p = [float(x) for x in wall.group(1).split(":")]
        secs = p[-1] + (p[-2]*60 if len(p) > 1 else 0) + (p[0]*3600 if len(p) > 2 else 0)
    rows.append(dict(
        n=n, points=n*n,
        psi_axis=float(z["psi_axis"]), psi_bndry=float(z["psi_bndry"]),
        picard=j.get("stage2_iterations"), status=j.get("picard_status"),
        wall_s=secs, peak_rss_mb=(int(rss.group(1))/1024.0 if rss else None)))

# RICHARDSON, on the two finest rungs, in the order the successive differences
# actually show rather than an assumed one.
out = dict(case="H_limited_circular", boundary="direct Greens integral, cached",
           rungs=rows)
if len(rows) >= 3:
    a, b, c = rows[-3]["psi_axis"], rows[-2]["psi_axis"], rows[-1]["psi_axis"]
    d1, d2 = a - b, b - c
    if d2 != 0.0:
        order = np.log(abs(d1/d2))/np.log(2.0)
        # THE CORRECTION CARRIES THE SIGN OF ( c - b ) AND NOT OF ( b - c ).
        # This sequence DECREASES toward its limit, so the extrapolate is BELOW
        # the finest rung; the first version added instead of subtracting and
        # put it above. Caught by CLAUDE.md's independently recorded 9.3014e-02,
        # which is what an existing number in the tree is for.
        rich = c + ( c - b )/(2.0**order - 1.0)
        out["order_in_h"] = float(order)
        out["psi_axis_richardson"] = float(rich)
        out["psi_axis_error_at_513"] = float(abs(c - rich)/abs(rich))
    for k in ("psi_bndry",):
        a, b, c = (R[k] for R in rows[-3:])
        d1, d2 = a - b, b - c
        if d2 != 0.0:
            o = np.log(abs(d1/d2))/np.log(2.0)
            out[k + "_richardson"] = float(c + ( c - b )/(2.0**o - 1.0))
            out[k + "_error_at_513"] = float(
                abs(c - out[k + "_richardson"])/abs(out[k + "_richardson"]))
print(json.dumps(out, indent=2))
