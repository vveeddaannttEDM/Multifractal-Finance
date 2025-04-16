import numpy as np

def compute_mfdfa(series, scale_range, q_vals):
    N = len(series)
    Fq = np.zeros((len(q_vals), len(scale_range)))

    for i, s in enumerate(scale_range):
        n_segments = N // s
        F2 = []

        for j in range(n_segments):
            seg = series[j * s:(j + 1) * s]
            poly = np.polyfit(range(s), seg, 1)
            trend = np.polyval(poly, range(s))
            F2.append(np.mean((seg - trend)**2))

        F2 = np.array(F2)
        for k, q in enumerate(q_vals):
            if q == 0:
                Fq[k, i] = np.exp(0.5 * np.mean(np.log(F2)))
            else:
                Fq[k, i] = (np.mean(F2**(q / 2)))**(1 / q)

    return Fq
