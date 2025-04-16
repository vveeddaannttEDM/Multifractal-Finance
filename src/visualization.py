import matplotlib.pyplot as plt

def plot_multifractal_spectrum(q_vals, Hq):
    tq = q_vals * Hq - 1
    alpha = np.diff(tq) / np.diff(q_vals)
    f_alpha = q_vals[:-1] * alpha - tq[:-1]

    plt.plot(alpha, f_alpha)
    plt.xlabel("α (Holder exponent)")
    plt.ylabel("f(α)")
    plt.title("Multifractal Spectrum")
    plt.grid(True)
    plt.show()
