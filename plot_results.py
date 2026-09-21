import json
import numpy as np
import matplotlib.pyplot as plt

# Load JSON
with open("resultados/sweep-power_legacy.json", "r") as f:
    data = json.load(f)

powers = []
ber_x = []
ber_y = []

for point in data["points"]:
    p = point["power"]
    powers.append(10 * np.log10(p * 1000)) # to dBm
    b_x = point["ber"][0][0][0]
    b_y = point["ber"][0][0][1]
    ber_x.append(b_x)
    ber_y.append(b_y)

plt.figure(figsize=(8, 6))
plt.semilogy(powers, ber_x, 'o-', label='Pol X')
plt.semilogy(powers, ber_y, 's-', label='Pol Y')
plt.xlabel('Launch Power (dBm)')
plt.ylabel('BER')
plt.title('BER vs Launch Power')
plt.grid(True, which="both", ls="--")
plt.legend()
plt.savefig('ber_plot.png', dpi=150)
plt.close()

# Load NPZ
npz = np.load("resultados/sweep-power_legacy.npz")
cx = npz["constellation_x"]
cy = npz["constellation_y"]

# Filter out initial samples and plot constellation
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
subset_cx = cx[-8000:] 
plt.scatter(np.real(subset_cx), np.imag(subset_cx), s=1, alpha=0.5, c='b')
plt.title("Constellation - Pol X (10 dBm)")
plt.xlabel("In-phase")
plt.ylabel("Quadrature")
plt.axis('square')
plt.grid(True, ls="--")

plt.subplot(1, 2, 2)
subset_cy = cy[-8000:]
plt.scatter(np.real(subset_cy), np.imag(subset_cy), s=1, alpha=0.5, c='r')
plt.title("Constellation - Pol Y (10 dBm)")
plt.xlabel("In-phase")
plt.ylabel("Quadrature")
plt.axis('square')
plt.grid(True, ls="--")

plt.tight_layout()
plt.savefig('constellation_plot.png', dpi=150)
plt.close()
