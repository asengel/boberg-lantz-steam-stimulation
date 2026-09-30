"""Run Green & Willhite Example 8.6 with the Boberg-Lantz model."""

from pathlib import Path
import sys

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from boberg_lantz import ModelInputs, run_model


def main():
    inputs = ModelInputs()
    results, geometry = run_model(inputs)

    print("\nHeated-zone geometry")
    print("--------------------")
    print(f"Hs   = {geometry['H_s']:.3f} Btu/lbm")
    print(f"tD   = {geometry['tD_injection']:.6f}")
    print(f"G(tD)= {geometry['G_tD']:.6f}")
    print(f"Ah   = {geometry['A_h']:.2f} ft^2")
    print(f"rh   = {geometry['r_h']:.3f} ft")
    print(f"z    = {geometry['z']:.3f} ft")

    print("\nProduction response")
    print("-------------------")
    print(results.to_string(index=False))

    results_dir = ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    results.to_csv(results_dir / "example_8_6_results.csv", index=False)

    # Oil production rate
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(results["t - ti (days)"], results["qo (STB/D)"], marker="o")
    ax.set_xlabel("Time after steam injection, days")
    ax.set_ylabel("Oil production rate, STB/D")
    ax.set_title("Boberg-Lantz production response")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(results_dir / "production_rate.png", dpi=180)
    plt.close(fig)

    # Produced-fluid temperature
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(results["t - ti (days)"], results["Tp (F)"], marker="o")
    ax.set_xlabel("Time after steam injection, days")
    ax.set_ylabel("Produced-fluid temperature, degF")
    ax.set_title("Produced-fluid temperature decline")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(results_dir / "produced_temperature.png", dpi=180)
    plt.close(fig)

    # Hot-oil viscosity
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(results["t - ti (days)"], results["mu_oh (cp)"], marker="o")
    ax.set_xlabel("Time after steam injection, days")
    ax.set_ylabel("Hot-oil viscosity, cp")
    ax.set_title("Oil-viscosity recovery during cooling")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(results_dir / "oil_viscosity.png", dpi=180)
    plt.close(fig)

    # Water-oil ratio
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.plot(results["t - ti (days)"], results["Fwo (bbl/bbl)"], marker="o")
    ax.set_xlabel("Time after steam injection, days")
    ax.set_ylabel("Water-oil ratio, bbl/bbl")
    ax.set_title("Water-oil ratio after steam stimulation")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(results_dir / "water_oil_ratio.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
