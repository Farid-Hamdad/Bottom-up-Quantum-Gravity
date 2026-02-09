# bottom-up-quantum-gravity

**Emergent Space, Time and Gravity from Quantum Entanglement**  
Farid Hamdad — 2026

---

## 🧠 Core Idea (One Sentence)

**Space, time and gravity are not fundamental.  
They emerge collectively from the entanglement structure of a global quantum state.**

---

## ❓ Why This Project?

Modern physics successfully describes:

- **Quantum physics** (quantum field theories),
- **Classical gravity** (general relativity),

yet leaves a deep question unanswered:

> **Why does spacetime exist at all, and why does gravity have a geometric and thermodynamic nature?**

This project explores a minimal but radical hypothesis:

> **Spacetime is not the stage of physics,  
> but an emergent object reconstructed from quantum entanglement.**

---

## 🔑 Minimal Postulate

The entire framework rests on a single assumption:

- There exists a **global pure quantum state** \( |\Psi\rangle \),
- defined on a **finite set of degrees of freedom** (qubits),
- with **no predefined space, no time, no background metric**.

Everything else —  
**time, space, dimension, geometry, gravity** —  
must emerge **solely from the internal structure of entanglement in \( |\Psi\rangle \)**.

---

## ⏱️ How Time Emerges

Time is **not an external parameter**.

It is identified with the **modular flow** associated with reduced density matrices:

\[
K_A = -\log \rho_A
\]

- Modular flow defines an intrinsic notion of evolution,
- determined by correlations between subsystems.

👉 **Time is interpreted as a relational, information-theoretic notion**,  
linked to relative decoherence between subsystems.

---

## 📐 How Space Emerges

1. Compute **mutual information** \( I(i:j) \) between qubit pairs.
2. Define an **information-theoretic distance**:

\[
d_{ij} = -\log\left(\frac{I(i:j)}{I_{\max}}\right)
\]

3. Reconstruct an effective geometry using embedding techniques  
   (graph reconstruction, MDS, spectral methods).

👉 **Space is the geometric map of quantum correlations.**

---

## 🧩 Key Result #1 — Dimension Is Emergent

### ▶️ Case: **N = 9 qubits**

- **Locally entangled states**
  - Stable reconstruction of a **2D spatial geometry**,
  - Spectral dimension consistent with \( d \approx 2 \).

- **Non-local entanglement patterns**
  - Forced emergence of an **additional effective dimension**,
  - Clear increase in spectral dimension.

👉 **Spatial dimension is not postulated — it is imposed by entanglement.**

---

### ▶️ Extension: **N = 16 qubits (New)**

- Improved numerical stability,
- Reduced finite-size artifacts,
- Clearer and smoother dimensional transitions.

Key observation:

- Dimension depends on **entanglement structure**,  
  **not** on the number of qubits alone.

👉 **Dimensional emergence persists and strengthens with increasing system size.**

---

## 🌀 Key Result #2 — ER = EPR Becomes Measurable

The **ER = EPR** correspondence is tested operationally:

- Two qubits are **strongly entangled**,
- but **topologically distant** in the interaction graph,
- yet become **geometrically close** in the reconstructed emergent space.

👉 **Maximal entanglement manifests as a geodesic shortcut**  
(*wormhole-like signature*).

This effect is observed for **N = 9** and confirmed for **N = 16**.

---

## 🌡️ Key Result #3 — Gravity as Thermodynamics

A Jacobson-type relation is tested numerically:

\[
\delta S \simeq \beta \, \delta E
\]

Results show:

- Variations of entanglement entropy correlate with
- an effective energy variation defined from correlation reorganization.

👉 **Gravity emerges as a thermodynamic equation of state**,  
not as a fundamental force.

---

## 🚫 What This Work Does NOT Claim

To be explicit:

- ❌ Not a full derivation of General Relativity,
- ❌ Not a final cosmological model,
- ❌ Not a standard unification framework.

This is a **proof of principle**:

- Conceptual,
- Numerical,
- Reproducible,

demonstrating how **space, time and gravity can emerge from a single quantum substrate**.

---

## 🎯 Target Audience

- MSc / PhD students in physics
- Theoretical physicists
- Quantum information researchers
- Scientists interested in emergence and information

👉 **No prior expertise in quantum gravity is required.**

---

## 📁 Repository Structure

```text
paper/      → Main article (LaTeX + PDF)
figures/    → Key figures
appendix/   → Extended analyses
code/       → Geometry reconstruction scripts
README.md   → This document
