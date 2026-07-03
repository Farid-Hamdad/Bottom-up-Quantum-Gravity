# Paper 16 — Roadmap

## Objective

Paper 16 aims to turn the numerical spectral evidence of Paper 15 into a controlled continuum statement:

\[
c_N\frac{D_N-W_N}{\epsilon_N}
\longrightarrow
-\Delta_g.
\]

## Current numerical results

### Circle \(S^1\)

\[
c_NL_N\to-\Delta_{S^1}
\]

At \(N=1024\):

- modes compared: 12
- mean relative error: 0.0051
- median relative error: 0.0051
- max relative error: 0.0107
- \(\lambda_1^{scaled}=1.0062\)
- \(\lambda_1^{target}=1.0000\)

### Flat torus \(T^2\)

\[
c_NL_N\to-\Delta_{T^2}
\]

At \(N=1024\):

- modes compared: 20
- mean relative error: 0.0266
- median relative error: 0.0263
- max relative error: 0.0522
- \(\lambda_1^{scaled}=40.9613\)
- \(\lambda_1^{target}=39.4784\)

## Next steps

1. Derive the normalization \(c_N\).
2. Test larger \(N\) using sparse eigensolvers.
3. Test eigenfunction convergence.
4. Test the sphere spectrum:
   \[
   \lambda_\ell=\ell(\ell+1).
   \]
5. Replace synthetic kernels by true mutual-information graphs.
6. Connect spectral convergence to the BuP spectral action.
