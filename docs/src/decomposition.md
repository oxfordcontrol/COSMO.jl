# Chordal Decomposition

For very large sparse SDPs it often makes sense to analyse the sparsity structure of the PSD constraints. In many cases one large PSD constraint can be decomposed into several smaller constraints which can significantly reduce the solver time.

Please note that this feature is still experimental and not active by default. If you want to activate this feature adjust the settings of COSMO accordingly:

```julia
settings = COSMO.Settings()
settings.decompose = true
```