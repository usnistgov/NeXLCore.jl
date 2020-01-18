# Base system
Pkg.activate()
Pkg.gc()
Pkg.update()
Pkg.precompile()
# NeXLUncertainties
Pkg.activate(joinpath(homedir(),".julia/dev/NeXLUncertainties"))
Pkg.update()
Pkg.precompile()
# NeXLCore
Pkg.activate(joinpath(homedir(),".julia/dev/NeXLCore"))
Pkg.update()
Pkg.precompile()
# NeXLMatrixCorrection
Pkg.activate(joinpath(homedir(),".julia/dev/NeXLMatrixCorrection"))
Pkg.update()
Pkg.precompile()
# NeXLSpectrum
Pkg.activate(joinpath(homedir(),".julia/dev/NeXLSpectrum"))
Pkg.update()
Pkg.precompile()
# NeXLParticle
Pkg.activate(joinpath(homedir(),".julia/dev/NeXLParticle"))
Pkg.update()
Pkg.precompile()

Pkg.activate()
