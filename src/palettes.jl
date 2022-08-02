using Colors
using FileIO

# All these palettes have 256 colors to simplify saving as 8-bit palettized PNG.

"""
A colorblind friendly palette suited to use with log base-10 transformed data on the range [0.0, 1.0].
"""
const Log3BandColorblind = [
    colorant"rgb(0,0,0)",
    range(colorant"rgb(0%,38%,100%)", stop=colorant"white", length=120)[84:-1:1]...,
    range(colorant"rgb(75%,0%,0%)", stop=colorant"white", length=120)[85:-1:1]...,
    range(colorant"rgb(100%,90%,0%)", stop=colorant"black", length=85)...,
    colorant"yellow"
]

"""
Dave Bright's palette suited to use with log base-10 transformed data on the range [0.0, 1.0].
The 256ᵗʰ entry is yellow for an error condition.
"""
const Log3BandBright = [
    colorant"rgb(0,0,0)",
    range(colorant"rgb(0,0,60)", stop=colorant"rgb(200,200,255)", length=84)..., # 0–200/0–200/60–255
    range(colorant"rgb(0,60,0)", stop=colorant"rgb(200,255,200)", length=85)..., # 0–200/60–255/0–200
    range(colorant"rgb(60,0,0)", stop=colorant"rgb(255,200,200)", length=85)..., # 60–255/0–200/0–200
    colorant"yellow"
]

"""
An RGB gray-scale palette with entry 256 being yellow.
"""
const GrayScale = [
    range(colorant"rgb(0,0,0)", stop=colorant"rgb(255,255,255)", length=255)...,
    colorant"yellow"
]

linearscale(f::AbstractFloat)::Int = isnan(f) ? 256 : trunc(Int, 254.0 * clamp(f, 0.0, 1.0)) + 1
log3scale(f::AbstractFloat)::Int = isnan(f) ? 256 : trunc(Int, (254.0 / 3.0) * (log(10.0, clamp(f, 1.0e-3, 1.0)) + 3.0)) + 1

"""
Tranforms numbers on the range [1.0e-3, 1.0] into a colorblind friendly palette using a log base-10 transform.
NaNs are plotted in yellow.
"""
Log3BandC(f::AbstractFloat) = Log3BandColorblind[log3scale(f)]

"""
Tranforms numbers on the range [1.0e-3, 1.0] into David Bright's Log3-band palette using a log base-10 transform.
NaNs are plotted in yellow.
"""
Log3Band(f::AbstractFloat) = Log3BandBright[log3scale(f)]

"""
Tranforms numbers on the range [1.0e-3, 1.0] onto a Log base-10 gray scale palette. NaNs are plotted in yellow.
"""
LogScale(f::AbstractFloat) = GrayScale[log3scale(f)]

"""
Tranforms numbers on the range [0.0, 1.0] onto a linear gray scale palette. NaNs are plotted in yellow.
"""
LinearScale(f::AbstractFloat) = GrayScale[linearscale(f)]

"""
Loads a legend image from the package source directory. Some legends include "LinearScale.png", "LogScale.png",
"Log3BandBright.png", "Log3BandColorblind.png".
"""
loadlegend(fn::String) = FileIO.load(joinpath(@__DIR__, "..", "resources", fn))