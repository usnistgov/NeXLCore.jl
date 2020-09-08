using Colors
using FileIO

# All these palettes have 256 colors to simplify saving as 8-bit palettized PNG. (Maybe...)

"""
A colorblind friendly palette suited to use with log base-10 transformed data on the range [0.0, 1.0].
"""
const Log3BandColorblind = reduce(
    append!,
    (
        range(colorant"rgb(0%,38%,100%)", stop = colorant"white", length = 120)[1:85],
        range(colorant"rgb(75%,0%,0%)", stop = colorant"white", length = 120)[1:85],
        range(colorant"rgb(100%,90%,0%)", stop = colorant"black", length = 85),
        [colorant"yellow"],
    ),
)

"""
Dave Bright's palette suited to use with log base-10 transformed data on the range [0.0, 1.0].
The 256ᵗʰ entry is yellow for an error condition.
"""
const Log3BandBright = reduce(
    append!,
    (
        range(colorant"rgb(60,0,0)", stop = colorant"rgb(255,200,200)", length = 85), # 60–255/0–200/0–200
        range(colorant"rgb(0,60,0)", stop = colorant"rgb(200,255,200)", length = 85), # 0–200/60–255/0–200
        range(colorant"rgb(0,0,60)", stop = colorant"rgb(200,200,255)", length = 85), # 0–200/0–200/60–255
        [colorant"yellow"],
    ),
)

"""
An RGB gray-scale palette with entry 256 being yellow.
"""
const GrayScale =
    reduce(append!, (range(colorant"rgb(255,255,255)", stop = colorant"rgb(0,0,0)", length = 255), [colorant"yellow"]))


"""
Tranforms numbers on the range [1.0e-3, 1.0] into a colorblind friendly palette using a log base-10 transform.
NaNs are plotted in yellow.
"""
Log3BandC(f::AbstractFloat)::Colorant = Log3BandColorblind[isnan(f) ? 256 :
                       min(max(trunc(Int, 255.0 * (log(10.0, max(f, 1.0e-3)) + 3.0) / 3.0), 0), 254) + 1]

"""
Tranforms numbers on the range [1.0e-3, 1.0] into David Bright's Log3-band palette using a log base-10 transform.
NaNs are plotted in yellow.
"""
Log3Band(f::AbstractFloat)::Colorant =
    Log3BandBright[isnan(f) ? 256 : min(max(trunc(Int, 255.0 * (log(10.0, max(f, 1.0e-3)) + 3.0) / 3.0), 0), 254) + 1]

"""
Tranforms numbers on the range [1.0e-3, 1.0] onto a Log base-10 gray scale palette. NaNs are plotted in yellow.
"""
LogScale(f::AbstractFloat)::Colorant =
    GrayScale[isnan(f) ? 256 : min(max(trunc(Int, 255.0 * (log(10.0, max(f, 1.0e-3)) + 3.0) / 3.0), 0), 254) + 1]

"""
Tranforms numbers on the range [0.0, 1.0] onto a linear gray scale palette. NaNs are plotted in yellow.
"""
LinearScale(f::AbstractFloat)::Colorant = GrayScale[isnan(f) ? 256 : min(254, max(0, trunc(Int, 254.0 * f))) + 1]

"""
Loads a legend image from the package source directory. Some legends include "LinearScale.png", "LogScale.png",
"Log3BandBright.png", "Log3BancColorblind.png".
"""
loadlegend(fn::String) = FileIO.load(joinpath(dirname(pathof(@__MODULE__)),fn))
