using Weave

cd(@__DIR__)
weave("backscatter_energy_depth.ipynb", fig_ext=".svg")
weave("backscatter_mc_browning1991.ipynb", fig_ext=".svg")
weave("backscatter_mc_browning1994.ipynb", fig_ext=".svg")
weave("backscatter_mc.ipynb", fig_ext=".svg")
weave("phirhoz_mc.ipynb", fig_ext=".svg")
#weave("visualize_mc.ipynb", fig_ext=".svg")