---
title: SeaFreeze
nav:
  order: 5
  tooltip: Open-source thermodynamics of water and ices
---

{% include section.html size="full" %}

{% capture text %}
Thermodynamic and elastic properties of pure water, ice polymorphs (Ih, II, III, V, VI and VII/X) up to 100 GPa and 10,000 K, and aqueous NaCl solutions up to 8 GPa and 2,000 K: the conditions found in all the hydrospheres of our solar system, and beyond.

{% include button.html link="https://github.com/Bjournaux/SeaFreeze" text="SeaFreeze on GitHub" icon="fa-brands fa-github" %}
{% include button.html link="https://seafreeze.streamlit.app/" text="Try it online (beta)" icon="fa-solid fa-play" style="highlight" %}
{% endcapture %}

{% include band.html image="images/hero.jpg" position="center 70%" eyebrow="Open-source software" title="SeaFreeze" text=text %}

{% include section.html %}

## What it does

SeaFreeze evaluates Gibbs energy representations built with the [Local Basis Function](https://github.com/jmichaelb/LocalBasisFunction) approach developed by Dr. J. Michael Brown (University of Washington), constructed to reproduce thermodynamic measurements.
From a pressure and temperature (and concentration for solutions), it returns the thermodynamic and elastic properties of each phase, and the stable phase of the water phase diagram.

The formalism is described in [Journaux et al. (2020)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JE006176), and the liquid water representation in [Bollengier, Brown and Shaw (2019)](https://aip.scitation.org/doi/abs/10.1063/1.5097179).
SeaFreeze is open source (GPL-3.0) and available as Python and MATLAB packages.
You can also try it directly in your browser with the [SeaFreeze online app](https://seafreeze.streamlit.app/) (beta; for research use, please install the Python or MATLAB package).

{% include section.html %}

## Explore water and ices in 3D

Density, sound speed, heat capacity, thermal expansivity and bulk modulus of liquid water and ices across the phase diagram, as computed by SeaFreeze.

{% include seafreeze-3d.html %}

{% include section.html %}

## Phase diagram

{% include figure.html image="images/seafreeze/phase-diagram.png" link="images/seafreeze/phase-diagram.png" caption="Phase transitions predicted by SeaFreeze compared with experimental data. Triple point predictions are also indicated. Click to enlarge." %}

{% include section.html %}

## Credits

SeaFreeze is developed by Baptiste Journaux (main developer) with J. Michael Brown.
It is based on the [Local Basis Function (LBF)](https://github.com/jmichaelb/LocalBasisFunction) code developed by J. Michael Brown (University of Washington).
Penny Espinoza and Ula Jones co-developed the first Python version, with Penny Espinoza also contributing to the equation-of-state development.
SeaFreeze was developed with support from the NASA Astrobiology Institute (Icy Worlds and Titan nodes), the NASA Postdoctoral Program and the NASA Solar System Workings program.

If you use SeaFreeze, please cite [Journaux et al. (2020)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019JE006176) and the references given in the [README](https://github.com/Bjournaux/SeaFreeze#readme).
