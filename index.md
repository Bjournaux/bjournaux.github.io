---
title: Home
---

{% include section.html size="full" %}

{% capture text %}
The UW Planetary Mineral Physics Laboratory studies water, ices and minerals at the extreme conditions of planets and moons to understand the evolution and habitability of icy and ocean worlds, support their exploration, and inform fields from environmental science to biomedical engineering.

{% include button.html link="research" text="Explore our research" icon="fa-solid fa-arrow-right" flip=true %}
{% include button.html link="join" text="Join the lab" style="outline" %}
{% endcapture %}

{%
  include hero.html
  image="images/hero.jpg"
  eyebrow="University of Washington · Earth &amp; Space Sciences"
  title="Mineral Physics of water-rich worlds"
  text=text
%}

{% include section.html %}

{% include section-head.html eyebrow="Research" title="What we study" link="research" link-text="All research" %}

{% capture text %}
We study ice polymorphs, salt hydrates and clathrates with diamond anvil cells, THz Raman spectroscopy and synchrotron X-ray diffraction, and are discovering new crystalline phases, such as [three new NaCl hydrates reported in PNAS](https://www.pnas.org/doi/10.1073/pnas.2217125120) in 2023.
{% endcapture %}
{% include feature.html image="images/research/ice7-crystal.jpg" link="research#mineral-physics" title="Water-rich mineral physics" text=text %}

{% capture text %}
Ultrasonic sound-speed, density and phase-equilibrium measurements feed fundamental thermodynamic representations of liquids and solids, released in our open-source code [SeaFreeze](seafreeze).
{% endcapture %}
{% include feature.html image="images/research/phase-3d.jpg" link="research#thermodynamics" title="Thermodynamics at extreme conditions" text=text flip=true %}

{% capture text %}
Our data constrain interior models through the [PlanetProfile](https://github.com/vancesteven/PlanetProfile) project, advancing Titan's [geodynamics](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021GL097602) and [seismic response](https://iopscience.iop.org/article/10.3847/PSJ/ac787e) and [exoplanet habitability](https://www.nature.com/articles/s41467-022-35187-4).
{% endcapture %}
{% include feature.html image="images/research/interior.jpg" link="research#geophysical-modeling" title="Icy moons & ocean exoplanets" text=text %}

{% capture text %}
Our measurements support the [NASA Dragonfly](https://dragonfly.jhuapl.edu/) seismometer, Europa Clipper and JUICE observations, and organ cryopreservation research with the [Public Thermo Lab](https://publicthermo.com/) at Texas A&M.
{% endcapture %}
{% include feature.html image="images/research/dragonfly-titan-dunes.jpg" link="research#applications" title="From space missions to cryobiology" text=text flip=true %}

{% include section.html size="full" %}

{% capture text %}
Thermodynamic and elastic properties of water, ice polymorphs (Ih to VII/X) up to 100 GPa and 10,000 K, and aqueous NaCl up to 8 GPa and 2,000 K. Python and MATLAB packages, and an online app.

{% include button.html link="https://github.com/Bjournaux/SeaFreeze" text="Get SeaFreeze" icon="fa-brands fa-github" %}
{% include button.html link="https://seafreeze.streamlit.app/" text="Try it online" icon="fa-solid fa-play" style="highlight" %}
{% include button.html link="seafreeze" text="Learn more" style="outline" %}
{% endcapture %}

{%
  include band.html
  image="images/home/clipper.jpg"
  eyebrow="Open-source software"
  title="SeaFreeze: thermodynamics of water and ices"
  text=text
%}

{% include section.html %}

{% include section-head.html eyebrow="News" title="Latest from the lab" link="news" link-text="All news" %}

{% include news-cards.html count=3 %}

{% include section.html %}

{% include section-head.html eyebrow="Publications" title="Featured papers" link="publications" link-text="All publications" %}

{% include featured-papers.html %}

{% include section.html %}

{% include section-head.html eyebrow="People" title="The team" link="people" link-text="Meet everyone" %}

{% include list.html data="members" component="portrait" filter="role == 'principal-investigator'" %}
{% include list.html data="members" component="portrait" filter="role == 'grad-student'" %}

{% include section.html size="full" %}

{% capture text %}
We are looking for talented candidates at the postdoctoral level with a background in mineral physics, physical chemistry, planetary sciences, geophysics and astrobiology. Postdoctoral opportunities are available through the [51 Pegasi b Fellowship](https://www.hsfoundation.org/programs/science/51-pegasi-b-fellowship/) of the Heising-Simons Foundation.

{% include button.html link="join" text="How to join" %}
{% endcapture %}

{%
  include band.html
  image="images/home/eclipse.jpg"
  position="right center"
  eyebrow="Opportunities"
  title="Join us"
  text=text
%}
