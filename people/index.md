---
title: People
nav:
  order: 2
  tooltip: Meet the team
---

# People

{% include section.html %}

## Principal Investigator

{% include list.html data="members" component="portrait" filter="role == 'principal-investigator'" %}

## Graduate Students

{% include list.html data="members" component="portrait" filter="role == 'grad-student'" %}

## Undergraduate Researchers

{% include list.html data="members" component="portrait" filter="role == 'undergrad'" %}

{% include section.html %}

## Alumni

### Former Graduate Students

{% include alumni-list.html level="graduate" %}

### Former Undergraduate Researchers

{% include alumni-list.html level="undergraduate" %}

{% include section.html size="full" %}

{% capture text %}
We are looking for talented candidates at the postdoctoral level with a background in mineral physics, physical chemistry, planetary sciences, geophysics and astrobiology.

{% include button.html link="join" text="How to join" %}
{% endcapture %}

{% include band.html image="images/home/eclipse.jpg" position="right center" eyebrow="Opportunities" title="Join us" text=text %}
