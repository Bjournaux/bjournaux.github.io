---
title: People
nav:
  order: 2
  tooltip: Meet the team
---

# People

{% include section.html %}

## Principal investigator

{% include list.html data="members" component="portrait" filter="role == 'principal-investigator'" %}

## Graduate students

{% include list.html data="members" component="portrait" filter="role == 'grad-student'" %}

## Undergraduate researchers

{% include list.html data="members" component="portrait" filter="role == 'undergrad'" %}

{% include section.html %}

## Alumni

### Former graduate students

{% include alumni-list.html level="graduate" %}

### Former undergraduate researchers

{% include alumni-list.html level="undergraduate" %}

{% include section.html size="full" %}

{% capture text %}
We are looking for talented candidates at the postdoctoral level with a background in mineral physics, physical chemistry, planetary sciences, geophysics and astrobiology.

{% include button.html link="join" text="How to join" %}
{% endcapture %}

{% include band.html image="images/home/eclipse.jpg" position="right center" eyebrow="Opportunities" title="Join us" text=text %}
