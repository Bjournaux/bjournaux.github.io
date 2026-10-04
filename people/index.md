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

{% include alumni-list.html %}

{% include section.html size="full" %}

{% capture text %}
We are looking for talented graduate students and postdocs.

{% include button.html link="join" text="How to join" %}
{% endcapture %}

{% include band.html image="images/home/eclipse.jpg" position="right center" eyebrow="Opportunities" title="Join us" text=text %}
