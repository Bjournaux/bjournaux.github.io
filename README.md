Visit **[bjournaux.github.io](https://bjournaux.github.io)** 🚀

# UW Planetary Mineral Physics Laboratory website

Built with the [Lab Website Template](https://greene-lab.gitbook.io/lab-website-template-docs) (Jekyll + GitHub Pages).
Every change pushed to the `main` branch rebuilds and republishes the site automatically, in about 2 minutes.

## Editing without installing anything

On github.com, open a file, click the **pencil icon**, edit, then **Commit changes**.
To add an image, open its folder and use **Add file → Upload files**.

| I want to…                         | Edit this                                   |
| ---------------------------------- | ------------------------------------------- |
| Add or update a lab member         | `_members/firstname-lastname.md` + photo in `images/people/` |
| Move someone to the alumni list    | delete their `_members/` file, add them to `_data/alumni.yaml` |
| Update the PI's CV (PDF)           | replace `files/CV_Baptiste_Journaux.pdf`    |
| Add gallery photos                 | photo in `images/gallery/lab`, `field` or `sky` + one line in `gallery/index.md` |
| Post a news item                   | `_posts/YYYY-MM-DD-short-title.md`          |
| Add a publication                  | `_data/sources.yaml`                        |
| Change homepage text               | `index.md`                                  |
| Change another page                | `research/`, `people/`, `publications/`, `news/`, `seafreeze/`, `gallery/`, `join/` → `index.md` |
| Change colors or fonts             | `_styles/-theme.scss`                       |
| Change layout details              | `_styles/zz-lab-design.scss`                |
| Footer links (email, Scholar…)     | `_config.yaml` → `links:`                   |
| Logo / browser-tab icon            | `images/logo.svg`, `images/icon.svg`        |

### Add a lab member

Create `_members/jane-doe.md`, and put a square-ish photo at `images/people/jane-doe.jpg`:

```markdown
---
name: Jane Doe
image: images/people/jane-doe.jpg
role: grad-student        # principal-investigator, postdoc, grad-student, undergrad
links:
  email: jdoe@uw.edu
  orcid: 0000-0000-0000-0000
---

One or two sentences about Jane's research.
```

When someone leaves, delete their file in `_members/` and add a few lines to `_data/alumni.yaml` (name, role, years, and where they went).
New members without a photo yet can use an initials placeholder: just replace the image file later, keeping the same name.

### Post a news item

Create `_posts/2026-10-15-new-paper-in-nature.md` (the date in the file name is the post's date):

```markdown
---
title: New paper in Nature
icon: fa-solid fa-newspaper    # any Font Awesome icon: fa-solid fa-trophy, fa-solid fa-microphone…
---

The first paragraph is shown on the homepage card. Links work like [this](https://example.com).
```

If only the year or month is known, add `date-precision: year` (or `month`) so no fake day is shown.
The 3 most recent posts appear on the homepage automatically.

### Add a publication

Add an entry at the top of the right section in `_data/sources.yaml`:

```yaml
- id: doi:10.1038/s41586-025-09818-x
  link: https://www.nature.com/articles/s41586-025-09818-x   # optional
  type: paper                                               # paper, preprint or book
```

Title, authors, journal and year are filled in automatically from the DOI when the site rebuilds.
Add `group: featured` and an `image:` to show it on the homepage.
Entries without a DOI (e.g. under review) need `title`, `authors`, `publisher` and `type` written out.

## Preview on your own computer (optional)

Needs Ruby from Homebrew (`brew install ruby`). Then, in this folder:

```bash
./preview.sh
```

and open <http://localhost:4000>. Changes appear as you save. Press Ctrl+C to stop.
To regenerate citations locally after editing `sources.yaml`: `pip install -r _cite/requirements.txt`, then `python _cite/cite.py`.

## Notes

- Dark mode is the default (`_scripts/dark-mode.js`); visitors can switch with the toggle in the footer.
- `Gemfile.local` and `preview.sh` are only for local previews; GitHub uses `Gemfile` / `Gemfile.lock`.
