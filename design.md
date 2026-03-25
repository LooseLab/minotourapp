Design System Strategy: Obsidian & Ether

## 1. Overview & Creative North Star

**Creative North Star: The Digital Curator**  
This design system is built for high-fidelity genomic analysis, where clarity, precision, and professional focus are paramount. It moves away from generic “tech-blue” aesthetics toward a sophisticated, editorial-inspired interface that feels both scientific and premium.

## 2. Visual Language

**Color palette:** Deep “Obsidian” greens (`#004B44`) serve as the primary brand anchor, providing a sense of depth and stability. “Ether” neutrals and high-clarity whites ensure maximum readability and a clean, laboratory-grade feel.

**Typography:** Manrope is used for its modern, geometric properties and excellent legibility across dense data. Tight tracking and bold weights are used in headers to create a strong visual hierarchy.

**Shape & form:** Round four (4px corner radius) provides a subtle, professional softness without feeling “bubbly.” Card-based layouts with soft shadows create clear separation of concerns.

**Spacing:** A generous, grid-based spacing system ensures that even complex genomic data has “room to breathe,” reducing cognitive load for researchers.

## 3. Core components (Tailwind CSS)

Design tokens below are expressed in Tailwind-style utility language. In this Django app, **logged-in data pages** often implement the same tokens with **scoped plain CSS** (see §5) so we do not load Tailwind’s preflight alongside AdminLTE/Bootstrap.

### Navigation (header & footer)

- **Top app bar:** `bg-white/80` `dark:bg-emerald-950/80` `backdrop-blur-md` `border-b` `border-emerald-50/50`. Clean horizontal navigation with high-contrast active states using `border-b-2` `border-emerald-900`.
- **Footer:** `bg-emerald-50` `dark:bg-emerald-950` `py-8` `border-t` `border-emerald-200/20`. Multi-column link layout with `text-emerald-700/60` and a clear copyright line.

**minoTour implementation note:** The live app shell uses a dark Obsidian gradient nav (`obsidian-chrome.css`) for brand continuity with existing pages. A lighter “registry” style (white TopAppBar) is approximated on the **flowcells** content area (breadcrumbs, toolbar, alerts) per the reference mock.

### Main content & data

- **Page header:** Large, bold titles using `text-5xl` `font-extrabold` `tracking-tighter` `text-[#00322d]`.
- **Data table (flowcells registry):**
  - **Container:** `bg-white` `rounded-lg` `shadow-sm` `border` `border-slate-200/50` `overflow-hidden`.
  - **Typography:** Small, uppercase labels for column headers to maximise data density.
  - **Status badges:** `rounded-full` `px-3` `py-1` `text-[10px]` `font-bold` `uppercase` with semantic colouring (e.g. `bg-emerald-100` `text-emerald-800` for online).
- **Action buttons:**
  - **Primary:** `bg-[#004B44]` `hover:bg-[#00322d]` `text-white` `font-bold` `py-2` `px-4` `rounded-md` `shadow-sm` `transition-all`.
  - **Secondary / ghost:** `bg-slate-100` `text-slate-600` `hover:bg-slate-200`.

### Information architecture (dashboards)

- **Status alerts:** Subtle banner-style alerts using `bg-emerald-50` `border-l-4` `border-emerald-500` for system status updates.
- **Analytical cards:** Feature sections (e.g. “Capacity forecasting”) may use high-quality background imagery with glassmorphism overlays and clear primary CTAs.

## 4. User experience principles

**Clarity first:** Every element must serve a functional purpose. Remove non-essential visual noise.

**Immediate feedback:** Interactive elements (buttons, links) use subtle scale or colour shifts to acknowledge user intent.

**Professional trust:** A refined, non-standard palette (deep greens) signals that minoTour is a specialised tool for serious research.

## 5. Visual reference (flowcells registry mock)

A high-fidelity reference for the **flowcells index** layout (breadcrumbs, toolbar, status strip, dense table, badges, insight cards) is kept in-repo as:

- `docs/design-reference-flowcells-registry.png`

Use it alongside §3 when evolving the flowcells page. Product names in the mock (“GenomicCurator”) are illustrative; minoTour branding applies in code.

## 6. Implementation notes (codebase)

**Login (`/login`)**  
- Template: `web/templates/registration/login.html` extends `web/templates/web/template_login.html`.  
- Uses Tailwind via CDN with theme tokens for Obsidian (`#004B44`, `#00322d`) and Manrope.  
- Intentionally avoids the main AdminLTE shell so the page stays lightweight.

**Logged-in shell (`template_private.html`)**  
- All pages that extend `web/templates/web/template_private.html` render their markup inside `<div class="private-app-shell">`.  
- Child templates define **`{% block private_main %}`** (not `{% block content %}` directly — the parent wraps `private_main` in the shell `div`).  
- **Stylesheet:** `web/static/web/css/obsidian-private-pages.css`, linked from `{% block extra_head %}` in `template_private.html`. It styles generic AdminLTE cards, prose, tables, DataTables chrome, forms, and alerts on About, Client (minFQ), Help, Messages, remote control, flowcell manager, profile, password change, reference/primer managers, etc.  
- **Does not override** `.flowcells-page` or `.flowcell-detail-page`; those use dedicated CSS below.  
- Pages that add their own CSS in `{% block extra_head %}` should call **`{{ block.super }}`** first so `obsidian-private-pages.css` stays loaded.

**Flowcells index (`/web/private/flowcells`)**  
- Template: `web/templates/web/flowcells.html` (extends `template_private.html`).  
- Styles: `web/static/web/css/obsidian-flowcells.css`, scoped under `.flowcells-page`.  
- Manrope is loaded globally in `template_base.html`.  
- **Why not Tailwind on data pages:** AdminLTE, Bootstrap, and DataTables are already loaded; Tailwind’s preflight would fight those. Scoped CSS mirrors the Tailwind tokens from §3.  
- **DataTables:** `web/static/web/js/controllers/FlowcellTableController.js` renders badge-style cells (active/archived/runs/permission) as HTML where helpful; row click behaviour unchanged.

**Flowcell detail (`/web/private/flowcells/<id>/`)**  
- Template: `web/templates/web/flowcell_index.html`; loads `obsidian-flowcells.css` (shared buttons/breadcrumbs) and `obsidian-flowcell-detail.css` scoped under `.flowcell-detail-page` (page header, badges, Obsidian tab bar, card shell).  
- Tab-specific styling in `obsidian-flowcell-det8c300af59ea2fb4216c2435f749610317924727f ail.css` is scoped with IDs such as `#tab-summary-data`, `#tab-basecalled-data`, `#tab-notifications`, `#tab-reads`, `#tab-tasks`, `#tab-sharing` (cards, tables, charts, settings-style panels).

**Settings drawer (header bar, user menu)**  
- AdminLTE **control sidebar** in `template_private.html` (`obsidian-settings-sidebar`). The panel body uses **Tailwind** via CDN in `{% block extra_head %}` with **`corePlugins.preflight: false`** and **`prefix: 'tw-'`** on all Tailwind classes so utilities never collide with Bootstrap class names such as **`collapse`** on `#navbar-minotour-content` (Tailwind’s `collapse` utility would otherwise break the main nav). Obsidian tokens and Manrope match `template_login.html`.  
- `obsidian-chrome.css` includes a **defensive `@media (min-width: 768px)`** rule for `#navbar-minotour-content` so the expanded navbar stays visible if the CDN ever emits a conflicting `.collapse` rule. Shell-only tweaks for the drawer: border, shadow, z-index. Navbar user trigger: `obsidian-nav-user` (gear + username).

**Adding `extra_head` to new private pages**  
- Put page-scoped CSS in `{% block extra_head %}` (defined in `web/templates/web/template_base.html`).  
- If the page extends `template_private.html`, include `{{ block.super }}` unless you intentionally replace the default stylesheets.

**Global header & footer (logged-in and legacy public shell)**  
- **Templates:** `{% block header %}` in `web/templates/web/template_private.html` or `template_public.html`; footer in `web/templates/web/template_base.html`.  
- **Styles:** `web/static/web/css/obsidian-chrome.css` after `minotour.css`.  
- **Nav:** `#tn` with `obsidian-site-nav` — Obsidian gradient, white links, logo + wordmark; active link uses a subtle bottom emphasis (see §3 navigation).  
- **Footer:** `#sticky-footer` with `site-footer` `obsidian-footer` — white bar, Ether text, Obsidian links, logo.

## 7. Theme development guardrails (required)

When adding or updating a theme, always validate against flowcell detail tabs and DataTables-heavy views. Most regressions came from theme files styling the shell/nav/cards but not overriding high-specificity tab/table selectors from `obsidian-flowcell-detail.css`.

### Required compatibility tokens

Each theme should define (directly or via compatibility layer) values equivalent to:

- `surface` / `surface-alt` (panel backgrounds)
- `header` (table/card header background)
- `text` / `muted` (foreground hierarchy)
- `border`
- `pass` / `fail` row backgrounds

### Required selector coverage

Every theme must explicitly style these areas:

- `#tab-basecalled-data` (including `.bc-data-table` and pass/fail rows)
- `#tab-summary-data` (`.summary-data-table`)
- `#tab-notifications` (`.notif-data-table`)
- `#tab-reads` (`table.dataTable`)
- `#tab-tasks` (`.task-history-table`)
- `#tab-sharing` (`.sharing-data-table`)
- DataTables controls inside flowcell detail tabs (`length`, `filter`, `info`, `paginate`, input/select)

### Compatibility layer

- File: `web/static/web/css/obsidian-theme-compat.css`
- Loaded last in `template_base.html`, after all theme files.
- Purpose: enforce minimum readable contrast in tab cards/tables/controls across *all* themes, including newly added ones, while allowing theme-specific files to keep distinct aesthetics.

### QA checklist before shipping a theme

1. Flowcell detail tabs: Summary, Basecalled, Read Data, Tasks, Sharing, Notifications.
2. At least one table with striped rows + hover.
3. DataTables filter input and pagination readability.
4. Highcharts title, axis labels, legend, tooltip text contrast.
5. Registry table flowcell-name readability (`.fc-name` / `fc-col-name`).
6. Settings drawer remains clickable (z-index/pointer-events check).

If a theme fails any of the above, fix theme variables/selectors first; do not rely on ad-hoc per-page patches.
