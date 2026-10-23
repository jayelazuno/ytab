# Step 14F.1 UI reference research

## Local RNAcross files inspected

- `docs/codex/RNAcross/R/10_ui.R`
- `docs/codex/RNAcross/R/02_constants_themes.R`
- `docs/codex/RNAcross/R/07_visualization_core.R`
- `docs/codex/RNAcross/R/08_visualization_heatmaps.R`

## Useful RNAcross UI patterns found

- Centralized theme, palette, and plot-style configuration.
- Compact control panels separated from visualization outputs.
- User-facing plot appearance controls for size, labels, palette, and heatmap display.
- Dynamic heatmap sizing and explicit label visibility controls.
- Clear separation of generated outputs/download actions from primary interactive views.
- Collapsible advanced/technical sections rather than always-visible diagnostics.
- Responsive visual framing so plots do not stretch to unreadable browser widths.

## External references searched or cloned

No external repositories were cloned. No external runtime references were added.

RNAcross plus existing Shiny/bslib conventions provided enough implementation guidance for this display-only polish pass. The reserved reference directory `codex/ui_references/` is ignored in `.gitignore` for future research clones if explicitly needed later.

## Ideas adapted into YTAB-native code

- Shared UI component helpers for cards, metric cards, warning banners, download cards, technical details, and two-column layouts.
- Shared display-only plot customization controls for plot size, plot width, text size, label mode, label angle, grid visibility, bar orientation/value labels, point size/opacity, and heatmap label controls.
- Left-control/right-visualization layouts for plot-heavy tabs.
- Static generated image cards that preserve aspect ratio and expose filenames only in collapsible details.
- App-rendered plot frames that constrain width by default to reduce horizontal overstretching.

## Ideas rejected

- No RNAcross scientific logic, expression-data assumptions, species assumptions, or data loaders were copied.
- No draggable/custom JavaScript plot editor was adopted; it would add complexity beyond the current release need.
- No external package or reference repository was added as a runtime dependency.
- No code with unclear external license was copied.

## Runtime dependency confirmation

The YTAB app does not source, import, or read runtime data from `codex/RNAcross` or `codex/ui_references`.

