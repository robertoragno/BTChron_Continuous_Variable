# dataset_2 — Clay tobacco pipe stem bore dating

Binford's (1962) pipe stem formula is a regression of date on a measured value:

```
date = 1931.85 - 38.26 * (mean stem bore in 64ths of an inch)
```

This is the time-as-response direction of the paper, and it has been standard
practice in historical archaeology for sixty years. Its calibration units are
excavated contexts dated to a *range* by documents, stratigraphy or diagnostic
artifacts, which the formula collapses to a midpoint: Case 2, independent windows of
unequal width, one per context.

## Source

McMillan, Lauren K. 2010. *Put This in Your Pipe and Smoke It: An Evaluation of
Tobacco Pipe Stem Dating Methods*. MA thesis, East Carolina University.
<https://thescholarship.ecu.edu/handle/10342/2903>

The thesis reports, for 64 contexts across 26 sites in Maryland, Virginia, North
Carolina and South Carolina: counts of stems at each bore size from 9/64" to 4/64"
(Tables 4.2–4.30), the context's date, and the dates returned by the Binford, Hanson
and Heighton & Deagan formulas (Appendix A). Its inclusion criteria keep the dates
independent of the pipes — sites "had to be dated through the historical record or
with other temporally diagnostic artifacts than pipes, such as ceramics" (p. 26).

64 contexts, 17,931 measured stems, 25 to 4,751 per context. Dates span 1620–1820
with windows of 0 to 60 years, median 20; mean bore runs from 8.53 to 4.50.

## Data

| File | Contents |
|---|---|
| `data/pipe_stems_mcmillan2010.csv` | one row per context: counts by bore size, date bounds, mean bore, the three published formula dates |
| `scripts/00_transcribe.py` | rebuilds the CSV from the thesis |

The thesis PDF is not redistributed here. To rebuild the CSV, download it from the
link above, run `pdftotext -layout` on it to `data/mcmillan2010_pdf_text.txt`, then
run the script.

`date_type` records the form in which the thesis prints each date, not the method
behind it: 47 contexts have an explicit year range, 14 a decade or span of decades,
and 3 a single year with no window. St John's Circular Pit appears twice, as
"Phase I" and as "total" over a longer span, so the two rows are nested. 14 of the 64
contexts come from Clifts Plantation.

Columns in the PDF tables are positional, so each count is assigned to the bore size
whose header centre is nearest. Recomputing Binford's formula from the transcribed
counts reproduces the thesis's published Binford date to within 0.75 years for 61 of
the 64 contexts, which confirms the column assignment. The three exceptions are
inconsistencies within the thesis rather than transcription errors, and are flagged
`binford_check = DIFFERS`.

## First look

A midpoint OLS of date on mean bore over all 64 contexts gives
`date = 1864.8 - 25.41 * bore`, against Binford's `1931.85 - 38.26 * bore`. The
shallower slope is not the midpoint collapse, and not measurement error in the
context mean, whose reliability is 0.994. It is the regression direction: with
`cor(date, bore) = -0.817`, regressing bore on date and inverting the slope returns
38.07, Binford's coefficient almost exactly. Binford did not fit a least-squares line
at all: 38.26 is "the interval between the means of the Harrington time periods"
(McMillan 2010:20), so averaging had already removed the within-period scatter and
his slope is the undiluted one. Hole (1980:287) objected to that averaging at the
time, on the grounds that the regression rests on "averages of diameters (from
Harrington), not on separate individual points", so that "the variability in the data
is hidden in the averaging process" — an objection about the spread of bore
diameters within each period, not about the width of the periods.

The comparison to run is latent dates against midpoints in both directions, with
Binford's line as the reference.

## Limitations

- The context mean bore is a sample mean over 25 to 4,751 stems, so the predictor
  carries error that varies widely across contexts. The counts per bore size are
  retained so the stems can be modelled individually instead.
- The relation is not linear across the whole span: Hanson answered Binford's single
  line with ten, and Heighton and Deagan with a curve. McMillan sampled 1600–1800
  deliberately, while noting that Noël Hume (1969:300) and Higgins (1999) hold the
  formula to work only between 1680 and 1760.
- The method applies to English ball clay pipes and does not transfer outside that
  industry.
