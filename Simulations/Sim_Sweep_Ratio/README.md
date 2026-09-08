# Sweep: dating uncertainty relative to the spread of the material

**This is a side study. Nothing in Cases 1-4 depends on it, and the whole
folder can be deleted without breaking anything.** It reads the shared models
and shared scripts; it writes only inside itself.

## The question

Cases 1 and 2 each answer "does modelling the date distribution help?" at one
setting of the dating. Put side by side they appear to say that radiocarbon
benefits more than typochronology, but that reading is an artefact of how the
two were set up:

    Case 2 typochronology : windows   8-23% of the studied period
    Case 1 radiocarbon    : envelopes 30-66% of the studied period

Case 1 is simply a harder problem. This sweep varies that one quantity and
holds everything else fixed, so the axis can be read on its own:

    ratio = typical dating uncertainty / how far apart the samples truly are

## Why the period varies and the windows do not

The intuition to avoid is that the ratio is something an analyst chooses. It is
not. Both halves are properties of a dataset: how precisely each object can be
dated, and how spread out the objects are in time. A single-site Roman study
and a continental radiocarbon compilation can use identical dating methods and
sit at opposite ends of this axis purely because of the second half.

So window widths are held at a realistic typological range in YEARS, and the
studied period is what varies. That mimics the real contrast - the same dating
practice applied to material spanning two centuries or three millennia.

## Keeping the signal constant

A longer period with the same slope means a much larger total change in y, so
attenuation would be confounded with signal-to-noise. Instead each dataset is
given a total RISE across its period, and the slope is derived as
`rise / span`. Signal strength is then the same in every cell and only the
dating ratio moves.

## The oracle column

`time_ref_range` grows with the period, and the `normal(0, 10)` prior on the
normalised slope therefore means a slightly different thing in each cell. To
show that this is not what produces the trend, every cell also fits an oracle
model on the TRUE dates. The oracle sees the same axis and the same prior but
no dating error, so it must return c = 1 in every cell. If it does, any
attenuation in the other two models is the dating ratio and not the prior.

## Running it

    Rscript Simulations/Sim_Sweep_Ratio/scripts/01_design.R
    RECOVERY_WORKERS=24 Rscript Simulations/Sim_Sweep_Ratio/scripts/02_sweep.R
    Rscript Simulations/Sim_Sweep_Ratio/scripts/03_plots.R

`SWEEP_REPS` and `RECOVERY_LIMIT` cut the work down while developing.
