JEP-G
* Look at balance between experiment and theoretical contribution
* Maybe look at Cognition as well to look at differences between how to frame

Goodness of fit

RT data and model

Plot of mean RT for Accept Reject per cell. No dot for an individual if they had fewer than X responses. Group average in summary plot from those participants who had values in individual plots.

Same for proportion of responses (Accept/Reject)


Work on something similar for veridical task

One that could work well - Dan Little and Mario and Chris Donkin have been on papers where have a 3 x 3 grid with carved out 2 x 2 grid for the double target cells.
gical Rules and the Classification of Integral-Dimension Stimuli

Rule based classification is quite like Veridical Choice specific rule (possible citation)

## Simulation strategy, 

Submitted a new estimation job. Contaminants come from uniform distribution between 0.35 (minimum) and 10 (maximum) from the filter steps.
One thing I thought of is response times simulated from the LBA won't be restricted to the same range, what strategy is appropriate here? Regenerate until values fall inside the range (re-simulate), post generation filter the same as in raw data (filter) or leave values outside the raw data range (leave).

**Guy's response**

Since so few responses fell outside the filtered window (i.e., >10s), a model that provides a good fit to data shouldn’t predict many responses outside that window. I tend to use your “leave” strategy. If this seems insufficient when inspecting the model predictions, we can revisit the decision.


## Notes from Guy about checking for recovery problems

Look a chains for random effects, see if they go to prior
Look at code, make sure not simulating from 98% instead of 2%

## Note from Paul OzMathPsych2022

Watching out for statistical facilitation for coactive models??

* Visual inspection of posteriors random effects, to see whether they are stationary

# Note from Rield (2008) reading

The coactive decision process here is really just a reflection of the EQW (equal weights strategy) where each attribute contributes equally to a decision. Perhaps an MAU (multiattr utility) - not equal weights, not sure about discoverability though
