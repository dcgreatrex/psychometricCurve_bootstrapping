
# Fitting psychometric functions using a bootstrapping procedure
Author:	David Greatrex, PhD Candidate, University of Cambridge.  
Date: 01/03/2016.
Language: R.
Modifications:

# Summary:
This package fits sigmoid functions to psychoacoustic data collected during psychological
experimentation at the University of Cambridge by David Greatrex.

# Data:
The data is taken from an experiment investigating the effects of rhythmic expectation on
complex decision-making. This used the psychophysical method of constant stimuli, in which
the target stimulus was varied systematically around a detection threshold.

# Purpose:
The purpose of these analytical scripts is to fit psychometric functions to each participants
data, individually by each condtion of interest. Slope and threshold values are then extracted
from each fit and compared across condition and group using pairwise t-tests.

# Method:
The process of fitting models to each participant’s data and then testing group differences using 
null hypothesis significance testing on the estimated parameters, rather than fitting one model to 
the entire dataset, has been termed the ‘summary statistics approach’ to computational modelling and 
is a widely used procedure - see Daw (2011) for a thorough review. The approach is closely related 
to hierarchical regression modelling in that it treats each parameter estimate as a random variable 
(random effect) by drawing a participant from the population at random and then running the entire 
experiment and analysis on them. In contrast, were only one model fit to the entire dataset, the 
estimated parameters would be treated as fixed effects and as a result, all between participant 
variability would be neglected (Daw, 2011, pp. 7-8).

# Bootstrapping:
To extract the most unbiased fit for each model, a bootstrapped distribution of model-fitting 
parameters was generated for each psychometric function (as in Knoblauch and Maloney, 2012; 
James et al., 2013; Piazza et al., 2013; Gold et al., 2015). This involved resampling the data with 
replacement 10,000 times to produce 10,000 bootstrapped samples and fitting each bootstrapped sample 
with the same model used on the original data. This resulted in a bootstrapped distribution of parameter 
values B0 and B1 from which the average parameter values (mean(B0) and mean(B1)) were retained. Threshold
and slope where then defined by:
Threshold = -mean(B0)/mean(B1)
Slope = mean(B1)
where threshold represented the location of the probe tone corresponding to 50% of left responses (also 
the mean of the underlying cumulative distribution function) and slope to the steepness of the psychometric 
curve. Higher slopes represented greater discrimination of the average sequence location, whereas thresholds 
smaller or larger than zero represented a left or rightward bias in the perceived location of the average
relative to the probe tone.

# To run:
Download the content of the folder locally on to your hardrive including the folder structure. 
Open the r script 'main.R'. 
Change to path on line 18 to represent the folder location in which you saved the downloaded content.
Run 'main.R'. All plots will then be saved into the downloaded filestructure.
Note - due to bootstrapping the script will take a long time to complete. ~25 minutes. You can check progress
by looking at how many plots have been output into the plot folder. There are 24 participants in total.

# References:
*[Daw, N. (2011). Trial-by-trial data analysis using computational models. In. Delgado, M., Phelps, E., & Robbins, T., (Eds). Decision making, affect, and learning: Attention and Performace XX111. Oxford, Oxford University Press, 3-38.](http://www.cns.nyu.edu/~daw/d10.pdf)

*[Knoblauch, K. and Maloney, L. T. (2012). Modeling psychophysical data in R, volume
32. Springer Science & Business Media.](http://www.springer.com/gp/book/9781461444749)

*[James, G., Witten, D., Hastie, T., and Tibshirani, R. (2013). An introduction to
statistical learning, volume 112. Springer.](http://www-bcf.usc.edu/~gareth/ISL/ISLR%20First%20Printing.pdf)

*[[Piazza, E. A., Sweeny, T. D., Wessel, D., Silver, M. A., and Whitney, D. (2013). Humans
use summary statistics to perceive auditory sequences. Psychological science,
24(8):1389–1397].(http://pss.sagepub.com/content/24/8/1389.short)

*[[Gold, J. R., Nodal, F. R., Peters, F., King, A. J., and Bajo, V. M. (2015). Auditory
gap-in-noise detection behavior in ferrets and humans. Behavioral neuroscience,
129(4):473.](http://psycnet.apa.org/journals/bne/129/4/473/)