# nonlinear-factor-attribution
The accompanying is research code illustrating the methodology proposed in the paper "nonlinear factor attribution" by De Boer (2019), 
provided for evaluation and testing purposes only; no warranty on accuracy or applicability of any kind is provided.
The (modified) Barra data set sourced from the PA library is for illustration purposes only: 
there appear to be some obvious data issues which aren't addressed here, 
but the template is reprentative of that commonly used by different vendors of attribution tools.

The main function defined here is run.attribution(), as described in its definition.
The test code first illustrates how the function might be run on a realistic data set,
and then validates the accuracy and power of the proposed statistical tests on simulated data where the true return distribution is known.

NOTE: the attribution provided is arithmetically "linked" across different subperiods (simple sums) since the tests for factor nonlinearities
only hold on that basis. However, it is straightforward to add some linking method (frongello, for instance) across time for both 
the linear and nonlinear factor contributions (and even on the stock level) to get an attribution that adds up to the active return
of the properly compounded portfolio returns versus the benchmark.
Also note: we recommend running the attribution on local currency returns, while separately reporting and assessing the currency contributions.
The reason is to avoid collinearity issues with any country effects added.
Note that an intercept to capture market returns is added by default to the attribution.
The estimation includes a constraint that the return across each segmentation variable (sectors, etc.) equals zero.
It allows weighted least squares estimation of the factor returns.

