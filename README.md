# R estimation code for "Incorporating Time-Varying Covariates in a Simple Mixture Model for Discrete-Time Duration Data"

Pete Fader and Bruce Hardie provide a technical note https://brucehardie.com/notes/037/ on incorporating time-varying covariates for a discrete-time duration mixture model. It is an extention to the (shifted) Beta-Geometric Distribution, now with covariates. They call it **G2G+covariates model**.

However, the code on how to estimate this model is not available.  I created the MLE estimation code using R.

## How to use

The user needs to prepare the data in _person-period_ format.  For example:



A very nice tutoral on how to prepare your data can be found here: https://www.rensvandeschoot.com/tutorials/discrete-time-survival/

Specifically, the R package '''discSurv''' is a great package to prepare your dataset.

To run the model, you first need to source the background functions.

'''
source("~/R/G2G/G2G_background_functions.r")
'''

Then after ensuring the data is prepared correctly, you just need to run this command:

'''
G2G_varying_MLE(Surv(exit,event) ~ sex + immigrant + foodprices, data=Scania_PersonPeriod_Train, subject="id") 
'''

The left-hand-side of formula follows the standard R convention for survival analysis, '''Surv(period, occur)'''.  Remember that the second variable '''occur''' should be 0 unless the event of interest _occurs_ at that period, which then takes the value of 1.  The right-hand-side of the formula consists of the time-varying and also time-invarying covariates. The '''data=''' indicates the data frame.  Finally, the '''subject=''' indicates the column representing the _person_, and note that it should be a text argument.   


Please contact me (kaloklee@gmail.com) if you have any questions or comments.

Some references for the (shifted) Beta-Geometric Distribution:

Fader, Peter S. and Bruce G.S. Hardie (2007), "How to Project Customer Retention," Journal of Interactive Marketing, 31 (1), 76-90.

Lee, Ka Lok, Peter S. Fader, and Bruce G.S. Hardie (2007), “How to Project Patient Persistency,” Foresight: The International Journal of Applied Forecasting, Issue 8, 31-35.
