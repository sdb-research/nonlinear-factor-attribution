################################################################################################################################################
#The below is research code illustrating the methodology proposed in the paper "nonlinear factor attribution" by De Boer (2019), 
#provided for evaluation and testing purposes only; no warranty on accuracy or applicability of any kind is provided.
#The (modified) Barra data set sourced from the PA library is for illustration purposes only: 
#there appear to be some obvious data issues which aren't addressed here, 
#but the template is reprentative of that commonly used by different vendors of attribution tools.
#
#The relevant functions are specified below the lines of testing code, and those definitions should be run first.
#The main function defined below is run.attribution(), as described in its definition.
#The test code first illustrates how the function might be run on a realistic data set,
#and then validates the accuracy and power of the proposed statistical tests on simulated data where the true return distribution is known.
#
#NOTE: the attribution provided is arithmetically "linked" across different subperiods (simple sums) since the tests for factor nonlinearities
#only hold on that basis. However, it is straightforward to add some linking method (frongello, for instance) across time for both 
#the linear and nonlinear factor contributions (and even on the stock level) to get an attribution that adds up to the active return
#of the properly compounded portfolio returns versus the benchmark.
#Also note: we recommend running the attribution on local currency returns, while separately reporting and assessing the currency contributions.
#The reason is to avoid collinearity issues with any country effects added.
#Note that an intercept to capture market returns is added by default to the attribution.
#The estimation includes a constraint that the return across each segmentation variable (sectors, etc.) equals zero.
#It allows weighted least squares estimation of the factor returns.
#################################################################################################################################################

#packages to be sourced: "pa" package has test data
install.packages("pa"); library(pa);
install.packages("fastDummies"); library(fastDummies)
install.packages("quadprog"); library(quadprog)
install.packages("lubridate"); library(lubridate)
install.packages("doParallel"); library(doParallel)
install.packages("doRNG");library("doRNG");
options(digits=10)
 
 
#set parameter names: see definition of run.attribution() function for more detail
ret.var = "return" #forward returns
reg.var = c("sector","value", "growth") #for attribution, variables to include
benchmark.weight = "benchmark" #BM weights
portfolio.weight = "portfolio" #portfolio weights
reg.wt.var = "reg.weight"  #regression weights
reg.con.wt.var = "benchmark" #bm-weighted return across each segmentation variable is constrained to equal zero
sec.id.var = "barrid" #security ID, for stock-level output
nl.realloc.var = c("value","growth") #choice of interaction effects for nonlinear attribution
nl.realloc.intercept = 1 #recommended baseline intercept value for nonlinear attribution
date.var = "date" #date column in the input data
parallel = FALSE #the initial 3 period use case is too small to benefit
pooled = FALSE #best default setting, most robust (see notes to function)

#########################################################################################################################################	
#first test: show how attribution function matches linear attribution for PA package when unweighted regression and a single segmentation variable.
#note that 'quarter' is the modified Barra data set (template is right, but data seems scrambled, since actual data proprietary)
#from the PA package while 'jan' is just the January data from that.
#this also illustrates how the run.attribution function can be used for standard linear attribution.
#########################################################################################################################################	
data(quarter); data(jan);
regr = run.attribution(x.all=jan, date.var, sec.id.var="barrid",ret.var, reg.var = c("sector","value","momentum"), portfolio.weight, benchmark.weight)
regr$linear.attribution.aggr #aggregated over all periods (just one here) and the segment variables
r1 <-regress(x = jan, date.var = "date", ret.var = "return", reg.var = c("sector","value","momentum"), benchmark.weight = "benchmark",   portfolio.weight = "portfolio")
summary(r1)	

#########################################################################################################################################	
#second test: for more than one segmentation variable, the method provided here applies a constraint on all estimated returns while the PA library 
#method starts omitting factors to bypass the multi-collinearity problem. Thus the attributions will differ but the residual will be the same;
#Most importantly, we illustrate the nonlinear attribution here, which flags the aggregate residual is borderline significant (p=11%) 
#and how small-cap stocks with large absolute 'size' exposure are disproportionately responsible for that (highly statistically significant).
#########################################################################################################################################	
reg.var = c("sector","country","momentum" ,  "value"  ,    "size"  ,     "growth","yield")
reg.wt.var = ""  #regression weights
reg.con.wt.var = "" #bm-weighted return to segmentation variables equals zero
sec.id.var = "barrid" #security ID, for stock-level output
nl.realloc.var = c("momentum" ,  "value"  ,    "size"  ,     "growth","yield") #choice of interaction effects for nonlinear attribution
nl.realloc.intercept = 1 #recommended baseline intercept value for nonlinear attribution

#run the nonlinear attribution
regr = run.attribution(x.all=quarter, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var, reg.con.wt.var,nl.realloc.var,nl.realloc.intercept, nl.test.type="actual", pooled=FALSE, parallel=FALSE)

#display the results: see run.attribution() function definition for more detail
regr$linear.attribution.aggr #aggregated attribution over time and for all segmentation variables
(regr$periodical.attributions$linear.attribution)[regr$periodical.attributions$linear.attribution$factors=="size",] #size contributions: truly outsize
(regr$periodical.attributions$linear.attribution)[regr$periodical.attributions$linear.attribution$factors=="stock.specific",]	#residuals: mostly offset size, already suggesting the benefit of size exposure to portfolio didn't materialize

regr$attribution.stats.aggr #test statistics for the linear attribution (residual significance) and the nonlinear attribution
regr$non.linear.attribution.aggr	#the nonlinear attribution: would reallocate most of the residual to the size factor, leaving just a small remaining stock specific contribution

#The below illustrates how the residual return reallocation can be fully audited to even the individual stock level;
#We identify top 5 contributors (stocks/month pairs) to negative residual, which contribute about half of it;
#Note how all have extremely negative standardized exposure to size, which explains the re-allocation ;
#(The data are likely scrambled since for other months these stocks' size exposures are 0: this is just for illustration purposes);
#Also note how the "minimal" residual return contribution can be non-zero even for zero active weight stocks since 
#in an unconstrained portfolio underweight size, they would have been held, and the "residual weights" are relative to that;
#Keep in mind: the attribution residual reflects constrained induced inefficiency relative to this hypothetical optimal portfolio with
#the same factor exposures as the actual optimal active portfolio!
#For the alternative residual return contrbutions decomposition per active weight times residual return, this is not the case;

q=which(rank(regr$periodical.attributions$stk.contributions$resid.stk.contributions.minimal)<=5)
sum(regr$periodical.attributions$stk.contributions$resid.stk.contributions.minimal)
sum(regr$periodical.attributions$stk.contributions$resid.stk.contributions.minimal[q])
merge(regr$periodical.attributions$stk.contributions[q,],quarter,by.x=c("sec.id","date"),by.y=c("barrid","date"),all.x=TRUE,all.y=FALSE)

#For validation: the residuals under the standard linear attribution as coded in PA package match
r1 <-regress(x = quarter, date.var = "date", ret.var = "return", reg.var = reg.var, benchmark.weight = "benchmark",   portfolio.weight = "portfolio")
summary(r1)	
	



#########################################################################################################################################	
#NOTE: the below is validation code for the statistical tests and NOT relevant to simply seeing how the nonlinear attribution might be run and used!
#We now create a large number of random test data sets to show the accuracy and the power of the statistical nonlinearity tests provided here.
#We can change the 'test type' to get different types of violations of the assumptions of the linear factor model.
#The larger the number of iterations, the better the estimates of the true p statistics.
#Note: data and fctors remain the same, only returns are randomly re-generated for every test iteration, per the test scenario assumptions.
#########################################################################################################################################	

x.all = create.data(type = "random", nof.periods = 2 , nof.stocks = 100, nof.factors = 2, nof.segments = 1)
reg.var = c("F1", "F2","S1")
nl.realloc.var = c("F1","F2")
nl.realloc.intercept = 1

#Test 1: assumptions of linear factor model met exactly, and we run pooled estimation of level of stock-spec risk: all test stats mathematically precise
#We see empirically, all tests for 95% confidence level indeed only reject the linear model ~5% of the time.
run.validation.tests(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="reg.weight", reg.con.wt.var="benchmark",nl.realloc.var,nl.realloc.intercept, nl.test.type="test1", pooled=TRUE, parallel=FALSE, nof.iterations = 100)

#reject.aggr.resid is the % of simulations the attribution residual was statistically significant at 95% confidence: about 5% in theory
#reject.linearity is the % of simulations the reallocation was statistically significant at 95% confidence: about 5% in theory
#sign.nl.interactions is the % of simulation the reallocation for each factor was statistically significant, about 5% in theory
#sign.periodical.resid is the average across simulations of the % of periods that the residual was significant, about 5% in theory
#the attribution.residual.distribution shows the similated PF, BM, active returns as well as residual from linear attribution and the
#result from the nonlinear reallocation, showing how the adjusted residual goes down in magnitude

#test 1: assumptions of linear factor model met exactly, and we run non-pooled estimation of level of stock-spec risk: 
#the tests in this case are an approximation based on period-specific estimates of level of specific risk, but it's robust:
#roughly all tests for 95% confidence level indeed only reject the linear model ~5% of the time
run.validation.tests(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="reg.weight", reg.con.wt.var="benchmark",nl.realloc.var,nl.realloc.intercept, nl.test.type="test1", pooled=FALSE, parallel=FALSE, nof.iterations = 100)


#test 2: the assumptions of the mathematical derivations in the paper are violated due to nonstationary periods
#in particular, the first period highly volatile specific risk yet very little active share so essentially irrelevant for attribution, 
#but "pooled" approach thrown off by it and interaction p-stats no longer hold 
#(we get H0 of linearity rejected incorrectly near 0% of simulations, which is conservative in this particular case); 
#the non-pooled approach is robust, p values hold approximately (H0 rejected incorrectly 5% of the time, which at 95% confidence is expected )
#therefore, we recommend pooled=FALSE as default setting: specific risk level estimated separately for each period
run.validation.tests(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="reg.weight", reg.con.wt.var="benchmark",nl.realloc.var,nl.realloc.intercept, nl.test.type="test2", pooled=TRUE, parallel=FALSE, nof.iterations = 100)

run.validation.tests(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="reg.weight", reg.con.wt.var="benchmark",nl.realloc.var,nl.realloc.intercept, nl.test.type="test2", pooled=FALSE, parallel=FALSE, nof.iterations = 1000)


#test3: we create a factor to which the multi-factor portfolio is exposed, but missing from the attribution
#we expect the residual to fail the linearity test more frequently (meaning to be significant), which it does here about 25%-30% of simulations
#we also expect the nonlinear interactions not to pick up on this, which indeed it doesn;t (or rather, about 5% of the time, at the stated p-value)

x.all = create.data(type = "random", nof.periods = 2 , nof.stocks = 100, nof.factors = 3, nof.segments = 1)
run.validation.tests(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="reg.weight", reg.con.wt.var="benchmark",nl.realloc.var,nl.realloc.intercept, nl.test.type="test3", pooled=FALSE, parallel=FALSE, nof.iterations = 1000)

				
#test4: we create a nonlinearity in pay-off from factor 1: piecewise linear	, by adding abs(F1) as a hidden factor	
#if the specific risk is sufficiently diversified (1000 stocks, 3 periods here), in 71% of simulations the aggregate residual is identified as 
#not just being noise, and in about 50% of simulations the reallocation is deemed statistically significant
#importantly, about half of the original attribution residual gets reallocated to factor 1, and a neglible part to factor 2
#the remainder stays as the new stock-specific risk
x.all = create.data(type = "random", nof.periods = 3 , nof.stocks = 1000, nof.factors = 2, nof.segments = 0);reg.var = c("F1", "F2");
run.validation.tests(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="reg.weight", reg.con.wt.var="benchmark",nl.realloc.var,nl.realloc.intercept, nl.test.type="test4", pooled=FALSE, parallel=TRUE, nof.iterations = 1000)


			
#test5: create a interaction in factor 1 with factor 2; about 90% -95% of the simulations, the estimation rejects the linear atttribution being noise,
#and also rejects the reallocation being shuffling around noise; importantly, it re-allocated the residual about 50/50 to F1 and F2
#the remaining residual after the reallocation is about 25% of what it was.

run.validation.tests(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="reg.weight", reg.con.wt.var="benchmark",nl.realloc.var,nl.realloc.intercept, nl.test.type="test5", pooled=FALSE, parallel=TRUE, nof.iterations = 1000)

#we now run the same, but with lower level of the intercept, meaning stocks w/ more moderate values of F1 and F2 will be classified 
#as factor driven positions. this reallates about 90% of the residual but slightly
#lowers the % of times that the hypothesis of not capturing any interactions gets rejected to about 85%
#an intercept of 2 would re-allocate less since only positions with very large exposures (at the tails) 
#are classified as factor driven positions, but the % of times the interaction there is deemed significant remains at about 85%, so 1 is better choice
run.validation.tests(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="reg.weight", reg.con.wt.var="benchmark",nl.realloc.var,nl.realloc.intercept=0.5, nl.test.type="test5", pooled=FALSE, parallel=FALSE, nof.iterations = 1000)

run.validation.tests(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="reg.weight", reg.con.wt.var="benchmark",nl.realloc.var,nl.realloc.intercept=2, nl.test.type="test5", pooled=FALSE, parallel=TRUE, nof.iterations = 1000)



 
####################################################################################
# function definitions: run first, even before including libraries
####################################################################################

rm(list=ls()) #start from clean slate

run.attribution <- function(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var = "", reg.con.wt.var="",nl.realloc.var=c(),nl.realloc.intercept=1, nl.test.type="actual", pooled=FALSE, parallel=FALSE)
{
	#####################################################################################################################################
	#this is the main function defined here, that runs both linear and the proposed nonlinear attribution method
	#specified returns need to be "forward", for the upcoming period, while holdings and factor exposures are beginning-of-period
	#
	#the function will loop over all dates included in the data object and run an attribution for each
	#standard linear attribution runs a restricted weighted least squares fit of forward returns on beginning-of-period factor exposures,
	#with the weighted sum of returns over all levels for each segmentation variable is restricted to be zero .
	#An intercept is added by default in this regression, which captures the market return.
	#the approach then runs the nonlinear factor attribution framework to reallocate the residual from the linear attribution
	#
	#the function returns the linear attribution (arithmetic sum across periods, geometric isn't coded) as well as the nonlinear attribution
	#plus test statistics for the linear attribution residual and the nonlinear attribution as described in the paper
	#it is recommended to run the attribution using local currency returns and then separately report currency contributions in the attribution
	#it would be desirable to geometrically link the linear attributions over time, as well as the nonlinear attributions for each time period,
	#but that isn't coded here: the test statistics and t-stats reported are only valid mathematically for arithmetic linking
	#
	#input parameter definitions:
	#'x.all': a data frame with all the necessary input data
	#'benchmark.weight': a string that has the name of the benchmark weights column in x.all (itself a numeric column)
	#'portfolio.weight': a string that has the name of the portfolio weights column in x.all (itself a numeric column)
	#'ret.var': a string that has the name of the (forward) returns column in x.all (itself a numeric column)
    #'date.var': a string that has the name of the dates column in x.all, specifying the attribution sub-period (dates need to be Date object)
    #'sec.id.var': a string that has the name of the security identifier column in x.all, all strings (for instance, BARRID)
    #'reg.var': a vector of the names of the attribution variables (factors, segmentation variables) to include (as indep vars in the regression)
    #'ret.wt.var': a string that has the name of the weight variable for the weighted least squares estimation underlying the attribution
    #if missing, set to equal-weighted regression
    #'ret.con.wt.var': a string that has the name of the variable for the weighted return on segmentation variables to be zero constraint 
    #if missing, set to equal-weighted average of segment returns to equal zero. 
    #Example: cap-weighted sum of estimated returns to each sector equals 0, rendering an interpretation as active sector returns instead.
    #'nl.realloc.var': a vector of the names of the interaction factors for the nonlinear attribution, if any are desired 
    #(if empty, only linear attribution is done)
    #'nl.realloc.intercept': the level of the intercept that will be included in the stock classification scheme underlying the nonlinear attribution 
    #'parallel': a boolean whether parallel processing should be applied (faster when many periods, but more difficult to debug if issues)
    #'nl.test.type': a string variable, possibly blank, that indicates if the return data should be distorted for testing purposes 
    #(list of test options is in single.period.attribution.date.filter function); ONLY USEFUL WHEN VALIDATING CODE, NOT FOR TRUE ATTRIBUTIONS
    #'pooled': a boolean indicating whether the stock specific vol level should be estimated on a pooled basis for all attribution periods or 
    #separately for each period (pooled=TRUE is what the mathematical derivations for test statistics are valid for, 
    #but pooled=FALSE is approximately the same and much more robust when attribution periods are very different in terms of #stocks, 
    #vol by period, length of window, etc.; FALSE is recommended)
    #'parallel': boolean indicating the attributions for each period should be done in parallel: faster when many periods, harder to debug
    #
    #the output is a list with the following components:
    #periodical.attributions: detailed attribution for each period (both at factor level and separately, individual stock level)
    #note how the "minimal" residual return contribution can be non-zero even for zero active weight stocks since 
	#in an unconstrained portfolio underweight size, they would have been held, and the "residual weights" are relative to that;
	#keep in mind: the attribution residual reflects constrained induced inefficiency relative to this hypothetical optimal portfolio with
	#the same factor exposures as the actual optimal active portfolio!
	#for the alternative residual return contrbutions decomposition per active weight times residual return, this is not the case;
    #some of the information in this object is just for the purposes of calculating test statistics, and not relevant to report
    #linear.attribution.aggr: arithmetically linked aggregate attribution over time and all segmentation variables
    #non.linear.attribution.aggr: the nonlinear residual reallocation from the prposed method, to the interaction factors passed on 
    #attribution.stats.aggr: test statistics on statistical significance of the linear attribution residual, as well as the nonlinear reallocation
    #low p values indicate statistically significant impact of nonlinear interactions; see the paper for details
    #a statistically significant attribution residual indicates some latent factor exposure of the portfolio
    #a statistically significant "interaction effect" means the residual reallocation is unlikely to just be shifting around stock specific risk
    #and one of the included interaction terms disporportionately explains the attribution residual    
    #####################################################################################################################################
	
    	
	#input data checks, though relies on user to specify valid segments, factor exposures (non-missing), returns, etc.
    stopifnot(is.data.frame(x.all))
    stopifnot(length(benchmark.weight) == 1)
    stopifnot(length(portfolio.weight) == 1)
    stopifnot(length(ret.var) == 1)
    stopifnot(length(date.var) == 1)
    stopifnot(is.Date(x.all[[date.var]]))
    stopifnot(all(sapply(c(sec.id.var, reg.var, portfolio.weight, benchmark.weight, reg.wt.var, reg.con.wt.var),  is.character)))
    stopifnot(all(c(date.var, sec.id.var, reg.var, portfolio.weight, benchmark.weight) %in% names(x.all)))
    
    if(length(nl.realloc.var)>0)
    {
    	stopifnot(all(sapply(c(nl.realloc.var),  is.character)))
    	stopifnot(all(c(nl.realloc.var) %in% names(x.all)))
    }
    
    if(reg.wt.var != "")
    {
    	stopifnot(all(c(reg.wt.var) %in% names(x.all)))
        stopifnot(is.numeric(x.all[[reg.wt.var]]))
        stopifnot((x.all[[reg.wt.var]] > 0)) #need to be strictly non-zero, since use for specific risk estimates
	}
	
    if(reg.con.wt.var != "")
    {
    	stopifnot(all(c(reg.con.wt.var) %in% names(x.all)))
    	stopifnot(is.numeric(x.all[[reg.con.wt.var]]))
    	stopifnot((x.all[[reg.con.wt.var]] >= 0))
    }
      
    stopifnot(is.numeric(nl.realloc.intercept))
    stopifnot(is.numeric(x.all[[portfolio.weight]]))
    stopifnot(is.numeric(x.all[[benchmark.weight]]))
    stopifnot(is.numeric(x.all[[ret.var]]))
	stopifnot((nl.realloc.intercept >= 0))
		
	#initialize attribution output
	attribution = list(dates = c())
	no_cores <- detectCores()   
	dates = unique(x.all[[date.var]]); dates = dates[order(dates)];
	stopifnot((length(dates) > 0))

	if(parallel & no_cores > 1)
	{
		#parallel processing should be faster when larger number of periods
		#for small number of periods, the overhead of parallel processing offsets the benefits
		start.time = Sys.time()
		cl <- makeCluster(no_cores, type="FORK")  
		registerDoParallel(cl)

		clusterEvalQ(cl, {
    	library(quadprog)
    	library(fastDummies)})
		
		attribution <- foreach(date = dates,.combine =  concatenate.attributions) %dopar% single.period.attribution.date.filter(x.all, date,date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var, reg.con.wt.var,nl.realloc.var,nl.realloc.intercept, nl.test.type) 
		
		stopCluster(cl);
		showConnections()
		end.time = Sys.time()
		end.time - start.time
	}else
	{	
		#simple loop over all dates: easier to debug if need be
		start.time = Sys.time()
		for(date in dates)
		{
			#need to write it such that output is a regression object, and we merge them
			regr = single.period.attribution.date.filter(x.all, date, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var, 		reg.con.wt.var,nl.realloc.var,nl.realloc.intercept , nl.test.type) 

			attribution = concatenate.attributions(attribution, regr )
			rm(regr);
		}
		end.time = Sys.time()
		end.time - start.time
		
	}
	
	#get arithmetic attribution; for proper geometric linking, can use reshape function to get "wide" format and then use Frongello.R library available online
	#note: for this, can geometrically link the nonlinear estimated contributions at same time, OR already aggregate them with the linear contributions of same factors
	#note: the statistical testing of nonlinear attributions only applies to arithmetic linking: see it as testing of periodical averages
	linear.attribution.aggr = tapply(as.matrix(attribution$linear.attribution[,"estimated.factor.contributions"]), data.frame(attribution$linear.attribution[,"factor.category"]),sum)
	
	#first, test aggregate.residual for whether it plausibly reflects true stock-specific noise
	#note: pooled test assumes stationary specific risk in all periods (so periods of equal length too) and estimates this from pooled residuals
	#theoretically this is when test statistics on nonlinearity are exact; when not pooled, they are a close approximation using specific vol estimated period by period
	attribution.stats.aggr = colSums(attribution$attribution.stats[,-which(colnames(attribution$attribution.stats)=="date")])
	attribution.stats.aggr["sigma.resid"] = sqrt(sum(attribution$attribution.stats$sigma.resid^2))
    
	if(pooled){attribution.stats.aggr["sigma.resid"] = sqrt(attribution.stats.aggr["v.resid"] * (attribution.stats.aggr["wss1"] /  attribution.stats.aggr["dof1"]))}
	attribution.stats.aggr["t.resid"] = attribution.stats.aggr["attribution.residual"] / attribution.stats.aggr["sigma.resid"]
	attribution.stats.aggr["p.resid"] = pf(attribution.stats.aggr["t.resid"]^2,1,attribution.stats.aggr["dof1"],lower.tail = FALSE) #p value of F test (=two sided t)
	
	#then, nonlinearity test: do the selected interaction effects identify factors disproportionately responsible for the attribution residual
	#per the logic in De Boer (2019) that classifies each stock's by the type of factor (non-directionally) per squared factor exposures
	#interaction effects could be set of main alpha factors OR some factors related to PF construction inputs (volatility, size, liquidity, etc.)
	#the former will identify if any alpha factors appeared to fall short (or exceed) linear estimated return contrubution
	#the latter would aim to characterize for what type of stocks the alpha factor performance deviated from linear estimate
	if(length(nl.realloc.var) > 0)
	{
		
	non.linear.attribution.aggr = aggregate((attribution$non.linear.attribution[,3:ncol(attribution$non.linear.attribution)]), data.frame(attribution$non.linear.attribution[,"factors"]),sum)
	colnames(non.linear.attribution.aggr)=c("factors",colnames(attribution$non.linear.attribution[-c(1:2)]))
	
	cov.delta = non.linear.attribution.aggr[,paste0("delta.vs.",non.linear.attribution.aggr[,"factors"])]
	rownames(cov.delta) = non.linear.attribution.aggr[,"factors"]

	if(pooled)
	{
		cov.delta = non.linear.attribution.aggr[,paste0("delta.v.",non.linear.attribution.aggr[,"factors"])]
		rownames(cov.delta) = non.linear.attribution.aggr[,"factors"]
		ss.pooled = (attribution.stats.aggr["wss.inter"] /  attribution.stats.aggr["dof.inter"])
		cov.delta = cov.delta * as.numeric(ss.pooled)
	}
	
	#set t-values of deviations from uninformative nonlinear attribution
	non.linear.attribution.aggr$delta.t = non.linear.attribution.aggr$delta / sqrt(diag(as.matrix(cov.delta))) 

	#set F-test for deviation from uninformative; since cov.delta degenerate, can leave our *any* of the factors here, which one is irrelevant
	omit.k = 1
	cov.delta.test = as.matrix(cov.delta[-omit.k,-omit.k])
	f.test = t(non.linear.attribution.aggr$delta[-omit.k]) %*% solve(cov.delta.test) %*% non.linear.attribution.aggr$delta[-omit.k] / length(nl.realloc.var)
	
	attribution.stats.aggr["p.nl.attrib"] = pf(f.test,length(nl.realloc.var),attribution.stats.aggr["dof.inter"],lower.tail = FALSE) 

	non.linear.attribution.aggr=non.linear.attribution.aggr[,c("factors", "contribution", "uninformative.reallocation.wts","uninformative.reallocation","delta", "delta.t")]
	non.linear.attribution.aggr[,c("uninformative.reallocation.wts")] = non.linear.attribution.aggr[,c("uninformative.reallocation.wts")] / length(dates)
	}else
	{
		#no interaction effects selected
		attribution.stats.aggr[c("p.nl.attrib")] = NA
		non.linear.attribution.aggr = data.frame()
	}
	
	#create the test statistic
	attribution.stats.aggr = attribution.stats.aggr[c("pf.return","bm.return", "pf.active.return", "attribution.residual","sigma.resid","p.resid","p.nl.attrib")]

	#create output of attribution function
	attribution.aggregate = list(periodical.attributions = attribution, linear.attribution.aggr= linear.attribution.aggr, non.linear.attribution.aggr = non.linear.attribution.aggr, attribution.stats.aggr=attribution.stats.aggr)
	
	#return attribution output
	return(attribution.aggregate)
}#end of run.attribution function


 
single.period.attribution.date.filter <- function (x.all, date, date.var, sec.id.var, ret.var, reg.var , portfolio.weight, benchmark.weight, reg.wt.var="", reg.con.wt.var="", nl.realloc.var = c(), nl.realloc.intercept = 1, nl.test.type = "actual") 
{
	##################################################################################################################
	#this function runs factor attribution for a single period as indicated by the 'date' variable
	#it is intended to be used as an auxiliary function only by the run.attribution() function, not directly
	#the output is an attribution object of information that will be concatenated over time  
    #for background on the input data, see the description of the run.attribution() function
    #the "test type" variable allows for different randomization tests to validate the generated test statistics
	##################################################################################################################

	x = x.all[x.all[[date.var]] == date,]

	if(nl.test.type != "actual")
	{	
		#for test1, factor model holds exactly, with same level of stock-specific risk (stationary) each period
		#randomly generate specific returns and factor returns
		#for test2, introduce non-stationarity (different spec risk by period, or put otherwise, different period lengths)
		#fr aggregate tests, pooled versus nonpooled only much matters when #stocks varies by period then as well
		
		periodical.stock.spec.var = 10; 
		if(nl.test.type == "test2")
		{
			if(date==as.Date("1999-12-31"))
			{
				#in first period of test, very high variance but little active position, so should be ignored...
				#yet "pooled" would use it.
				periodical.stock.spec.var = 50;
				x[,"portfolio"] = x[,"benchmark"] + rnorm(nrow(x))/10000000; 
				x[,"portfolio"] = x[,"portfolio"] /sum(x[,"portfolio"] );
			}
		}
		
		regular.vars = reg.var[sapply(data.frame(x[,reg.var]), is.numeric)]
		vars = matrix(1, nrow(x),1); if(reg.wt.var != ""){vars = 1 / x[,reg.wt.var]}
		vars = vars / mean(vars); vars = vars * periodical.stock.spec.var
		x[,ret.var] = (rnorm(1)*5+7/12) + as.matrix(x[,regular.vars]) %*% as.matrix(rnorm(length(regular.vars))+1) + rnorm(nrow(x)) * sqrt(vars)
		
		if(nl.test.type == "test3")
		{
			#create a missing variable in the attribution to which the multi-factor portfolio is exposed
			if("F3" %in% colnames(x) & !("F3" %in% reg.var))
			{
				x[,ret.var] = x[,ret.var] + x[,"F3"] * rnorm(1) 
			}
		}
		
		if(nl.test.type == "test4")
		{
			#create a nonlinearity in factor 1
			if("F1" %in% colnames(x) & ("F1" %in% reg.var))
			{
				x[,ret.var] = x[,ret.var] + abs(x[,"F1"]) * abs(rnorm(1)) 
			}
		}
		
				
		if(nl.test.type == "test5")
		{
			#create a nonlinearity in factor 1
			if("F1" %in% colnames(x) & ("F1" %in% reg.var) & "F2" %in% colnames(x) & ("F2" %in% reg.var))
			{
				x[,ret.var] = x[,ret.var] + (x[,"F1"]*x[,"F2"]) * (1 +rnorm(1))
			}
		}
		
		

	}
	
	#call the single period attribution function
	return(single.period.attribution(as.data.frame(x), date.var, sec.id.var, ret.var, reg.var , portfolio.weight, benchmark.weight, reg.wt.var, reg.con.wt.var, nl.realloc.var , nl.realloc.intercept ));
} #end of single.period.attribution.date.filter function

single.period.attribution <- function (x, date.var, sec.id.var, ret.var, reg.var , portfolio.weight, benchmark.weight, reg.wt.var="", reg.con.wt.var="", nl.realloc.var = c(), nl.realloc.intercept = 1) 
{
	##################################################################################################################
	#this function runs factor attribution for a single period
	#it is intended to be used as an auxiliary function only by the run.attribution() function, not directly
	#the output is an attribution object of information that will be concatenated over time  
    #for background on the input data, see the description of the run.attribution() function
	##################################################################################################################

	#rescale weights so that for sure they sum to 1 for the period
	#do so for both the regression weights and the weights for the regression segment constraints
	if(reg.wt.var== ""){reg.wt.var= "reg.weight"; x$reg.weight = 1;}
	if(reg.con.wt.var== ""){reg.con.wt.var= "reg.con.weight"; x$reg.con.weight = 1;}
	x[,reg.wt.var] = x[,reg.wt.var] / sum(x[,reg.con.wt.var])
	x[,reg.con.wt.var] = x[,reg.con.wt.var] / sum(x[,reg.con.wt.var])

	#create the matrix of independent variables for the linear regression, including dummies for the non-numeric classification variables
	segmentation.vars = reg.var[sapply(data.frame(x[,reg.var]), is.factor)]
	regular.vars = reg.var[sapply(data.frame(x[,reg.var]), is.numeric)]
	if(length(segmentation.vars)>0)
	{
		lm.mat <- dummy_cols(x[,reg.var],select_columns=segmentation.vars)
		lm.mat = lm.mat[,-which(sapply(data.frame(lm.mat), is.factor))]
	}else{lm.mat = as.data.frame(x[,reg.var])}	
	
	#add a market intercept and create weighted version of the matrix (multiply each column by weight vector, both regr wts and regr constraint wts)
	lm.mat = as.matrix(cbind(1,lm.mat)); colnames(lm.mat)[1]="market"
	lm.mat.w = lm.mat * matrix(x[,reg.wt.var],nrow(lm.mat),ncol(lm.mat),byrow=FALSE)
	lm.mat.cw = lm.mat * matrix(x[,reg.con.wt.var],nrow(lm.mat),ncol(lm.mat),byrow=FALSE)

	#Set up the restricted weighted least squares estimation using QuadProg module
	#This routine implements the dual method of Goldfarb and Idnani (1982, 1983) for solving quadratic programming problems 
	#of the form #min(−dT b + 1/2bT Db) with the constraints AT b >= b0. 

	epsilon = 0.0000001 #for the optimization to have unique solution, apply minor RIDGE penalty term
	Dmat = t(lm.mat.w) %*% lm.mat + epsilon*diag(ncol(lm.mat.w))
	dvec = t(lm.mat.w) %*% x[,ret.var] 
	meq = length(segmentation.vars)
	Amat = matrix(0,meq,ncol(lm.mat)); colnames(Amat)=colnames(lm.mat);
	bvec = matrix(0,meq,1)

	
	#set the constraint that for each segmentation variable, weighted (by regr constraint weights) sum of estimated coefficients equals 0
	#so for country dummies, sum of weighted country returns equals zero, giving estimated coefficients interpretation of active returns
	#see Menchero (2010), "Characteristics of Factor Portfolios" for alternatives, e.g., cap-weighting the constraint on segment returns
	for(k in min(1,meq):meq)
	{
		dummy.vars.k = which(colnames(lm.mat) %in% paste0(segmentation.vars[k],"_",unique(x[,segmentation.vars[k]])))
		Amat[k,dummy.vars.k] = colSums(lm.mat.cw)[dummy.vars.k]
	}
	

	#now run the restricted regression that estimates factor returns and evaluate the fit
	#note that weighted segment returns are constrained to sum to zero
	coef = data.frame(estimated.factor.returns = solve.QP(Dmat, dvec, t(Amat), bvec, meq)$solution)
	rownames(coef) = colnames(lm.mat.w)
	fit = lm.mat %*% as.matrix(coef);  #estimated systematic part of stock returns
	resid = as.matrix(x[,ret.var]  - fit);         #estimated specific stock returns
	colnames(resid) = "residual.returns"
	dof = nrow(lm.mat.w) - ncol(lm.mat.w) + meq;   #degrees of freedom
	wss=sum(resid^2 * x[,reg.wt.var])              #weighted sum of squared residuals
		
		
	#now get portfolio active weight and corresponding factor exposures
	#from that and the estimated factor returns, get estimated factor contributions as the product between those,
	#and also calculate the attribution residual (return left unexplained by attribution)
	active.weights = as.matrix(x[[portfolio.weight]] - x[[benchmark.weight]])
	active.factor.exposures = data.frame(active.factor.exposures = t(lm.mat)%*%active.weights)
	estimated.factor.contributions = data.frame(estimated.factor.contributions = active.factor.exposures*coef)
	colnames(estimated.factor.contributions)="estimated.factor.contributions"
	pf.return = sum(x[[portfolio.weight]] * x[[ret.var]])
	bm.return = sum(x[[benchmark.weight]] * x[[ret.var]])
    pf.active.return = sum(active.weights * x[[ret.var]])
	attribution.residual = pf.active.return - sum(estimated.factor.contributions ,na.rm=TRUE)
	
	#now store the final attribution and add the unexplained residual
	#create a row for each part of the attribution, which decomposes portfolio active return
	linear.attribution = data.frame(date = (x[1,date.var]),factors = c(paste0(colnames(lm.mat.w)),"stock.specific"))
	linear.attribution$factor.category = paste0(linear.attribution$factors)
	linear.attribution$active.factor.exposures = as.matrix(rbind(active.factor.exposures,1)	)	
	linear.attribution$estimated.factor.returns = as.matrix(rbind(coef,attribution.residual ))
	linear.attribution$estimated.factor.contributions = as.matrix(rbind(estimated.factor.contributions,attribution.residual ))
	for(k in min(1,meq):meq)
	{
		#add categorization of the factors for later grouping of output
		dummy.vars.k = which(linear.attribution$factors %in% paste0(segmentation.vars[k],"_",unique(x[,segmentation.vars[k]])))
		linear.attribution[dummy.vars.k,"factor.category"] = paste0(segmentation.vars[k])
	}
	
	#remove rownames of the columns since duplication thereof will get in the way of combining attributions over time
	rownames(linear.attribution[,4]) = c(); rownames(linear.attribution[,5]) = c();	rownames(linear.attribution[,6]) = c();
	
   #also store output of regression; difficult to get t-stats for individual factor returns due to the regression constraint
   attribution.stats = data.frame(date = (x[1,date.var]),pf.return = pf.return, bm.return = bm.return, pf.active.return = pf.active.return, attribution.residual = attribution.residual)
   attribution.stats$dof = dof #degrees of freedom
   attribution.stats$wss = wss #weighted sum of squares
   attribution.stats$sigma = sqrt(wss/(dof))
   
	#free up some memory
	rm(active.factor.exposures); rm(Amat); 
	rm(bvec); rm(coef); rm(Dmat); rm(dvec); rm(estimated.factor.contributions);
	rm(lm.mat); rm(lm.mat.w); rm(lm.mat.cw);  gc();

	#now get the vector of stock weights that generates the attribution residual that allows calculcating its t-stat, even for single-period attribution
	#if the residual is statistically significant, this is evidence of miss-specification of the linear model in a way that affected the attribution
	#if the residual is "large" but not significant, this tends to suggest the portfolio is highly concentrated and poorly diversified
	#see De Boer (2019): the return to this "residual portfolio" is the residual, allowing us to estimate its distribution
	#we calculate this residual generating weights as the residual (transformed) from regressing portfolio active weights (transformed) 
	#on the regression variables for the attribution; we use recursion here and call the function itself, albeit with the ret.var labeled as ""transformed.active.weight""
	#first, initialize the outputs of this regression to zero, in case not run (when we already are inside the recursion call)
	stk.contributions = data.frame(date = (x[1,date.var]),sec.id = x[[sec.id.var]], residual.returns = resid); 
	if(ret.var != "transformed.active.weight" & !("w.resid.scaled" %in% reg.var) & !("w.resid.scaled1" %in% reg.var))
	{
		#add transformed portfolio active weights to the regression
		x1 = cbind(x,active.weights / x[[reg.wt.var]])
		colnames(x1)[ncol(x1)] = "transformed.active.weight"
		regr1 = single.period.attribution(x1, date.var, sec.id.var, ret.var="transformed.active.weight", reg.var , portfolio.weight, benchmark.weight,reg.wt.var,reg.con.wt.var ) 

		#sanity check: return should equal the attribution residual
		#mathematically, getting the residual is same as putting scaled dependend variable in constraint regression, and scaling back the resultant residuals
		w.resid = as.matrix(regr1$stk.contributions$residual.returns * x[[reg.wt.var]] )
		c(sum(w.resid * x[[ret.var]]), attribution.residual )
		
		#now, again run this attribution but with residual weights added as regressor so can estimate its volatility
		#recall from De Boer (2019) that estimated contribution from this added variable will subsume the original residual
		x1 =  cbind(x,w.resid / x[[reg.wt.var]] / (sum(w.resid^2 / x[[reg.wt.var]])))
		colnames(x1)[ncol(x1)] = "w.resid.scaled"
		regr1 = single.period.attribution(x1, date.var, sec.id.var,ret.var, c(reg.var, "w.resid.scaled"), portfolio.weight, benchmark.weight,reg.wt.var,reg.con.wt.var ) 

		#sanity check 2: residual from this attribution ("stock specific") should be 0 and contribution of added variable ("w.resid.scaled") equals the original residual
		#exposure to the added factor equals 1 due to proper scaling
		regr1$linear.attribution[regr1$linear.attribution$factors%in%c( "w.resid.scaled","stock.specific"),]

		#store a minimal decomposition of stocks' contributions to the residual
		#this has lower ex-ante "magnitude" than the standard decomposition of residual return contributions 
		#the former looks at PF active positions in excess of the factor mimicking portfolios, the latter looks at positions in excess of BM
		stk.contributions$act.stk.contributions = matrix(active.weights * x[[ret.var]])
		stk.contributions$resid.stk.contributions.minimal = matrix(w.resid * resid)
		stk.contributions$resid.stk.contributions.standard = matrix(active.weights * resid)
		stk.contributions$w.resid = matrix(w.resid)
		
		#now store the stats from this regression, allows calculating t-stat of residual this period and in aggregate over time	
		#recall that we assume regression weights are inversely proportional to stock specific variances of all stocks
		attribution.stats$dof1 = regr1$attribution.stats$dof #degrees of freedom of regression with added variable
  		attribution.stats$wss1 = regr1$attribution.stats$wss  #weighted sum of squares
  		attribution.stats$v.resid = (sum(w.resid^2 / x[[reg.wt.var]]))  #variance multiplier (with specific risk level) of residual return
		attribution.stats$sigma1 = regr1$attribution.stats$sigma #estimated stock specific vol level, statistically indep of return to residual portfolio
		attribution.stats$sigma.resid = sqrt(attribution.stats$v.resid) * attribution.stats$sigma1 #estimated volatility of attribution residual
				
		#lastly, calculate the t-stat of the residual return in the attribution, even for a single period, 
		#based on assumptions of attribution regression model
		attribution.stats$t.resid = attribution.residual / attribution.stats$sigma.resid

		#free up some memory
		rm(x1); rm(regr1); gc(); rm(w.resid);
	}
	
	
	#now start including the nonlinear reallocation if so desired (otherwise just return the linear residual)
	#this follows the nonlinear residual reallocation rule proposed in De Boer (2019) that reallocates stocks'
	#contributions to the attribution residual in proportion to their squared exposures to each of the "interaction factors"
	#that have been passed on by the user, with an "intercept" siphoning off the residual contributions of stocks with moderate to low
	#such factor exposures. if residual return contributions get allocated disproportionately to certain factors, this
	#is evidence of nonlinear interactions of the factor returns (with portfolio construction constraints), 
	#which is heuristically reflected in the nonlinear attribution
	#that reallocated the original attribution residual from the linear regression; a better approach, but very time-consuming and manual,
	#would be to identify what interactions caused problems and include those in the attribution model (see De Boer and Jeet [2016] for examples)
	non.linear.attribution = data.frame(date = (x[1,date.var]),factors = c("nl_stock.specific"), contribution = attribution.residual)
	
	if(length(nl.realloc.var) > 0 & "w.resid" %in% colnames(stk.contributions) )
	{
		#create the matrix for the NL residual reallocation
		#first standardize the exposures and add the "intercept" (use simple z-scoring, but weighted should be fine too)
		z =  as.matrix(x[,nl.realloc.var]);
		z = z - matrix(colMeans(z),nrow(z),ncol(z),byrow=TRUE)
		z = z / matrix(apply(z,2,sd),nrow(z),ncol(z),byrow=TRUE)
		z = cbind(z, nl.realloc.intercept ) 
		colnames(z) = c(paste0("nl_",nl.realloc.var),"nl_stock.specific")

		#then,create the weight matrix for reallocating the residual & premultiply this by residual portfolio weights
		#in effect, allocating residual portfolio weights across the reallocation factors
		w.realloc = z^2 / matrix(rowSums(z^2),nrow(x),ncol(z),byrow=FALSE)
		x.realloc = matrix(stk.contributions$w.resid,nrow(stk.contributions),ncol(w.realloc),byrow=FALSE) * w.realloc
		
		#sanity check: sum of return on reallocation matrix should equal the attribution residual
		t(x.realloc) %*% x[[ret.var]]; sum(t(x.realloc) %*% x[[ret.var]]);		
		
		#then to make "minimal", back out "factor exposures" by taking residual from regression x.realloc on the attribution factors
		#this gives the projection matrix P.realloc that created the residual reallocation per De Boer (2019)
		#mathematically, getting the residual is same as putting scaled dependend variable in constraint regression, and scaling back the resultant residuals
		p.realloc = x.realloc; regr1 = c();
		for(i in 1:ncol(x.realloc))
		{
			x.realloc1 = cbind(x,x.realloc[,i] / x[[reg.wt.var]])
			colnames(x.realloc1)[ncol(x.realloc1)] = "transformed.active.weight"
			regr1 = single.period.attribution(x.realloc1, date.var, sec.id.var, ret.var="transformed.active.weight", reg.var , portfolio.weight, benchmark.weight,reg.wt.var, reg.con.wt.var) 
			p.realloc[,i] = as.matrix(regr1$stk.contributions$residual.returns * x[[reg.wt.var]] )			
		}
		
		#sanity check: return on adjusted reallocation matrix should equal the attribution residual
		#also, return on reallocation matrix p should equal the attribution stock level residuals reallocated proportionately to squared factor exposures
		t(p.realloc) %*% x[[ret.var]]; sum(t(p.realloc) %*% x[[ret.var]]); t(w.realloc) %*% stk.contributions$resid.stk.contributions.minimal 
		
		#as is noted in De Boer (2019), this approach will always reallocate the residual in full, which may not be justified
		#to test whether there was some pattern in the residual stock contributions where they predominatly arose in a particular segment of factor bets 
		#(as identified by large squared factor exposures)
		#we should compare the nonlinear reallocation to what would be the default reallocation when there is no such pattern, meaning the residual return contributions
		#simply get redistributed per the average reallocation weight; the matrix p0.realloc captures that non-informative reallocation
		p.realloc0 = matrix(stk.contributions$w.resid,nrow(stk.contributions),ncol(w.realloc),byrow=FALSE) * matrix(colMeans(w.realloc),nrow(p.realloc),ncol(p.realloc),byrow=TRUE)

		#sanity check: this too should reallocate the residual in full
		t(p.realloc0) %*% x[[ret.var]]; sum(t(p.realloc0) %*% x[[ret.var]])
		
		#to run a statistical test to see if the nonlinear reallocation from p.realloc is different from the simple non-informative reallocation p.realloc0,
		#we need an estimate of the stock specific noise level that is statistically indep from this difference
		#for this, we run a regression with all reallocation "portfolios" (columns of p.realloc) added
		x.inter = p.realloc / matrix(x[[reg.wt.var]],nrow(x),ncol(p.realloc),byrow=FALSE) 
		x.inter = x.inter / as.numeric((sum(stk.contributions$w.resid^2 / x[[reg.wt.var]]))) #scaling to make exposures to added variables sum to 1
		x.inter = cbind(x,x.inter)
		colnames(x.inter)[-c(1:ncol(x))] = c("w.resid.scaled",paste0("w.resid.scaled",c(1:(ncol(p.realloc)-1))))
		regr.inter = single.period.attribution(x.inter, date.var, sec.id.var,ret.var, c(reg.var, colnames(x.inter)[-c(1:ncol(x))]), portfolio.weight, benchmark.weight,reg.wt.var,reg.con.wt.var ) 

		#sanity check: attribution resulting from this will give zero remaining residual
		#however, while sum of exposures to added variables equals 1, the estimated contributions do not match the residual reallocation rule
		#due to multi-collinearity of the reallocation variables: while residual is zero, this is less robust approach than the residual reallocation heuristic
		#therefore, we just use this regression to allow statistical testing
		regr.inter$linear.attribution[regr.inter$linear.attribution$factors%in%c( colnames(x.inter)[-c(1:ncol(x))],"stock.specific"),]
		t(w.realloc) %*% stk.contributions$resid.stk.contributions.minimal 
		
		#now store the regression output
		attribution.stats$dof.inter = regr.inter$attribution.stats$dof #degrees of freedom of regression with added variable
  		attribution.stats$wss.inter = regr.inter$attribution.stats$wss  #weighted sum of squares
  		attribution.stats$sigma.inter = regr.inter$attribution.stats$sigma #estimated stock specific vol level, statistically indep of return to residual portfolio
		
		#now create the nonlinear residual reallocation, the baseline non-informative version thereof, the difference between them, and the
		#covariance matrix of the difference between them, so that can test if the differences are statistically significant
		#if so, this indicates the residual return contributions aren't white noise, but rather have some nonlinear relationship to one or more factor bets
		
		non.linear.attribution = data.frame(date = (x[1,date.var]),factors = colnames(p.realloc), contribution = t(p.realloc) %*% x[[ret.var]],
			uninformative.reallocation.wts = colMeans(w.realloc),uninformative.reallocation = t(p.realloc0) %*% x[[ret.var]] )
		non.linear.attribution$delta = non.linear.attribution$contribution - non.linear.attribution$uninformative.reallocation 	
	
	 	#get the differences in "portfolio weights" between the respective informative and uninformative reallocations and
	 	#then get specific risk thereof (inversely proportional to regression weights)
		p.realloc.delta = p.realloc - p.realloc0
		p.realloc.delta.cov = t(p.realloc.delta *  matrix(1/x[[reg.wt.var]],nrow(x),ncol(p.realloc),byrow=FALSE) ) %*% p.realloc.delta
		non.linear.attribution$delta.t = non.linear.attribution$delta / sqrt(diag(p.realloc.delta.cov ) ) / attribution.stats$sigma.inter 
		non.linear.attribution[,paste0("delta.v.",colnames(p.realloc))] = p.realloc.delta.cov #(degenerate!) variance of the differences, so that can aggregate over time
		non.linear.attribution[,paste0("delta.vs.",colnames(p.realloc))] = p.realloc.delta.cov * attribution.stats$sigma.inter^2 #same, already multiplied by SS variance
		rownames(non.linear.attribution)=c()
		#free up some memory
		rm(regr1); rm(x.inter); rm(regr.inter); rm(p.realloc); rm(p.realloc0); rm(p.realloc.delta); rm(p.realloc.delta.cov); rm(w.realloc); rm(x.realloc); rm(z); gc();
	}
	
	#final attribution output object 
	regr.output = list(dates = c(x[1,date.var]), linear.attribution= linear.attribution,attribution.stats = attribution.stats,stk.contributions=stk.contributions, non.linear.attribution= non.linear.attribution)

	#free up some memory
	rm(active.weights); rm(resid); rm(fit); rm(linear.attribution); rm(attribution.stats); rm(stk.contributions); rm(non.linear.attribution); gc();

	#return final object
	return(regr.output)
} #end of single.period.attribution



concatenate.attributions <- function(attribution1, attribution2)
{
	##########################################################################
	#this function combines two attribution objects from different periods
	##########################################################################
	
	if(length(attribution1$dates)==0){return(attribution2);}
	
	attribution2$dates = c(attribution1$dates, attribution2$dates)
	
	attribution2$linear.attribution = rbind(attribution1$linear.attribution, attribution2$linear.attribution)
	
	attribution2$attribution.stats = rbind(attribution1$attribution.stats, attribution2$attribution.stats)
	attribution2$stk.contributions = rbind(attribution1$stk.contributions, attribution2$stk.contributions)
	attribution2$non.linear.attribution = rbind(attribution1$non.linear.attribution, attribution2$non.linear.attribution)
	
	return(attribution2);
}


create.data <- function(type = "barra.sample", nof.periods = 3, nof.stocks = 3000, nof.factors = 5, nof.segments = 2)
{
	####################################################################################################################
	#this function creates the attribution data object
	#when type equals "barra.sample", we return the dataset from the PA library with sample Barra data
	#otherwise, it randomly creates a time series of attribution data of the desired dimensions for testing purposes
	####################################################################################################################
	
	if(type == "barra.sample")
	{
		#retrieve data from the PA library, named "quarter" which contains 3 months		
		x.all = quarter
		x.all$cap.usd[which(x.all$cap.usd > 1000000000000)] = 0 #over a trillion is data error for this historical period
		x.all$cap.usd[which(x.all$cap.usd < 10000000)] = 10000000 #floor at 10 MM
				
		#set regression weights; for the time period at hand, any MCAP of over a trillian USD is data error
		#sqrt(cap) or ln(cap) are useful; note: since we assume stock specific risk inversely proportional to regression weights,
		#they need to be non-zero

		x.all$reg.weight = sqrt(x.all$cap.usd)
		x.all$reg.weight[x.all$reg.weight < sqrt(100000000)] = sqrt(100000000)
		x.all$reg.weight[x.all$reg.weight > sqrt(500000000000)] = sqrt(500000000000)
	}else
	{
		#randomly create test data for the specified number of periods, number of stocks, and number of factors and segmentation variables
		x.all = data.frame()
		dates = seq(as.Date("1999/12/31"), by = "day", length.out = nof.periods)

		for(t in 1:nof.periods)
		{
			x.all.t = data.frame(date = dates[t],barrid = paste0("B",c(1:nof.stocks)))
			x.all.t[,paste0("F",c(1:nof.factors))] = matrix(rnorm(nof.factors * nof.stocks),nof.stocks,nof.factors)
			if(nof.segments > 0)
			{
				for(i in 1:nof.segments)
				{
					x.all.t[,paste0("S",i)] = as.factor(paste0("S",i,"_",matrix(floor(runif(nof.stocks)*(i+2))+1,nof.stocks,1)))
				}		
			}	
			
			#randomly create market caps as well as (independently, though could also have just taken sqrt(cap)) regression weights
			x.all.t[,"cap.usd"] = as.matrix(runif(nof.stocks) * 1000000000)
			x.all.t[,"reg.weight"] = 1 / (5 + runif(nof.stocks)*95)
			
			#randomly create stock returns
			x.all.t[,"return"] = as.matrix(rnorm(nof.stocks) / 100)
			
			#create a long-only "smart beta multi-factor" portfolio of sorts by taking average over factor exposures, truncated at zero, and e/w benchmark
			x.all.t[,"portfolio"] = rowMeans(x.all.t[,paste0("F",c(1:nof.factors))]) 
			x.all.t[x.all.t[,"portfolio"]<0,"portfolio"] = 0; x.all.t[,"portfolio"] = x.all.t[,"portfolio"]  / sum(x.all.t[,"portfolio"])
			x.all.t[,"benchmark"] = 1 / nof.stocks
			
			if(t == 1){x.all = x.all.t}else{x.all = rbind(x.all, x.all.t)}	
		}	
	}
		
	return(x.all)
} #end of create.data function


concatenate.aggregated.attributions <- function(attribution.aggregate1, attribution.aggregate2)
{
	##########################################################################
	#this function combines two aggregated attribution objects for purposes of
	#randomized testing only: allows validating p-stats
	##########################################################################
		
	attribution.aggregate3 = list()
	attribution.aggregate3$periodical.attributions$attribution.stats  = rbind(attribution.aggregate2$periodical.attributions$attribution.stats,attribution.aggregate1$periodical.attributions$attribution.stats)
	attribution.aggregate3$non.linear.attribution.aggr = rbind(attribution.aggregate2$non.linear.attribution.aggr ,attribution.aggregate1$non.linear.attribution.aggr)
	attribution.aggregate3$attribution.stats.aggr = rbind(attribution.aggregate2$attribution.stats.aggr,attribution.aggregate1$attribution.stats.aggr)

	return(attribution.aggregate3);
}#end of concatenate.aggregated.attributions function


run.attribution(x.all=quarter, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var, reg.con.wt.var,nl.realloc.var,nl.realloc.intercept, nl.test.type="actual", pooled=FALSE, parallel=FALSE)



run.validation.tests<-function(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var="", reg.con.wt.var="",nl.realloc.var=c(),nl.realloc.intercept=1, nl.test.type="test1", pooled=FALSE, parallel=FALSE, nof.iterations = 100)
{
	########################################################################################################################################################################
	#this function runs random validation tests for the nonlinear attribution methodology
	#the parameter 'nof.iterations' indicates how many random tests to run to evaluate the p-stats and power of the proposed statistical tests 
	#note that the data all remains the same, except the returns are randomly generated per the test scenario
	#if 'parallel' is true the simulations are run in parallel processing mode, which is quicker but harder to debug when something goes wrong
	#the parameter 'nl.test.type' flags what type of test to run, which are all specified in the ''single.period.attribution.date.filter' function
	#the following tests are coded:
	# 
	#test1: all the assumptions of the nonlinear attribution method and linear factor model are met exactly: a test at 95% confidence level should only trigger 5% of the time
	#test2: first period highly volatile specific risk yet very little active share so essentially irrelevant for attribution, but "pooled" approach thrown off by it
	#test3: portfolio exposed to a factor omitted from the attribution
	#test4: nonlinear payoff to factor 1
	#test5: interaction effect between factors 1 and 2	######################################################################################################################################################################
	
	#initialize attribution output
	random.attributions = list()
	no_cores <- detectCores()   
	
	if( no_cores > 1 & parallel)
	{
		#parallel processing should be faster when larger number of periods
		#for small number of periods, the overhead of parallel processing offsets the benefits
		start.time = Sys.time()
		cl <- makeCluster(no_cores, type="FORK")  
		registerDoParallel(cl)
		
		clusterEvalQ(cl, {
    	library(quadprog)
    	library(fastDummies)})
		
		random.attributions <- foreach(i = c(1:nof.iterations),.combine =  concatenate.aggregated.attributions) %dorng% {run.attribution(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var, reg.con.wt.var,nl.realloc.var,nl.realloc.intercept, nl.test.type, pooled, parallel=FALSE);}		
		stopCluster(cl);
		showConnections()
		end.time = Sys.time()
		end.time - start.time
	}else
	{	
		#simple loop over all dates: easier to debug if need be
		start.time = Sys.time()
		for(i in c(1:nof.iterations))
		{
			#need to write it such that output is a regression object, and we merge them
			attr = run.attribution(x.all, date.var, sec.id.var,ret.var, reg.var , portfolio.weight, benchmark.weight,reg.wt.var, reg.con.wt.var,nl.realloc.var,nl.realloc.intercept, nl.test.type, pooled, parallel)
			print(i);flush.console();
			if(i==1){random.attributions = attr}else{random.attributions = concatenate.aggregated.attributions(random.attributions, attr )}
			rm(attr);
		}
		end.time = Sys.time()
		end.time - start.time
	}

	attribution.residual.distribution=c(0,0,0,0);
	if(length(nl.realloc.var)>0)
	{
		attribution.residual.distribution = rbind(data.frame(factors = "lin_ss",contribution =random.attributions$attribution.stats.aggr[,"attribution.residual"]), random.attributions$non.linear.attribution.aggr[,c("factors","contribution")])
		attribution.residual.distribution = rbind(attribution.residual.distribution,data.frame(factors = "pf.ret",contribution =random.attributions$attribution.stats.aggr[,"pf.return"]))
		attribution.residual.distribution = rbind(attribution.residual.distribution,data.frame(factors = "bm.ret",contribution =random.attributions$attribution.stats.aggr[,"bm.return"]))
		attribution.residual.distribution = rbind(attribution.residual.distribution,data.frame(factors = "act.ret",contribution =random.attributions$attribution.stats.aggr[,"pf.active.return"]))
		attribution.residual.distribution = cbind(tapply(attribution.residual.distribution$contribution,attribution.residual.distribution$factors,mean),
		tapply(attribution.residual.distribution$contribution,attribution.residual.distribution$factors,sd))
  		colnames(attribution.residual.distribution)=c("mean","sd")
	}
	
	random.test.stats = list(
	reject.aggr.resid = mean((random.attributions$attribution.stats.aggr[,"p.resid"] <= 0.05)),
	reject.linearity = mean((random.attributions$attribution.stats.aggr[,"p.nl.attrib"] <= 0.05)),
	sign.nl.interactions = tapply(abs(random.attributions$non.linear.attribution.aggr$delta.t)>1.96,random.attributions$non.linear.attribution.aggr$factors,mean),
	sign.periodical.resid = mean(abs(random.attributions$periodical.attributions$attribution.stats$t.resid)>1.96),attribution.residual.distribution=attribution.residual.distribution, 	run.time = end.time - start.time)
	
	return(random.test.stats)
}#end of run.validation.tests function





