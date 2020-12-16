# high-perf-measure-abundance

This Shiny App computes to give the users the accessibility to judge whether a pilot study contains good biomarkers and is worthy to be extended with larger sample size using different performance methods. 
Due to the affordability to study large number of biomarker candidates to find a good biomarker to diagnose or predict the target. However, because of the high cost and long time to collect and prepare the data in this type of study, researchers prefer to check the potential of the large set of biomarker candidates in a small pilot study. Users can compare the number of biomarker candidates with a better value than the given value to the corresponding number in the random dataset. If there are significantly more such biomarkers in the real dataset of the pilot study than in the random dataset, even if they lose their significance after the multiple testing, it might be worthwhile to extend the study. 


The users can use the defaults or choose from the options. 

performance metho: mAUC= multiclass AUC, entropy, AAC= Area above the cost curve, HUM= hypervolume under manifold, or misClassRate= misclassification rate 

simulate random data: Monte-Carlo test or permutation test

impute missing values: by median or random 

is positive: TRUE or FALSE the user should select TRUE for mAUC, AAC, and HUM or FALSE for entropy and misClassRate.


no.simulations: refers to the number of biomarker candidates that should be simulated to calculate the required p-values.

Conf.Interval: level for the confidence interval (i.e. 95%). 

positive class: in the case of AAC, the user has to specify the positive class in the data set. Otherwise, automatically the first class, alphabetically, will be chosen.

corrected p-value: it is only for the performance part. The p-value could be corrected either by controlling the FWER using Holm-Bonferroni correction or FDR using Benjamini and Hochberg method.




# Installation

R version should be R (≥ 4.0.2) 
In case Java are not installed, Users can see and install it in https://www.java.com/en/download/manual.jsp 

Use the function runGitHub() from the package shiny:

	library(shiny)
	runGitHub("high-perf-measure-abundance","Amani-Al-Mekhlafi",subdir = "R")

The following packages are downloaded or imported when starting this app:
"Hmisc"
"pROC"
"discretization"
"plotly"
"Biocomb"
"data.table"
"DT"
