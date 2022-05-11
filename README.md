# high-perf-measure-abundance

This Shiny App computes to give the users the accessibility to judge whether a pilot study contains good biomarkers and is worthy to be extended with larger sample size using different performance methods. 
Due to the affordability to study large number of biomarker candidates to find a good biomarker to diagnose or predict the target. However, because of the high cost and long time to collect and prepare the data in this type of study, researchers prefer to check the potential of the large set of biomarker candidates in a small pilot study. Users can compare the number of biomarker candidates with a better value than the given value to the corresponding number in the random dataset. If there are significantly more such biomarkers in the real dataset of the pilot study than in the random dataset, even if they lose their significance after the multiple testing, it might be worthwhile to extend the study. 


The users can use the defaults or choose from the options after Loading the data set, which must be saved in *.csv format, observation in rows, biomarker candidates in columns,  and contain the factor (states of disease) in the last column. 

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

The following packages are downloaded or imported when starting this app:
"Hmisc"
"pROC"
"discretization"
"plotly"
"Biocomb"
"data.table"
"DT"

Use the function runGitHub() from the package shiny:

	library(shiny)
	runGitHub("high-perf-measure-abundance","Amani-Al-Mekhlafi",subdir = "R")

New Windows will be opened; 
•	Load the data set, which must be saved in *.csv format, observation in rows, biomarker candidates in columns, and contain the factor (states of disease) in           the last column. 

•	The users have multiple choices;

•	Performance; to show the performance of each individual biomarker candidates and the p-value without/with correction using FDR or FWER

•	HiPerMAb table and curves; to show a table or to visualize the HiPerMAb

•	HiPerMAb power; to show the required number of biomarkers with a significant p-value

Users can use the default or changes the inputs. For performance and HiPerMAb table and curves, users must specify the number of simulation according to the number of biomarker candidates. The inputs are;

•	performance method: refers to the performance measurement method, which is one of the following five; mAUC for Multiclass AUC, entropy, AAC for Area above           the cost curve, HUM for Hypervolume under manifold, or misClassRate.

•	simulate random data: refers to how simulate the random data, either by Mon-te-Carlo or permutation test.

•	impute missing values: refers to simple imputation of missing values. The users can choose either median or random to impute missing values. 

•	is positive: refers to the relation between the value of the performance measure-ment and the quality of the marker. In the case of direct relationship, the         user should select TRUE such for Multiclass AUC, Area above the cost curve, and Hy-pervolume under manifold.  While in the case of inverse relationship, the         user should select FALSE such as entropy and misClassRate.

•	no.simulations: refers to the number of biomarker candidates that should be sim-ulated to calculate the required p-values, as described in the methods.

•	Conf.Interval: confidence level for the confidence interval (i.e. 95%). Feeding the probabilities for good values of the performance measure into a binomial         distribu-tion, one can also compute confidence intervals for the number of biomarker can-didates with better performance values as described by Klawonn et           al. (2016). 

•	positive class: in the case of Area above the relative cost curve, the user has to specify the positive class in the data set. Otherwise, automatically the         first class, alphabetically, will be chosen.

•	corrected p-value: it is only for the performance part. The p-value could be cor-rected either by controlling the FWER or FDR.
