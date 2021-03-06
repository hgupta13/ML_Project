#We are using the data from Romero, Mauricio; Sandefur, Justin; Sandholtz, Wayne, 2018,
#"Partnership Schools for Liberia", https://doi.org/10.7910/DVN/5OPIYU, Harvard Dataverse, V4.

#The data was used for the paper "Outsourcing Education: Experimental Evidence from Liberia" by the same set of authors.

#In our project we use the publically available data provided by the authors. The entire dataset is available at the above link 
#and is also available on Google Drive.

#The code for creating the sample.csv file is described below

#Import necessary packages
library(haven)
library(dplyr)


#Import raw data from the publically available data
#This data has been preprocessed
#See "Readme.pdf" for details

#Baseline data for student controls
Student_baseline <-read_dta("~/PSL Endline Dataverse Files/Analysis/usedata_public/baseline/Student.dta")

#Endline data for student outcomes
Student_endline <- read_dta("~/PSL Endline Dataverse Files/Analysis/usedata_public/endline/Student.dta")

#school characterisitcs
school <-  read_dta("~/PSL Endline Dataverse Files/Analysis/usedata_public/SchoolControls.dta")

#Funding data

funding <-read_dta("~/PSL Endline Dataverse Files/Analysis/usedata_public/baseline/FundingWide.dta")


#We select the relevant student and school controls
# The controls are the ones described in Table A.5 of the online appendix of the paper

#Student Controls
Student_baseline = Student_baseline %>% 
  select(age, gender, IndexAssets,grade_2015,schoolid_nopii, groupid_nopii,
           uniqueid_nopii, location, W)

#School level controls
school = school %>% 
  select(schoolid_nopii,groupid_nopii,Average2015,FacilitiesIndex,Rural, TimeToBank,
         treatment, turnedpsl)


#We choose the outcome data from the tests conducted at the end of the school year
#Outcome data
outcomes = Student_endline %>% select(uniqueid_nopii,schoolid_nopii,
                                      groupid_nopii,Interviewed, treatment,
                                      IRT_English, IRT_Abstract,
                                      IRT_math, IRT_Composite,)

#Join the three datasets to create our sample
sample = inner_join(Student_baseline,outcomes,
                    by=c("uniqueid_nopii", "schoolid_nopii","groupid_nopii"))

sample = left_join(sample,school, by = c("schoolid_nopii","groupid_nopii"))
sample = left_join(sample,funding)


#Select relevant covariates. See overleaf for details
sample = sample %>% select(uniqueid_nopii,schoolid_nopii,groupid_nopii, #Identifiers
                           age, gender, grade_2015,IndexAssets, #Student Characteristics
                           Average2015, FacilitiesIndex, Rural, TimeToBank, nTb_a_1, #school characterisitcs
                           treatment.x,IRT_Composite)

colnames(sample) = c("studentid","schoolid","groupid","age","gender","grade","student_wealth",
                     "enrollment_2015","school_facilities","rural","time_to_bank", "NGO_donations"
                     ,"W","IRT_score")


write.csv(sample, "sample.csv")
