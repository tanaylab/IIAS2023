# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,Rmd,R:light
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: R 4.2.1
#     language: R
#     name: r_4.2.1
# ---

# # Getting started with naryn

# Naryn is an implementation of a time-series database, for efficient storage, retrieval and analysis of electronic health records (EHR).

# ## Setup

# + [markdown] tags=[]
# ### Install packages if needed

# + vscode={"languageId": "r"}
required_packages <- c("naryn", "dplyr", "purrr", "ggplot2", "survminer", "survival", "cmprsk", "pROC", "xgboost")
for (pkg in required_packages){
    if (!(pkg %in% rownames(installed.packages()))){
        install.packages(pkg)
    }
}

# + [markdown] tags=[]
# ### Load packages

# + vscode={"languageId": "r"}
library(naryn)
library(dplyr)
library(ggplot2)
theme_set(theme_classic())
# -

# ### Download the example database 
#
# Towards this vignette we are going to use a small database which was simulated to include an example of a typical EMR database. It can be downloaded from [here](https://naryn.s3.eu-west-1.amazonaws.com/naryn_example_db.tar.gz) or using the following code:

# + vscode={"languageId": "r"}
if (!dir.exists("sample_db") && !dir.exists("/kaggle/input/simulated-ehr-dataset/sample_db")){
    emr_download_example_data()
}
# -

if (dir.exists("/kaggle/input/simulated-ehr-dataset/sample_db")){
    db_dir <- "/kaggle/input/simulated-ehr-dataset/sample_db"
} else {
    db_dir <- "sample_db"
}

# Note that although smaller than the real database, this example database is still quite large (~1.2GB) and will take a few minutes to download.

# ## Tracks

# The basic element of Naryn is a `track` - a single numerical data element (e.g. Creatinine lab measurment) that is recorded for many patients at various time points. A track can be thought of as a very sparse two-dimensional matrix with a row for each patient in the database, and a column for each timepoint in the resolution of hours. Another way to think of a track is as a table with triplets of patient, timepoint and value: 
#
# |id|time   |lab.CREATININE.BLOOD|
# |--|-------|--------------------|
# |1 |1228730|1.1803077           |
# |1 |1229711|1.1490350           |
# |1 |1230163|1.1682341           |
# |1 |1239512|1.1314427           |
# |1 |1248503|1.2582408           |
# |2 |1274663|0.7628091           |
# |2 |1288516|0.7240421           |
# |2 |1321502|0.7352935           |
# |3 |1222798|0.9173284           |
# |3 |1286821|0.9149746           |
#
# The value can be any numerical value, but it is usually a lab test result, or categorical variable representing a diagnosis. The time is a number representing the number of hours since 1/3/1867 00:00, and the id is a unique identifier for each patient which is defined in a special track called '*patients.dob*', which contains the time of birth for each patient (see below).

# ## Connect to database

# A naryn database is a folder containing a number of tracks. To connect to a database, use the `emr_db.connect` function:

# + vscode={"languageId": "r"}
emr_db.connect(db_dir)
# -

# Now we can use `emr_track.ls` to list all the tracks in the database:
#
#

# + vscode={"languageId": "r"}
head(emr_track.ls())

# + vscode={"languageId": "r"}
# show number of available tracks in database
length(emr_track.ls())
# -

# We can also use `emr_track.ls` to list all the tracks that match a certain pattern:
#
#

# + vscode={"languageId": "r"}
emr_track.ls("CREATININE")
# -

# > Note: Naryn supports connecting to multiple databases at the same time, by giving a vector of paths to `emr_db.connect`, see more in the 'Database' section of the [advanced vignette](https://tanaylab.github.io/naryn/articles/advanced.html).

# ## Extract data from tracks

# We can now go back to the Creatinine example and extract the data from the track. This can be done using the `emr_extract()` function which is the 'Swiss army knife' of the package. It can be used to extract data from a single track, or from multiple tracks, while applying various filters and transformations, but we will start with the simplest example of extracting the data from a single track:

# + vscode={"languageId": "r"}
creatinine_df <- emr_extract("lab.CREATININE.BLOOD")
head(creatinine_df)
# -

# We can see that the data is returned as a data frame with four columns: id, time, ref and 'lab.CREATININE.BLOOD' which contains the values. The id column contains the patient id, the time column contains the time of the measurement in hours since 1/3/1867 00:00, and the value column contains the Creatinine measurment of the patient at each timepoint. 
#
# For information regarding the ref column, see the 'Records and References' section in the [advanced vignette](https://tanaylab.github.io/naryn/articles/advanced.html).
#
# The data frame is sorted by id and time, so that the data for each patient is consecutive. 
#
# You can change the name of the value column using the `names` argument:

# + vscode={"languageId": "r"}
head(emr_extract("lab.CREATININE.BLOOD", names = "Creatinine"))
# -

# ### Track expression 

# The first argument of `emr_extract` was, in the example above, a string representing the name of the track. However, it can also be a *track expression*, which is a string that can contain functions that will be applied to the track after extracting the data. For example, we can extract the lab value multiplied by two:

# + vscode={"languageId": "r"}
creatinine_df1 <- emr_extract("lab.CREATININE.BLOOD * 2", names = "Creatinine")
head(creatinine_df1)
# -

# Note that the functions used in the track expression must be able to be applied to a vector of values and should return a vector of the same length.

# #### Changing time to date

# We can transform the `time` field from the internal `naryn` representation to year, month day and hour using the `emr_time2date()` function:

# + vscode={"languageId": "r"}
creatinine_df %>%
    mutate(emr_time2date(time)) %>%
    head()
# -

# ## Iterators and virtual tracks

# In the previous example we extracted the data from a single track, what happens if we want to extract data from multiple tracks at once? For example, we would want to extract the Creatinine and blood Glucose levels of patients which had Chronic Kidney Disease (CKD), 5 years before their diagnosis. In order to achieve this we need to introduce the concept of *iterators*. 
#
# An iterator is a set of points in the *patient-time space*[^1] that defines the way in which naryn traverses the database. In the Creatinine example, the iterator was simply the set of all *patient-time* points that were included in the track, but now - we would want our point-of-view to be the CKD diagnosis so we would set the iterator to the track of ICD9 code 585 `dx.icd9_585`: 

# + vscode={"languageId": "r"}
ckd_labs <- emr_extract(c("lab.CREATININE.BLOOD", "lab.GLUCOSE.BLOOD"), iterator = "dx.icd9_585")
head(ckd_labs)
nrow(na.omit(ckd_labs))
# -

# We can see that the "lab.CREATININE.BLOOD" and "lab.GLUCOSE.BLOOD" are all NA. This is because the iterator is set to the 
# CKD diagnosis, and the Creatinine and RBC tests were not conducted at the exact time the diagnosis was recorded. In order to get the Glucose/Creatinine test that was done 5 years prior to the diagnosis, we would have to: 
#
# 1. Tell naryn to look back 5 years before the diagnosis
# 2. Tell naryn what to do if there is more than one Creatinine test in the 5 years prior to the diagnosis, for example - take the latest test, or more generally - tell naryn which function to apply to values of the track. 
#
# This can be done using a *virtual track*. A virtual track is a way to tell naryn how to compute the value of a track when at a specific point in the *patient-time space*. It is created using the `emr_virtual_track()` function, which tells naryn how to shift the time for a given track (1, `time.shift`) and which function to apply to the values of the track (2, `func`). 
#
# Note that the time shift is always in reference to the *iterator*. So for example, in our case we would like to look at an Creatinine / Glucose measurement 5 years before the diagnosis. If our iterator is at the time of the CKD diagnosis, we want to look at a time window that starts 5 years before and ends at the time of diagnosis:
#

# + vscode={"languageId": "r"}
emr_vtrack.create("creatinine_5y", "lab.CREATININE.BLOOD", time.shift = c(-years(5), 0), func = "latest")
emr_vtrack.create("glucose_5y", "lab.GLUCOSE.BLOOD", time.shift = c(-years(5), 0), func = "latest")
# -

# Now, we can extract the data from the virtual tracks (instead of the tracks themselves):
#
#

# + vscode={"languageId": "r"}
ckd_labs <- emr_extract(c("creatinine_5y", "glucose_5y"), iterator = "dx.icd9_585", names = c("Creatinine", "Glucose"))
head(ckd_labs, n = 20)
# -

# We would like to also know how long before the diagnosis each test was performed, so we will create an additional virtual track, this time with a function that computes the difference between the time of the diagnosis and the latest blood test. Also, we would like the time difference to be in the resolution of months, so we will give `emr_extract` a track expression that divides the time difference by 30 * 24 (the number of hours in a month), which is equivalent to the `month` function:

# + vscode={"languageId": "r"}
emr_vtrack.create("creatinine_5y_d", "lab.CREATININE.BLOOD", time.shift = c(-years(5), 0), func = "dt2.latest")
emr_vtrack.create("glucose_5y_d", "lab.GLUCOSE.BLOOD", time.shift = c(-years(5), 0), func = "dt2.latest")
ckd_labs <- emr_extract(
    c("creatinine_5y", "glucose_5y", "creatinine_5y_d/month()", "glucose_5y_d/month()"),
    iterator = "dx.icd9_585",
    names = c("Creatinine", "Glucose", "Creatinine_d", "Glucose_d")
)
head(ckd_labs, n = 20)

# -

# Yay! we got what we wanted, but something is still weird - we can see that some patients (and in this sampled database - most patients) have more than one diagnosis of a CKD. Many times, this is indeed the case - a patient can be diagnosed with a disease multiple times, but more commonly in EMR data - the same diagnosis is recorded multiple times, and what we actually want is the earliest diagnosis. In order to achieve that we would have to use *filters*.
#
# [^1]: Actually, *patient-time-reference* space, but this is explained in the advanced vignette.

# ## Filters

# A *filter* is a logical condition that is applied to the *iterator* in order to decide which points to include and which to exclude. In our case we want to include only the earliest diagnosis of CKD, so we would use the `emr_filter.create` function to create a filter that would exclude all patient-time points that had a diagnosis of CKD prior to the current point:

# + vscode={"languageId": "r"}
emr_filter.create("ckd_in_past", "dx.icd9_585", time.shift = c(-years(120), -1))
# -

# We will now apply the filter to our query: 
#
#

# + vscode={"languageId": "r"}
ckd_labs <- emr_extract(
    c("creatinine_5y", "glucose_5y", "creatinine_5y_d/month()", "glucose_5y_d/month()"),
    iterator = "dx.icd9_585",
    names = c("Creatinine", "Glucose", "Creatinine_d", "Glucose_d"),
    filter = "!ckd_in_past"
)

head(ckd_labs)
# -

# Voila! every patient now has only the earliest diagnosis of CKD. 

# ### Value filters 

# Filters can not only be used to exclude or include points by the mere existence of a point in the *patient-time space* (like the previous example), but also by the value of the point. For example, we can create a filter that would include only points where the Glucose was abnormal (say, above 100):

# + vscode={"languageId": "r"}
emr_filter.create("abnormal_glucose", "glucose_5y", val = 100, operator = ">")
ckd_labs_abnormal_glucose <- emr_extract(
    c("creatinine_5y", "glucose_5y", "creatinine_5y_d/month()", "glucose_5y_d/month()"),
    iterator = "dx.icd9_585",
    filter = "!ckd_in_past & glucose_5y",
    names = c("Creatinine", "Glucose", "Creatinine_d", "Glucose_d")
)
head(ckd_labs_abnormal_glucose)

# -

# > Exercise: why do we still get Glucose values of below 100? (hint: `vtrack` function)

# Another example is to include or exclude based on the value of a categorical track. For example, in order to include only patients that were diagnosed with stage iv of CKD (ICD9 code 585.4) we filter "dx.icd9_585" to include only points with value of "14" (see note below):

# + vscode={"languageId": "r"}
emr_filter.create("severe", "dx.icd9_585", val = 14)
ckd_severe <- emr_extract(
    c("creatinine_5y", "glucose_5y", "dx.icd9_585"),
    iterator = "dx.icd9_585",
    filter = "!ckd_in_past & severe",
    names = c("Creatinine", "Glucose", "icd9_585")
)
head(ckd_severe)
# -

# > NOTE: since ICD9 diagnosis codes have a tree like structure, and X.0, X.00 are both valid codes and must be distinguishable, the diagnosis tracks all include a prefix of 1 for the minor code, so X.0 will be translated to icd9_X and a value of 10 will be stored instead of 0.

# ## Time range

# Until now we learned how to use filters to exclude or include datapoints by their value or existence, but what if we want to include only points that are within a certain time range? 
#
# For example, in many EHR databases the diagnoses in beginning of the database do not reflect the actual time the diagnosis was given, or we know that there was a change in the standard care at a certain year. In such cases we can use the `stime` and `etime` arguments which limit the query to points that are between the given start and end times: 

ckd_df <- emr_extract(
    c("creatinine_5y", "glucose_5y", "dx.icd9_585"),
    iterator = "dx.icd9_585",    
    names = c("Creatinine", "Glucose", "icd9_585"),
    stime = emr_date2time(1, 1, 2008), 
    etime = emr_date2time(1, 1, 2022),
)
head(ckd_df)

# ## Define virtual tracks for age and sex

# Given what we learned about iterators, virtual tracks and filters, how would we extract the *age* of a patient at a given timepoint? For example, say that we want to know the age of a patient at the time of their first diagnosis of CKD. We can do this by creating a virtual track that computes the difference between the time of the diagnosis and the time of birth. 
#
# A virtual track needs 4 things - name, source, time shift and function. We should set the name to "age", the source to "patients.dob" and at each point of the iterator we would like to go backward in time a maximal amount (say - 120 years) and compute the time difference between the point and the birth, so we will set the `time.shift` to `c(-years(120), 0)` and the function to `dt2.earliest`:

# + vscode={"languageId": "r"}
emr_vtrack.create("age", "patients.dob", time.shift = c(-years(120), 0), func = "dt2.earliest")
# -

# We can now extract the age of the patient at the time of their first diagnosis of CKD in the resolution of years: 

# + vscode={"languageId": "r"}
age_at_ckd_diag <- emr_extract("age/year()", iterator = "dx.icd9_585", name = "age")
head(age_at_ckd_diag)
# -

# Sex can be defined in a similar way going backward in time and taking the earliest value at the "patients.dob" track:
#
#

# + vscode={"languageId": "r"}
emr_vtrack.create("sex", "patients.dob", time.shift = c(-years(120), 0), func = "earliest")
# -

emr_extract(c("age/year()", "sex"), iterator = "dx.icd9_585", name = c("age", "sex")) %>% 
    mutate(sex = c("male", "female")[sex]) %>% 
    ggplot(aes(x=age, color = sex, group = sex)) + geom_density()

# We can define in the same way filters for other tracks, such as the time of death of a patient or the fact that they registered to the DB or left for good: 

emr_filter.create('born', 'patients.dob', time.shift=c(-120,0)*year())
emr_filter.create('dead', 'patients.dod', time.shift=c(-120,0)*year())
emr_filter.create('registered', 'patients.status.register', time.shift=c(-120,0)*year())
emr_filter.create('left_for_good', 'patients.status.lfg', time.shift=c(-120,0)*year())

# We can use the same approach to define a virtual track for the onset of CKD instead of using the filter we defined above: 

emr_vtrack.create('disease_onset', 'dx.icd9_585', time.shift=c(0,120)*year(), func='earliest.time')

# And finally, we can define a filter for the fact that the patient will have the disease in the future:

emr_filter.create('will_have_disease', 'dx.icd9_585', time.shift=c(0,120)*year())

# ## Incidence rate

# Next, we would want to compute the incidence rate of CKD, for each age and sex. The formula for incidence rate is: 
#
# $$ (New Cases_{a,s}) / (Population_{a,s} x Timeframe) $$
#
# $$a = age, s = sex$$

# The formula may appear straightforward, but it can be challenging to implement in traditional databases, particularly when computing the denominator, which requires counting the age and sex of all patients in the system. An alternative approach, known as the "index-date" method, is often used, in which the calculation is performed only for patients at a specific point in time, rather than the entire database. This approach can result in a small number of cases, making it difficult to stratify the computation based on age and sex.
#
# So how can we use `naryn` to compute this?

# + [markdown] tags=[]
# #### The enumerator 
# -

# We need to count the new cases of CKD in each age group and sex. To do that we will use the function `emr_dist` which can compute the distributions of values of a track. 
#
# `emr_dist` accepts pairs of track expression and breaks, and counts the number of values in each strata, so we can ask it to count the numbers of disease onsets in each age and sex. 

# +
age_breaks <- seq(35, 85, by = 5)
stime <- emr_date2time(1, 1, 2008) # taking 3 years of history
etime <- emr_date2time(1, 1, 2022) # end of db

disease_count <- emr_dist(
        "age/year()", 
        age_breaks, 
        "sex", 
        NULL, 
        iterator = "dx.icd9_585",    
        filter = "!ckd_in_past",
        right = FALSE, 
        dataframe = TRUE, 
        names = c("age", "sex"), 
        stime = stime,
        etime = etime
    ) %>% rename(n_sick = n)
tail(disease_count)
# -

# - Note again that we used the `"!ckd_in_past"` filter in order to make sure that we are going over the *onset* of CKD and not only diagnosis. 
#
# - For sex, we set the breaks to "NULL" in order to implicitly use the possible values of the track. 

# + [markdown] tags=[]
# #### The denominator 
# -

# In order to calculate the denominator, we would have to sum the 'man-months' we have in the system. To do that we would use another type of iterator called `beat iterator`. A beat iterator goes over all the patient in the system every "time beat", for example, in our case, every month. 
#
# We would have to count only patients that did not have the disease yet (although we allow them to get it in the future):

emr_filter.create("has_disease", "dx.icd9_585", time.shift = c(-120, 0) * year())

# In addition, we count only patients that were born, didn't die yet and are currently registered to the system. This would be enforced using filter parameter. 

full_pop <- emr_dist(
    "age/year()", 
    age_breaks, 
    "sex", 
    NULL,
    iterator = month(),
    filter = "born & !dead & registered & !left_for_good & !has_disease", 
    right = FALSE, 
    dataframe = TRUE, 
    names = c("age", "sex"), 
    stime = stime, 
    etime = etime
) %>%
    rename(man_months = n)
head(full_pop)

# #### Computing incidence rate

# Now that we have both the enumerator and the denominator we can compute the incidence rate per 100k:

# +
incidence <- disease_count %>%
    full_join(full_pop, by = c("age", "sex")) %>%
    mutate(
        p = n_sick / man_months,
        incidence = p * 12 * 1e5,
        lower_ci = floor(1000 * (p - qnorm(0.975) * sqrt(p * (1 - p) / man_months)) * 12 * 1e5) / 1000,
        upper_ci = floor(1000 * (p + qnorm(0.975) * sqrt(p * (1 - p) / man_months)) * 12 * 1e5) / 1000,
        sex = factor(c("male", "female")[sex])
    )

head(incidence)
# -

incidence %>% 
    ggplot(aes(x = age, y = incidence, ymin = lower_ci, ymax = upper_ci, group = sex, color = sex)) +
    geom_line() +
    ylab("Incidence per 100k") + 
    xlab("Age") + 
    geom_errorbar(width = 0.3) 

# ## Survival analysis

# Next, we want to estimate the survival of CKD patients between the ages of 60 and 65 (not inclusive). We have a track called "patients.dod" which indicates the death of patients. We will start by computing the time from the disease diagnosis to death:

# create a filter for male patients (val=1) between ages 60 to 65 (time.shift=c(-65, -60))
emr_filter.create("male_60_65", "patients.dob", time.shift = c(-65, -60) * year(), val = 1)

emr_vtrack.create("survival", "patients.dod", time.shift = c(0, 120) * 365 * 24, func = "dt1.earliest")
disease_survival <- emr_extract(c("age/year()", "sex", "survival"),
    iterator = "dx.icd9_585",    
    stime = stime, etime = etime,
    filter = "!ckd_in_past & male_60_65",
    names = c("age", "sex", "survival")
)
head(disease_survival)

# The `NA` values indicate that the patient did not die yet. 
#
# We now need to compute the survival of the "background" patients. We want to go over all the patients at ages 60-65, and compute their time until death ("survival" virtual track). To do that we can use an `extended beat iterator`, which, as its name suggests, is a variation on the beat iterator which aligns the beat to a set of points. In our case we will align the points on their time of birth ("patients.dob") and add the number of years for each patient (e.g. 60) to extract the data at each patient's birthday. 
#
# We will use the same filters we used in the incidence rate analysis above.

gen_pop_survival <- purrr::map_df(60:64, ~ 
    emr_extract(
        c("age/year()", "sex", "survival"),
        iterator = list(.x * year(), "patients.dob"),
        stime = stime, etime = etime,
        filter = "!has_disease & born & !dead & registered & !left_for_good & male_60_65",
        names = c("age", "sex", "survival")
)) 

# We prefer not to use a patient twice (in different ages):

gen_pop_survival <- gen_pop_survival%>%
    sample_n(n()) %>%
    distinct(id, .keep_all = T)
head(gen_pop_survival)

# We can now compute the Kaplan-Meier survival curve. 
#
# We will start by combining the two data frames and reshape it to the format fo the survival function. Note that we apply censoring to reflect the latest update of the database which was at Jan 1, 2021:

survival <- disease_survival %>%
    mutate(cohort = "disease") %>%
    bind_rows(gen_pop_survival %>% mutate(cohort = "general pop")) %>%
    mutate(
        follow_time = ifelse(!is.na(survival), survival, (emr_date2time(1, 1, 2021) - time)),
        status = ifelse(is.na(survival), 0, 1)
    ) %>%
    mutate(follow_time = follow_time / year())
head(survival)

# We can now use `survminer` and `survival` packages to compute and plot the Kaplan-Meier curves:

fit <- survminer::surv_fit(survival::Surv(follow_time, status) ~ cohort, data = survival)
survminer::ggsurvplot(fit, data = survival, xlab = "Time (years)", xlim = c(0, 10))

# ## Time to outcome

# Another type of analysis we can do is compute the distribution of time to outcome. Say, for example, that we want to check if patients with high Creatinine levels are more likely to be diagnosed with CKD. We will again focus our analysis on ages 60-65. 

# We will start by extracting the high and low creatinine tests. We will do it using `emr_screen` which instead of extracting the information from a track expression, returns the patient+time points that satisfy a condition:

creat_high <- emr_screen("lab.CREATININE.BLOOD > 1.5", filter = "!has_disease") %>% distinct(id, .keep_all = TRUE)
creat_low <- emr_screen("lab.CREATININE.BLOOD < 1", filter = "!has_disease") %>% distinct(id, .keep_all = TRUE)
head(creat_high)
head(creat_low)

cohort <- creat_high %>%
    mutate(cohort = "creatinine_high") %>%
    bind_rows(creat_low %>% mutate(cohort = "creatinine_low"))
head(cohort)

# We will now define a virtual track which will compute the time from the current iterator to the onset of the disease and extract it for every Creatinine lab test:

emr_filter.create("age_60_65", "patients.dob", time.shift = c(-65, -60) * year())
emr_vtrack.create(
    "time_to_disease", 
    "dx.icd9_585", 
    time.shift = c(0, 120) * year(), 
    func = "dt1.earliest", 
    filter = "!ckd_in_past"
)    
time_to_outcome <- emr_extract(
    c("age/year()", "sex", "survival", "time_to_disease"),
    iterator = cohort,
    stime = stime, etime = etime,
    filter = "age_60_65",
    names = c("age", "sex", "survival", "time_to_disease")
)
head(time_to_outcome)

# We will encode the event of CKD as `1` and the event of death as `2`:

time_to_outcome <- time_to_outcome %>% 
    mutate(
        follow_time = ifelse(!is.na(survival), survival, (emr_date2time(1, 1, 2021) - time)),
        follow_time = ifelse(!is.na(time_to_disease), time_to_disease, follow_time),
        status = ifelse(!is.na(time_to_disease), 1, ifelse(!is.na(survival), 2, 0))
    ) %>%
    mutate(follow_time = follow_time / year()) %>%
    inner_join(cohort %>% select(id, time, cohort), by = c("id", "time"))
head(time_to_outcome)

# We can now fit and plot a model with a competing risk of death which compares the probability for CKD for high and low Creatinine values:

# + tags=[]
fit <- cmprsk::cuminc(time_to_outcome$follow_time, time_to_outcome$status, group = time_to_outcome$cohort)
survminer::ggcompetingrisks(
    fit = fit, 
    multiple_panels = FALSE, 
    main = "", 
    xlab = "Time (years)", 
    ylab = "Cumulative incidence (%)", 
    xlim = c(0, 5), 
    ylim = c(0, 0.2)
)
# -

# ## Building a CKD predictor

# Now that we looked at the data a bit, we would like to create a model that would try to predict CKD from the EHR data we have 5 year before the diagnosis.  

# #### Sample controls

# One of the important decisions when building a predictor is defining the control. Minimally, we would like to sample patients that are similar on age and sex, but many times we would want to enforce other parameters such as the lack of a disease or number of available tests. 
#
# For this example, we would sample based on age, sex, and the calendaric year. The calendaric year is important in many cases as standard care / distribution of the population changes over time.
#
# We will compute the distribution of age, sex and calendaric year for the patient who had the disease, and then sample control using this distribution:

disease <- emr_extract("dx.icd9_585", filter = "!ckd_in_past") 
cohort_dist <- emr_dist(
    "age/year()", 
    20:90, 
    "sex", 
    NULL, 
    "emr_time2year(EMR_TIME)", 
    2008:2022, 
    iterator = disease, 
    names = c("age", "sex", "year"), 
    dataframe = TRUE, 
    right = FALSE
) %>%
    mutate(
        age = (20:90)[as.numeric(age)], # convert age and year from factor to numeric
        year = (2008:2022)[as.numeric(year)]
    ) %>%
    filter(n > 0)
head(cohort_dist)

# Go over each each combination of age, sex and calendaric year and sample patients:

# +
exclude_patients <- disease$id
cf <- 1
control <- data.frame(id = c(), time = c())

# go over each combination of age, sex and calendaric year
for (i in 1:nrow(cohort_dist)) {
    # exclude the patients with the disease and those that we already selected
    emr_filter.create("fexclude_patients", data.frame(id = exclude_patients, time = emr_date2time(1, 1, 2002)), time.shift = c(-120, 120) * year()) 
    
    # create a filter for the current sex and age
    emr_filter.create("fsex_age", "patients.dob", time.shift = c(-cohort_dist$age[i] - 1, -cohort_dist$age[i]) * year(), val = as.numeric(cohort_dist$sex[i]))
    
    # extract patients at the given age and sex 
    pool <- emr_extract("sex",
        iterator = cohort_dist$age[i] * year(), # this basically says - take only a single point at the current age
        filter = "fsex_age & !fexclude_patients",
        stime = emr_date2time(1, 1, cohort_dist$year[i]),
        etime = emr_date2time(1, 1, cohort_dist$year[i] + 1)
    )
    
    # sample patients+time and add to the control
    control <- rbind(control, pool %>% sample_n(cohort_dist$n[i] * cf))
    
    exclude_patients <- c(disease$id, control$id)
    emr_filter.rm("fexclude_patients")
    emr_filter.rm("fsex_age")
}
# -

head(control)

# #### Extract the data

# Now we can extract the features for the disease and control patient+time points. We would extract all the labs in the sampled database:

labs <- emr_track.ls("^lab.")
length(labs)
head(labs)

# Create for every lab a virtual track which would sample a value (if exists) 5 years prior to the current iterator:

purrr::walk(labs, ~ emr_vtrack.create(paste0("v_", .x), .x, time.shift = c(-5, -1) * year(), func = "sample"))

# Combine the disease and control to a single data frame:

options(emr_max.data.size = 1e09)
target <- disease %>%
    select(id, time) %>%
    mutate(class = 1) %>%
    bind_rows(control %>% select(id, time) %>% mutate(class = 0))
head(target)

# Extract the labs data + age and sex:

data <- emr_extract(c("age/year()", "sex", paste0("v_", labs)), iterator = target, names = c("age", "sex", labs))
head(data)

# #### Separate to train and test

# We will now separate our data to train and test, while maintaining age and sex:

train_test <- data %>%
    select(id, age, sex) %>%
    inner_join(target %>% select(id, class), by = "id") %>%
    mutate(age = cut(age, seq(20, 90, by = 5), right = FALSE, include.lowest = TRUE)) %>%
    group_by(age, sex, class) %>%
    mutate(train = runif(n()) < 0.8) %>%
    ungroup()
head(train_test)

# #### Run a XGBoost predictor

# +
target <- train_test
features <- train_test %>%
    select(id) %>%
    left_join(data) %>%
    select(-id, -time, -ref)

library(xgboost)
train_y <- target[target$train, ]
test_y <- target[!target$train, ]

train_x <- as.matrix(features[target$train, ])
test_x <- as.matrix(features[!target$train, ])

xgb_train <- xgb.DMatrix(data = train_x, label = train_y$class)
xgb_test <- xgb.DMatrix(data = test_x, label = test_y$class)

watchlist <- list(train = xgb_train, test = xgb_test)

# fit XGBoost model and display training and testing data at each round
model <- xgb.train(
    data = xgb_train, 
    booster = "gbtree", 
    objective = "binary:logistic", 
    max.depth = 3, 
    subsample = 0.5, 
    eta = 0.07, 
    min_child_weight = 1, 
    gamma = 0, 
    eval_metric = "auc", 
    watchlist = watchlist, 
    nrounds = 100
)
model
# -

# #### Evaluate the model

# Predict the test data:

test <- test_y %>% mutate(pred_score = predict(model, test_x))

# Calculate ROC:

pROC_obj <- pROC::roc(test$class, test$pred_score,
    smoothed = TRUE,
    # arguments for ci
    ci = TRUE, ci.alpha = 0.9, stratified = FALSE,
    # arguments for plot
    plot = TRUE, auc.polygon = TRUE, max.auc.polygon = TRUE, grid = TRUE,
    print.auc = TRUE, show.thres = TRUE
)
pROC_obj

# Plot the score difference between the classes:

test %>%
    mutate(score = cut(pred_score, c(0, 0.5, 1), right = FALSE, include.lowest = T)) %>%
    count(class, score) %>%
    ggplot(aes(x = score, y = n, fill = factor(class))) +
    geom_bar(stat = "identity", position = "dodge")

# ## Additional examples

# In the next few sections we will present a few additional examples of how to use `naryn` for analysis of EMR data. It is recommended to try and implement a solution to each example before looking at the code below it. 

# ### Extract survival time for patients with pancreatic cancer

# Extract the survival time (time until death) for patients with pancreatic cancer (ICD9 code 157.9):

# Define a filter for just Pancreatic cancer in the 157 icd9 code:

# + vscode={"languageId": "r"}
emr_filter.create("pancreatic_cancer", "dx.icd9_157", val = 19)
# -

# Define a filter for previous pancreatic cancer (sometime in the past):
#
#

# + vscode={"languageId": "r"}
emr_filter.create("pancreatic_cancer_in_past", "dx.icd9_157", val = 19, time.shift = c(-years(120), -1))
# -

# To find the first diagnosis of pancreatic cancer we will go over all 157 diagnosis, filter out those that are not pancreatic cancer and make sure there wasn't a prior diagnosis of pancreatic cancer:
#
#

# + vscode={"languageId": "r"}
pancreatic_cancer_survival <- emr_extract("survival", iterator = "dx.icd9_157", filter = "pancreatic_cancer & !pancreatic_cancer_in_past")
nrow(pancreatic_cancer_survival)
head(pancreatic_cancer_survival)
# -

# Note that NA in survival means that the patient still has not died.
#
# We can now compute the Kaplan-Meier survival curve for pancreatic cancer patients. 
# Censoring is applied to reflect latest update of database (Jan 6, 2022):
#

# + vscode={"languageId": "r"}
pancreatic_cancer_survival <- pancreatic_cancer_survival %>%
    mutate(
        follow_time = ifelse(!is.na(survival), survival, emr_date2time(6, 1, 2022) - time),
        status = ifelse(is.na(survival), 0, 1)
    ) %>%
    mutate(follow_time = follow_time / month())
# -

# Fit the survival curve: 
#
#

# + vscode={"languageId": "r"}
library(survminer)
library(survival)
fit <- survminer::surv_fit(survival::Surv(follow_time, status) ~ 1, data = pancreatic_cancer_survival)
survminer::ggsurvplot(fit, data = pancreatic_cancer_survival, xlab = "Time (months)")
# -

# ### Extract HGB and age for all patients between the ages 60 and 70

# Extract all Hemoglobin tests (HGB), age (in years) and sex for patients between the ages of 60 and 70:

# + vscode={"languageId": "r"}
emr_filter.create("age_60_70", "patients.dob", time.shift = c(-years(70), -years(60)))
hgb_60_70 <- emr_extract(c("age/year()", "lab.HGB", "sex"), iterator = "lab.HGB", filter = "age_60_70", names = c("age", "HGB", "sex"))
# -

# We can now plot the distribution: 
#
#

# + vscode={"languageId": "r"}
hgb_60_70 %>%
    mutate(sex = factor(c("male", "female")[sex])) %>%
    ggplot(aes(x = HGB, color = sex, group = sex)) +
    geom_density()
# -

# ### Extract patients age at time of diagnosis of heart disease that did not have diabetes in the past

# Extract patients age at the time of their first diagnosis of heart disease (diagnosis.411) that did not have diabetes (diagnosis.250) in the past.
#
# This is very similar to the example above, in which we used a filter to exclude patients that had a previous diagnosis of heart disease. We are now going to use a slightly different approach that would first extract a the patient-time points of the heart disease onset and then use an additional query with these points as the iterator to filter out the patients that had diabetes in the past. This approach is useful in cases where we want to use the onset for other purposes, and it is convenient to have it in a separate table.
#
# What would be our iterator? We can iterate the heart disease diagnosis track, and in many cases this would be the right choice. However, sometimes this track can be quite large, for example if it is a very common lab test. In such cases we can iterate over the “patients.dob” track instead, while shifting the time forward by the maximum amount (e.g. 120 years):

# + vscode={"languageId": "r"}
# Find onset of heart disease
emr_vtrack.create("heart_onset", "dx.icd9_411", time.shift = c(0, years(120)), func = "earliest.time")
emr_filter.create("has_heart_disease", "dx.icd9_411", time.shift = c(0, years(120)))

# Retreive the earliest time of heart disease for all patients in db that have a heart disease
heart_onset <- emr_extract("heart_onset", iterator = "patients.dob", filter = "has_heart_disease")
head(heart_onset)
# -

# Now we can use the data frame we created as an iterator:
#
#

# + vscode={"languageId": "r"}
# filter patients that already have diabetes before heart disease
emr_filter.create("has_diabetes", "dx.icd9_585", time.shift = c(-years(120), 0))
heart_onset_no_diabetes <- emr_extract("age/year()", iterator = heart_onset %>% select(id, time = heart_onset), filter = "!has_diabetes", names = "age")

nrow(heart_onset)
nrow(heart_onset_no_diabetes)

head(heart_onset_no_diabetes)
# -

# ### Extract median Hemoglobin for all males between ages 60 to 70

# Extract for each male patient the median Hemoglobin value between ages 60 and 70. 
#
# Again we would use "patients.dob" as our iterator - this is very efficient as every patient would be examined only once, and define a virtual track with a time shift of 60 to 70 years with a function of median: 

# + vscode={"languageId": "r"}
emr_vtrack.create("median_hgb", "lab.HGB", time.shift = c(years(60), years(70)), func = "quantile", params = c(0.5))
# -

# Define a filter for males only: 
#
#

# + vscode={"languageId": "r"}
emr_filter.create("is_male", "patients.dob", val = 1)
# -

# Add another filter for having an HGB test in the relevant ages: 
#
#

# + vscode={"languageId": "r"}
emr_filter.create("has_hgb", "lab.HGB", time.shift = c(years(60), years(70)))
# -

# Extract median_hgb for males only:
#
#

# + vscode={"languageId": "r"}
male_hgb_60_70_q50 <- emr_extract("median_hgb", iterator = "patients.dob", filter = "is_male & has_hgb")
# -

# Compare with female median hgb:
#
#

# + vscode={"languageId": "r"}
female_hgb_60_70_q50 <- emr_extract("median_hgb", iterator = "patients.dob", filter = "!is_male & has_hgb")
plot(density(male_hgb_60_70_q50$median_hgb, na.rm = TRUE), col = "blue", main = "median HGB in ages 60-70")
lines(density(female_hgb_60_70_q50$median_hgb, na.rm = TRUE), col = "red")
# -

# ### Count number of patients by age that were in the system in January 2020

# What does it mean to be "in the system"? Every EHR system would have its own definition of this, but in general - we want the patients that have already been born, have not died yet, have registered with the EHR system and haven't left the system (for good) yet.
#
# We will start by creating a set of filters that will define the above conditions:

# + vscode={"languageId": "r"}
emr_filter.create("born", "patients.dob", time.shift = c(-years(120), 0))
emr_filter.create("dead", "patients.dod", time.shift = c(-years(120), 0))
emr_filter.create("registered", "patients.status.register", time.shift = c(-years(120), 0))
emr_filter.create("left_for_good", "patients.status.lfg", time.shift = c(-years(120), 0))
# -

# We can now use an iterator of a single point while applying all the filters to count the number of patients that were in the system in January 2020:
#
#

# + vscode={"languageId": "r"}
age_dist_2020 <- emr_dist("age/year()", 0:120, "sex", NULL,
    iterator = 1,
    stime = emr_date2time(1, 1, 2020),
    etime = emr_date2time(1, 1, 2020),
    filter = "born & !dead & registered & !left_for_good",
    names = c("age", "sex"),
    dataframe = TRUE,
    right = FALSE
)
head(age_dist_2020 %>% filter(n > 0))

# -

# We can plot the distribution:
#
#

# + vscode={"languageId": "r"}
age_dist_2020 %>%
    mutate(sex = factor(c("male", "female")[sex])) %>%
    mutate(age = as.numeric(age) - 1) %>%
    ggplot(aes(x = age, y = n, colour = sex, fill = sex, group = sex)) +
    geom_col() +
    facet_wrap(~sex) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# -


