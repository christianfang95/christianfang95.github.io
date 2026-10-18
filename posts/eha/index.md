---
layout: post
title: Event History Analysis in R
subtitle: Kaplan–Meier and Discrete-Time Logistic Regression
date: '2026-10-18'
description: >-
  Learn how to run a basic survival/event-history analysis in R using
  Kaplan-Meier and logistic regression.
categories:
  - statistics
  - demography
  - R
toc: true
engine: markdown
code-fold: true
format:
  gfm:
    variant: +yaml_metadata_block
    wrap: preserve
execute:
  echo: true
  warning: false
  message: false
output-file: index.md
---


- [Introduction](#introduction)
  - [Setup](#setup)
- [Survival and the hazard](#survival-and-the-hazard)
- [Kaplan–Meier survival](#kaplanmeier-survival)
- [Person-Period files](#person-period-files)
- [Discrete-time EHA using logistic regression](#discrete-time-eha-using-logistic-regression)
  - [First model: the baseline hazard only](#first-model-the-baseline-hazard-only)
  - [Predicted hazards using `marginaleffects`](#predicted-hazards-using-marginaleffects)
  - [From predicted hazards back to survival](#from-predicted-hazards-back-to-survival)
  - [Comparing logistic regression with Kaplan–Meier](#comparing-logistic-regression-with-kaplanmeier)
  - [Adding covariates](#adding-covariates)
  - [Why adjusted survival curves?](#why-adjusted-survival-curves)
  - [Creating counterfactual prediction data](#creating-counterfactual-prediction-data)
  - [Predicting adjusted hazards with `marginaleffects`](#predicting-adjusted-hazards-with-marginaleffects)
  - [From adjusted hazards to adjusted survival](#from-adjusted-hazards-to-adjusted-survival)
  - [Unadjusted versus adjusted group differences](#unadjusted-versus-adjusted-group-differences)
  - [Comparing probabilities directly with `marginaleffects`](#comparing-probabilities-directly-with-marginaleffects)
  - [Differences at particular ages](#differences-at-particular-ages)
- [Modeling the baseline hazard differently](#modeling-the-baseline-hazard-differently)
  - [Logit versus complementary log-log](#logit-versus-complementary-log-log)
- [The logic of the analysis](#the-logic-of-the-analysis)

# Introduction

This tutorial introduces the basics of event history analysis (EHA)/survival analysis using the `scania` dataset from the `eha` package to study mortality after age 50.

The main idea is to approach the same survival process in two ways. We first estimate survival non-parametrically using Kaplan–Meier (i.e., we do not impose a functional form on the survival process and just examine what the raw data shows). We then restructure the data into person-period format and estimate the age-specific hazard using logistic regression. From those predicted hazards we reconstruct a survival curve and compare it directly with Kaplan–Meier.

After establishing this connection do we add explanatory variables and calculate adjusted survival curves.

## Setup

Install the packages if necessary:

``` {r}
#| eval: false

install.packages(c(
  "eha",
  "survival",
  "dplyr",
  "tidyr",
  "ggplot2",
  "marginaleffects"
))
```

Load the packages:

``` {r}
library(eha)
library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(marginaleffects)
```

Load the Scania data:

``` {r}
data(scania)

head(scania)
summary(scania)
```

The most important variables include:

- `id`: person identifier
- `enter`: age at which observation begins
- `exit`: age at death or censoring
- `event`: whether death occurred at `exit`
- `sex`: sex
- `ses`: socioeconomic status

In these data, observation starts at age 50.

Suppose a person has:

``` text
enter = 50
exit  = 59.24
event = 1
```

This means that the person is observed from age 50 and dies at age 59.24.

If instead:

``` text
enter = 50
exit  = 80
event = 0
```

the person is observed until age 80 but does not die while under observation. The observation is **right-censored**.

# Survival and the hazard

Two quantities are central to event history analysis.

The **survival function**

$$
S(t)=P(T>t)
$$

gives the probability of remaining event-free beyond time $t$.

In our example:

> What is the probability that somebody who has reached age 50 survives beyond age 50?

The **hazard** asks a different question. In discrete time:

$$
h_j =
P(\text{death in interval }j
\mid
\text{alive at the beginning of interval }j).
$$

For example:

> Among people who are alive at age 70, what is the probability that they die between ages 70 and 71?

The two quantities are directly connected.

If the hazards in successive intervals are

$$
h_1,h_2,\ldots,h_j,
$$

then survival through interval $j$ is

$$
S_j =
\prod_{k=1}^{j}(1-h_k).
$$

This relationship will become important when we estimate the hazard with logistic regression.

# Kaplan–Meier survival

We begin with the Kaplan–Meier estimator.

Create a survival object:

``` {r}
surv <- Surv(
  time = scania$enter,
  time2 = scania$exit,
  event = scania$event
)
```

Because we specify both entry and exit time, this takes the form

``` r
Surv(enter, exit, event)
```

and describes the period during which each person is at risk.

Estimate overall survival:

``` {r}
km <- survfit(
  surv ~ 1,
  data = scania
)

km
```

Plot the survival curve:

``` {r}
plot(
  km,
  xlab = "Age",
  ylab = "Survival probability",
  main = "Kaplan–Meier survival curve"
)
```

The curve starts at 1 at age 50 and declines as deaths occur.

We can request survival probabilities at selected ages:

``` {r}
summary(
  km,
  times = c(60, 70, 80)
)
```

For example, the survival probability at age 70 answers:

> Among people who survived to age 50, what proportion are estimated to survive to at least age 70?

Kaplan–Meier estimates this directly from the observed survival times. We have not yet specified a regression model.

# Person-Period files

To estimate a **discrete-time event history model**, we restructure the data.

The original data have one row per person.

For discrete-time analysis we instead want one row for every time interval during which a person is at risk.

We will use one-year age intervals.

Someone who dies at age 55.4 will therefore contribute observations approximately like:

``` text
id    age    death
1      50      0
1      51      0
1      52      0
1      53      0
1      54      0
1      55      1
```

The variable `death` now indicates whether the event occurs **during that particular interval**.

Thus:

- `death = 0` means that the person survives the interval.
- `death = 1` means that the person dies during the interval.

Create the person-period data:

``` {r}
scania_pp <- scania |>
  mutate(
    event = as.integer(event),
    n_intervals = if_else(
      event == 1,
      ceiling(exit) - 50,
      floor(exit) - 50
    )
  ) |>
  filter(n_intervals > 0) |>
  uncount(
    n_intervals,
    .id = "period",
    .remove = FALSE
  ) |>
  mutate(
    age = 49 + period,
    death = as.integer(
      event == 1 &
        period == n_intervals
    )
  )
```

Inspect the resulting file:

``` {r}
scania_pp |>
  select(id, age, death, sex, ses) |>
  head(20)
```

We can also inspect one person who dies:

``` {r}
example_id <- scania_pp |>
  filter(death == 1) |>
  slice(1) |>
  pull(id)

scania_pp |>
  filter(id == example_id) |>
  select(id, age, death)
```

There are multiple observations for this person but only one event, in the final interval.

# Discrete-time EHA using logistic regression

We now estimate the same mortality process using logistic regression.

Because each row represents one interval during which a person is at risk, our dependent variable is

$$
P(death_{it}=1).
$$

More precisely, the model estimates

$$
P(
\text{death during interval }t
\mid
\text{alive at the start of interval }t
).
$$

This is the **discrete-time hazard**.

## First model: the baseline hazard only

Before adding explanatory variables such as sex or socioeconomic status, we estimate a model containing only **age**.

For the first model, it is useful to treat age as categorical:

``` {r}
m0 <- glm(
  death ~ factor(age),
  data = scania_pp,
  family = binomial(link = "logit")
)

summary(m0)
```

Why use `factor(age)`?

Because we do not yet want to impose a particular shape on the baseline hazard.

Instead, the model estimates a separate mortality probability for each age interval.

Conceptually, this model is estimating

$$
h_{50}, h_{51}, h_{52}, \ldots
$$

where, for example, $h_{70}$ is the probability of dying between approximately age 70 and 71 among people who are alive at age 70.

## Predicted hazards using `marginaleffects`

The coefficients from a logistic regression are expressed as log odds.

Those are usually not the most useful quantities for interpretation.

We can instead use `marginaleffects` to obtain predicted probabilities.

``` {r}
hazard0 <- predictions(
  m0,
  newdata = datagrid(
    age = sort(unique(scania_pp$age))
  ),
  type = "response"
)

hazard0
```

The variable `estimate` now contains the estimated probability of dying during each age interval.

``` {r}
hazard0 |>
  as.data.frame() |>
  select(age, estimate)
```

An estimate such as

``` text
age      = 59
estimate = 0.033
```

would mean:

> Among people who have survived to age 59, approximately 3.33% are predicted to die during the following age interval.

This is a **hazard probability**, not a survival probability.

## From predicted hazards back to survival

We can now convert the predicted hazards into a survival curve.

Recall:

$$
S_j =
\prod_{k=1}^{j}(1-h_k).
$$

Using R, this cumulative product is calculated with `cumprod()`.

``` {r}
surv0 <- hazard0 |>
  as.data.frame() |>
  mutate(age = as.numeric(levels(age))) |>
  arrange(age) |>
  mutate(
    survival = cumprod(1 - estimate),
    age_end = age + 1
  )
```

Inspect the result:

``` {r}
surv0 |>
  select(
    age,
    estimate,
    age_end,
    survival
  )
```

For every interval:

- `estimate` is the probability of dying during the interval.
- `survival` is the probability of surviving all intervals up to its end.

Plot the model-based survival curve:

``` {r}
ggplot(
  surv0,
  aes(
    x = age_end,
    y = survival
  )
) +
  geom_line(linewidth = 1) +
  labs(
    x = "Age",
    y = "Survival probability",
    title = "Survival derived from discrete-time logistic regression"
  ) +
  theme_minimal()
```

## Comparing logistic regression with Kaplan–Meier

We can now compare the two approaches.

First extract the Kaplan–Meier estimates:

``` {r}
km_summary <- summary(km)

km_df <- data.frame(
  age = km_summary$time,
  survival = km_summary$surv
)
```

Then plot the Kaplan–Meier curve together with the survival curve reconstructed from the discrete-time logistic regression:

``` {r}
ggplot() +
  geom_step(
    data = km_df,
    aes(
      x = age,
      y = survival
    ),
    linewidth = 0.9
  ) +
  geom_line(
    data = surv0,
    aes(
      x = age_end,
      y = survival
    ),
    linewidth = 1,
    color = 'red',
    alpha = 0.7
  ) +
  labs(
    x = "Age",
    y = "Survival probability",
    title = "Kaplan–Meier (black) and discrete-time logistic (red) survival"
  ) +
  theme_minimal()
```

The curves essentially the same survival process, which illustrates an important point:

> Kaplan–Meier and discrete-time logistic regression are not estimating completely different phenomena.

Both describe the same underlying process.

Kaplan–Meier estimates survival directly. The logistic regression estimates the **interval-specific hazard**, after which we reconstruct survival by accumulating those hazards:

$$
\hat h_t
\quad\longrightarrow\quad
\prod_{k \leq t}(1-\hat h_k)
\quad\longrightarrow\quad
\hat S(t).
$$

We can also compare the Kaplan-Meier survival estimates at ages 60, 70 and 80 directly to those from the logistic regression.

``` {r}
summary(
  km,
  times = c(60, 70, 80)
)
```

``` {r}
surv0 |>
  filter(age_end %in% c(60, 70, 80)) |>
  select(age_end, survival)
```

As we see, the estimates are basically identical.

So why use regression at all?

Because regression allows us to introduce explanatory variables.

## Adding covariates

We now extend the model by adding sex and socioeconomic status.

``` {r}
m1 <- glm(
  death ~ factor(age) + sex + ses,
  data = scania_pp,
  family = binomial(link = "logit")
)

summary(m1)
```

Age still represents the baseline hazard.

But the model now also allows mortality to differ according to:

- sex;
- socioeconomic status.

For example, the coefficient for sex describes how the log odds of dying during an interval differ between men and women at the same age and holding SES constant.

We could convert coefficients into odds ratios:

``` {r}
exp(coef(m1))
```

However, predicted probabilities and survival curves are generally easier to interpret.

## Why adjusted survival curves?

Suppose men and women differ in their socioeconomic composition.

An ordinary Kaplan–Meier comparison between men and women reflects both:

1.  differences associated with sex itself;
2.  differences in the distribution of SES between the two groups.

Regression allows us to ask a different question:

> What would predicted survival look like for men and women after accounting for differences in SES?

We can answer this using model-based standardization.

Rather than choosing one arbitrary SES category, we retain the observed SES distribution and predict survival under two counterfactual scenarios:

- everyone is treated as male;
- everyone is treated as female.

Everything else remains as observed.

## Creating counterfactual prediction data

First create one baseline observation per person:

``` {r}
baseline <- scania_pp |>
  distinct(
    id,
    sex,
    ses
  )
```

Now create every combination of:

- person;
- age;
- hypothetical sex.

``` {r}
prediction_data <- baseline |>
  select(id, ses) |>
  crossing(
    age = sort(unique(scania_pp$age)),
    sex_scenario = levels(scania_pp$sex)
  ) |>
  mutate(
    sex = factor(
      sex_scenario,
      levels = levels(scania_pp$sex)
    )
  )
```

Each person now appears at every age under both sex scenarios.

Their observed SES is retained.

This is important because we want to compare sex while averaging over the same SES distribution.

## Predicting adjusted hazards with `marginaleffects`

Generate predicted interval-specific hazards:

``` {r}
adjusted_hazards <- predictions(
  m1,
  newdata = prediction_data,
  type = "response"
) |>
  as.data.frame()
```

We now have a predicted probability of death for:

- every individual;
- every age;
- each hypothetical sex.

For example:

``` {r}
adjusted_hazards |>
  select(
    id,
    age,
    sex_scenario,
    ses,
    estimate
  ) |>
  head()
```

The `estimate` column is again the discrete-time hazard.

## From adjusted hazards to adjusted survival

For each hypothetical person, convert their sequence of hazards into survival:

``` {r}
adjusted_survival_individual <- adjusted_hazards |>
  arrange(
    sex_scenario,
    id,
    age
  ) |>
  group_by(
    sex_scenario,
    id
  ) |>
  mutate(
    survival = cumprod(1 - estimate),
    age_end = age + 1
  ) |>
  ungroup()
```

We now have one predicted survival trajectory for every individual under each sex scenario.

To obtain the **standardized survival curve**, average these predicted survival probabilities across individuals:

``` {r}
adjusted_survival <- adjusted_survival_individual |>
  group_by(
    sex_scenario,
    age_end
  ) |>
  summarise(
    survival = mean(survival),
    .groups = "drop"
  )
```

Plot the resulting curves:

``` {r}
ggplot(
  adjusted_survival,
  aes(
    x = age_end,
    y = survival,
    linetype = sex_scenario
  )
) +
  geom_line(linewidth = 1) +
  labs(
    x = "Age",
    y = "Adjusted survival probability",
    linetype = "Sex",
    title = "SES-adjusted survival curves"
  ) +
  theme_minimal()
```

These curves can be interpreted as:

> The predicted survival curves for men and women if both groups had the same SES distribution as the population in the data.

## Unadjusted versus adjusted group differences

We can also look at the observed Kaplan–Meier curves separately by sex:

``` {r}
km_sex <- survfit(
  Surv(enter, exit, event) ~ sex,
  data = scania
)

plot(
  km_sex,
  xlab = "Age",
  ylab = "Survival probability",
  lty = 1:2
)
```

These are **unadjusted** survival curves.

They describe the actual survival experience of men and women in the data.

By contrast, the regression-based curves from the previous section are adjusted for SES.

Conceptually:

``` text
Kaplan–Meier by sex
        ↓
Observed / unadjusted difference

Regression + standardization
        ↓
Model-based difference after adjusting for SES
```

This distinction is one of the main reasons to move from Kaplan–Meier estimation to regression.

## Comparing probabilities directly with `marginaleffects`

The regression model can also be used to compare predicted hazards directly.

For example:

``` {r}
avg_comparisons(
  m1,
  variables = "sex",
  type = "response"
)
```

This converts the sex coefficient into an average difference in the predicted probability of death during an interval.

For instance, an estimate of

``` text
0.00161
```

would correspond to roughly a **0.2 percentage-point difference in the interval-specific probability of death**, averaged across the observations used in the calculation.

This is often easier to interpret than an odds ratio.

## Differences at particular ages

Because logistic regression is nonlinear, the same coefficient can correspond to different probability differences at different ages.

We can therefore estimate contrasts at specific ages:

``` {r}
comparisons(
  m1,
  variables = "sex",
  newdata = datagrid(
    age = c(60, 70, 79),
    ses = unique
  ),
  type = "response"
)
```

The general principle is:

> Do not assume that a constant logit coefficient implies a constant difference in probabilities.

# Modeling the baseline hazard differently

For the first model we used:

``` r
factor(age)
```

This is useful because it imposes essentially no functional form on the age pattern of mortality.

But it requires many parameters.

Once the basic logic is clear, we could use a smoother specification.

For example, a quadratic:

``` {r}
m_quad <- glm(
  death ~ age + I(age^2) + sex + ses,
  data = scania_pp,
  family = binomial(link = "logit")
)
```

Or, usually more flexibly, a spline:

``` {r}
m_spline <- glm(
  death ~ splines::ns(age, df = 4) + sex + ses,
  data = scania_pp,
  family = binomial(link = "logit")
)
```

The difference is that these models impose some smoothness on the baseline hazard.

## Logit versus complementary log-log

So far we have used:

``` r
family = binomial(link = "logit")
```

This gives a discrete-time logistic event history model.

Another common specification is:

``` r
family = binomial(link = "cloglog")
```

For example:

``` {r}
m_cloglog <- glm(
  death ~ factor(age) + sex + ses,
  data = scania_pp,
  family = binomial(link = "cloglog")
)
```

The complementary log-log model has a particularly close connection to continuous-time proportional hazards models.

For short intervals and relatively uncommon events, logit and complementary log-log models often produce similar predicted hazards.

# The logic of the analysis

The whole exercise can now be understood as one sequence:

$$
\text{survival data}
\rightarrow
\text{Kaplan–Meier}
$$

followed by

$$
\text{person-period data}
\rightarrow
\text{logistic model of the hazard}
\rightarrow
\hat h_t
\rightarrow
\prod_t (1-\hat h_t)
\rightarrow
\hat S(t).
$$

With no covariates, this gives us another way to describe the same basic survival process estimated by Kaplan–Meier.

Once covariates are introduced, the regression model allows us to move beyond description and calculate survival curves under different covariate scenarios.

The key insight is therefore not simply that logistic regression can be used for EHA.

It is that a discrete-time survival process can be decomposed into a sequence of conditional binary probabilities:

$$
P(event_t=1 \mid event_{t-1}=0).
$$

Model those probabilities, and we model the hazard.

Accumulate them, and we recover survival.
