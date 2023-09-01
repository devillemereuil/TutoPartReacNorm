conflicted::conflicts_prefer(dplyr::filter,
                             tidyr::unpack,
                             brms::ar,
                             dplyr::lag,
                             tidyr::extract,
                             dplyr::count,
                             purrr::set_names)

## Setting up the "working directory" using the here library
here::i_am("tuto_discrete_environment.R")
library(here)

## Loading the required packages
# We will use the commodity functions provided by the many packages of the tidyverse
# If you are not familiar with the tidyverse, a very nice introduction to it can be found here:
# https://r4ds.hadley.nz/
library(tidyverse)

# For the implementation of the models to study reaction norms, we will use brms and Bayesian statistics.
# This is, for several reasons:
# 1. It allows to run all the different models (including the non-linear ones) within the same framework
# 2. MCMC is a robust algorithm to run those models
# 3. Bayesian statistics allow for quantifying quite naturally the uncertainty around derived parameters
#    such as the variances we want to compute here
library(brms)
library(coda)   # Coda is a nice MCMC utility library

# We need the "reacnorm" branch of QGglmm (not released yet)
remotes::install_github("devillemereuil/QGglmm@reacnorm")
library(QGglmm)

# Finally, we need this library for the rowVars() function used below
library(matrixStats)

## Setting up a seed for reproducibility
# Set seed
seed <- 2023
# This will ensure the code yields the same results
set.seed(seed)

# NOTE Mention something about Bessel's correction somewhere

## ------------------ Loading the data ----

## Loading the simulated data using the Performance Curve scenario ----
tbl_data <- read_csv(here("Data_discrete.csv"),
                     col_types = cols(Name_Env   = col_character(),
                                      Env        = col_double(),
                                      Individual = col_character(),
                                      Phen       = col_double()))

## Exploring the data ----
# The data contains four columns
# - Name_Env: an arbitrary name for the environment
# - Env: the environmental value at which the Individual was measured
# - Individual: which individual the measured was performed on
# - Phen: phenotypic value
tbl_data

## Plotting the data ----
# A quick plot of the individual reaction norms
p_rn <-
    ggplot(tbl_data) +
    geom_line(aes(x = Env, y = Phen, group = Individual, colour = Individual)) +
    theme(legend.position = "none") +
    xlab("Environment") + ylab("Phenotype")

plot(p_rn)

## ------------------ Formatting and transforming the data ----

## Generating the squared env column
tbl_data <-
    tbl_data %>%
    mutate(Env_Sq = (Env - mean(Env))^2)
# NOTE mean-centring to avoid correlation with slope

## ------------------ Parameters for running the models ----

# Total number of iterations
n_iter <- 4000
# Number of iterations that will be discarded for the warm-up
n_warm <- 1000
# Thinning interval (we keep each n_thin iteration to save memory)
n_thin <- 1
# Number of independent MCMC chains
n_chains <- 4

## ------------------ A study of reaction norm using a polynomial function ----

# We will first focus on modelling the reaction norm using a polynomial (here quadratic) function,
# because it is the most widespread approach and can used to tell a lot of things regarding the shape
# of the reaction norm, e.g. even in absence of another justified parametric model.

## Running the Reaction Norm model ----
# Setting up the formula beforehand (e.g. to use sparse arg to spare memory)
# form_rn <- brmsformula(Phen ~ Env + Env_Sq + (1 + Env + Env_Sq | Individual),
#                        sparse = TRUE)
# prior <- prior(normal(0, 10), class = "b")
# model_p2 <-
#     brm(formula     = form_rn,
#         data        = tbl_data,
#         save_pars   = save_pars(group = FALSE),
#         chains      = n_chains,
#         cores       = n_chains,
#         seed        = seed,
#         prior       = prior,
#         iter        = n_iter,
#         warmup      = n_warm,
#         thin        = n_thin)
# # Now saving the model output
# saveRDS(model_p2, file = here("model_p2_ds.rds"))
# The previous commands can be commented out and the model just loaded with this command
model_p2 <- readRDS(here("model_p2_ds.rds"))

## Checking the model ----
# - The convergence of the MCMC can be check using Rhat (should be < 1.01)
# - The length (sufficient number of effective iterations) can be checked with
#   the ESS (Effictive Sample Size), which should be > 500
# For more information about Bayesian statistics and MCMC, you can look into this tutorial:
# https://devillemereuil.legtux.org/wp-content/uploads/2021/09/tuto_en.pdf
summary(model_p2)
# plot(model_p2)

# We can add the predictions of the model on top of the reaction norms
tbl_p2_mod <-
    tbl_data |>
    mutate(Predict = predict(model_p2, re_formula = NA) |>
                     as_tibble()) |>
    unpack(Predict) |>
    select(Env,
           Predict = Estimate,
           Predict_Low = Q2.5,
           Predict_Up  = Q97.5) |>
    summarise(across(starts_with("Predict"), mean),
              .by = Env)
# The reason we need to take the mean is that, for some reason, brms predicts slightly
# different values for individuals in the same environment, even though we told it to
# ignore the random effect (re_formula = NA).

p_rn_p2 <-
    p_rn +
    geom_ribbon(data = tbl_p2_mod,
                mapping = aes(x = Env, ymin = Predict_Low, ymax = Predict_Up),
                alpha = 0.2) +
    geom_line(data = tbl_p2_mod,
              mapping = aes(x = Env, y = Predict),
              linewidth = 1)

plot(p_rn_p2)
# The fit is not extremely great, we will be able to quantify that.
# Note however, that the coefficients can still be interpreted as the average slope and curvature of
# the reaction norms. The black line here is the function that results from those average parameters,
# but it does not reflect the average local phenotype due to the bad fit.


## Getting the information about the design matrix X of the model ----
# Computing the design matrix, it corresponds to the (fixed) realisation of the random vector x in Eq. 13
design_mat <- model.matrix(Phen ~ Env + Env_Sq, data = tbl_data)

# Computing the environmental variance-covariance matrix (used in Eq. 14)
cov_X <- cov(design_mat)

# Computing the average of powers of X (will be used to compute Eq. 23)
M <- (1 / nrow(design_mat)) * t(design_mat) %*% design_mat

## Extracting the estimates from the model ----
# Getting the fixed effects from the model (with the whole posterior distribution)
theta_p2 <- fixef(model_p2, summary = FALSE)
# Getting the error variance-covariance matrix S_theta
S_theta_p2 <- vcov(model_p2)
# Getting the G-matrix from the random effects variances-covariances (with the whole posterior again)
# We use VarCorr() to extract the full posterior of the G-matrix as a 3-dimension array,
# then format it as a list of G-matrices
G_Theta_p2 <-
    VarCorr(model_p2, summary = FALSE)[["Individual"]][["cov"]] |>
    apply(1, \(mat) { mat }, simplify = FALSE)


## Computing V_Plas and related components ----

# Setting up the correcting variance accounting for the uncertainty issue (Eq. 43)
var_uncert_p2 <- sum(cov_X * S_theta_p2)
# NB: In a Bayesian setting, such correction is not completely enough, because the uncertainty
# is impacting both due to the uncertainty for the given dataset (as for frequentist inference)
# and due to the fact that we will account for the uncertainty in the posterior distribution.
# As such, the unbiased V_Plas should be V(design_mat %*% theta) - c * var_correct, where c is between 1 and 2,
# depending on the prior on fixed effects. Problem is, we do not know c!
# However, this is not too much of an issue because:
# 1. If we correct by var_correct, we are sure to remove the issue due to accounting for the posterior,
#    which is already quite a big part of the problem.
# 2. Since var_correct tends toward 0 when the sample size increases, our estimator is "convergent",
#    i.e. V(design_mat %*% theta) tends toward V_Plas when the sample size increases.
# 3. A consquence of both 1. and 2. is that the remaining bias will be extremely small,
#    and thus inconsequential for almost all practical applications. This can be seen here:
var_uncert_p2

# Computing V_Plas (while accounting for the whole posterior distribution)
var_plas_p2 <- apply(theta_p2, 1, \(th) t(th) %*% cov_X %*% th) - var_uncert_p2
# Trace and posterior distribution of V_Plas
plot(as.mcmc(var_plas_p2))
# Posterior mean and credible intervals for V_Plas
mean(var_plas_p2)                  # 0.075
HPDinterval(as.mcmc(var_plas_p2))  # 0.063 - 0.088
# Note that var_uncert is extremely small compared to the value of var_plas (see point 3 above)
var_uncert_p2 / mean(var_plas_p2)

# We can do the same separately for the slope
var_sl_p2 <- (theta_p2[ , "Env"]^2) * var(design_mat[ , "Env"])
pi_sl_p2  <- var_sl_p2 / var_plas_p2
plot(as.mcmc(pi_sl_p2))
mean(pi_sl_p2)                 # 0.24
HPDinterval(as.mcmc(pi_sl_p2)) # 0.18 - 0.30

# We can do the same separately for the curvature
var_cv_p2 <- (theta_p2[ , "Env_Sq"]^2) * var(design_mat[ , "Env_Sq"])
pi_cv_p2  <- var_cv_p2 / var_plas_p2
plot(as.mcmc(pi_cv_p2))
mean(pi_cv_p2)                 # 0.76
HPDinterval(as.mcmc(pi_cv_p2)) # 0.70 - 0.82

# We can make a nice graph to represent those values
tbl_vplas_p2 <-
    tibble(Var_Plas = var_plas_p2,
           Var_Sl   = var_sl_p2,
           Var_Cv   = var_cv_p2) |>
    pivot_longer(starts_with("Var"),
                 names_to = "Parameter",
                 values_to = "Variance") |>
    mutate(Parameter = factor(Parameter, levels = c("Var_Sl", "Var_Cv", "Var_Plas")))

p_var_plas_p2 <-
    ggplot(tbl_vplas_p2) +
    geom_violin(aes(x = Parameter, y = Variance))
plot(p_var_plas_p2)

## Computing V_Gen and related components ----

# We can compute the total genetic variance using the Eq. 23 in the article.
# However, here, we will use a simpler version of this equation. It can be shown that using the M matrix
# of the average for each exponentiation of the environment defined above, we have :
# V_Gen = V_A = colMeans(design_mat)^t %*% G_Theta_p2 %*% colMeans(design_mat) + Tr(ϴ %*% cov_X)
#             = sum(M * G_Theta_p2)
# which greatly simplifies the computation !
# Accounting for the whole posterior distribution, this would look like the following:
var_gen_p2 <- map_dbl(G_Theta_p2, \(mat_g) { sum(M * mat_g) })
# NB: In the case of a polynomial function, V_Gen = V_A:
var_a_p2   <- var_gen_p2

# We can decompose the impact of each component of the G-matrix on the total genetic variance
gamma <- map2_dfr(G_Theta_p2, var_a_p2, \(mat_g, va) { diag(M * mat_g) / va })
gamma[["Int_Env_Sq"]] <-
    map2_dbl(G_Theta_p2, var_a_p2,
             \(mat_g, va) { 2 * M["(Intercept)", "Env_Sq"] * mat_g["Intercept", "Env_Sq"] / va })
# Note that the gamma can be negative, or superior to 1, yet, they do sum to 1:
colMeans(gamma)
apply(gamma, 1, sum)
# What we can get from the values of rg here is that the most important source of genetic variance
# comes from the variance in the intercept (e.g. "pure" genetic variance among individuals),
# while a strong constrain on the genetic variance arises from the negative covariance between the
# intercept and the curvature of the reaction norm.

# Together, the variance in the intercept & curvature accounts for 83% of the variance
gamma |> select(-Env) |> apply(1, sum) |> mean()
# The slope (which is neither correlated with the intercept, nor the curvature) accounts for 17%
mean(gamma[["Env"]])

# Those parameters are important, because they condition the parameters that most likely to have an
# impact on the evolution of the trait (here, rather the intercept or curvature than the slope)

## Computing V_Res and retrieving the total variance ----

# Getting the residual variance
# Do not forgot the ^2 at the end to compute a variance, not a standard deviation
var_res_p2 <- VarCorr(model_p2, summary = FALSE)[["residual__"]][["sd"]][ , 1]^2

# Computing the total variance
var_tot_p2 <- var_plas_p2 + var_gen_p2 + var_res_p2
mean(var_tot_p2)

# Comparison with the sample phenotypic variance
var(tbl_data[["Phen"]])
# Very close to V_tot

## ------------------ A study of reaction norm using the character-state model ----

## Running the Character-State model ---
# Setting up the formula beforehand (e.g. to use sparse arg to spare memory)
# form_cs <- brmsformula(Phen ~ 0 + Name_Env + (0 + Name_Env | Individual),
#                        sparse = TRUE)
# model_cs <-
#     brm(formula     = form_cs,
#         data        = tbl_data,
#         save_pars = save_pars(group = FALSE),
#         chains      = n_chains,
#         cores       = n_chains,
#         seed        = seed,
#         iter        = n_iter,
#         warmup      = n_warm,
#         thin        = n_thin)
# # Now saving the model output
# saveRDS(model_cs, file = here("model_cs_ds.rds"))
# The previous commands can be commented out and the model just loaded with this command
model_cs <- readRDS(here("model_cs_ds.rds"))

## Checking the model ----
# Everything should look OK
summary(model_cs)
# plot(model_cs)

# We can add the predictions of the model on top of the reaction norms
tbl_cs_mod <-
    tbl_data |>
    mutate(Predict = predict(model_cs, re_formula = NA) |>
                     as_tibble()) |>
    unpack(Predict) |>
    select(Env,
           Predict = Estimate,
           Predict_Low = Q2.5,
           Predict_Up  = Q97.5) |>
    summarise(across(starts_with("Predict"), mean),
              .by = Env)

p_rn_cs <-
    p_rn +
    geom_pointrange(data = tbl_cs_mod,
                    mapping = aes(x = Env, y = Predict, ymin = Predict_Low, ymax = Predict_Up))

plot(p_rn_cs)
# The fit is not extremely great, we will be able to quantify that.
# Note however, that the coefficients can still be interpreted as the average slope and curvature of
# the reaction norms. The black line here is the function that results from those average parameters,
# but it does not reflect the average local phenotype due to the bad fit.

## A word about Bessel's correction ---
# The variance is defined as the average of squared contrast to the mean: mean((x - mean(x))^2).
# However, when computing the variance, the mean itself is unknown. This creates a bias in the estimation
# of the mean, which is easily fixed by dividing by N - 1 rather than N when computing the mean.
# This is known as Bessel's correction V(x) = (N / (N - 1)) * mean((x - mean(x))^2).

# Let's illustrate this:
N <- 10
x <- rnorm(N)
mean((x - mean(x))^2)
(N / (N - 1)) * mean((x - mean(x))^2)
var(x)
# Note that var() is using Bessel's correction, a sensible default to avoid the bias mentioned above.

# However, we now have a problem. What happens when we compute the variance across all environments,
# or across all the individuals? Then, we won't have the same output, even though we should in theory:
N_ind <- 100
var(x)
var(rep(x, each = N_ind))
# Note that the values are different... But we can express the "shortest" one as the "longest one:
((N - 1)/N) * ((N * N_ind) / ((N * N_ind - 1))) * var(x)

# When we computed V_plas, we used var() and cov() on the design_matrix directly,
# so we nrow(design_mat) = 1000 values:
nrow(design_mat)
# But we have only 10 environments, and here this is what we will use:
ncol(fixef(model_cs, summary = FALSE))
# The gap between 10 and 1000 with Bessel correction will be a problem, so we need to compute everything
# with the same number of "points". There are two ways, remove Bessel correction systematically, or express
# the variances on the same number of "points" as above. We will use this second method, since it's the
# easiest to implement:
correct_Bessel_issue <- function(var,
                                 N1  = ncol(fixef(model_cs, summary = FALSE)),
                                 N2  = nrow(design_mat)) {
    ((N1 - 1)/N1) * (N2 / (N2 - 1)) * var
}

## Computing V_plas ---

# Getting the uncertainty on the parameters
var_uncert_cs <-
    vcov(model_cs) |>
    diag() |>
    mean()

# Computing V_plas
var_plas_cs <-
    fixef(model_cs, summary = FALSE) |>
    rowVars() |>
    correct_Bessel_issue()
# Correcting for the uncertainty
var_plas_cs <- var_plas_cs - var_uncert_cs

# Estimates of V_plas
mean(var_plas_cs)
HPDinterval(as.mcmc(var_plas_cs))

## Computing V_gen ---

# Getting the genetic variance in each environment
var_gen_cs <-
    model_cs |>
    as_tibble() |>
    select(starts_with("sd_Individual")) |>
    mutate(across(everything(), \(x) x^2)) |>
    as.matrix() |>
    rowMeans()

# Estimates of V_gen
mean(var_gen_cs)
HPDinterval(as.mcmc(var_gen_cs))

## Computing V_res and V_tot ---

# Getting the residual variance
var_res_cs <- VarCorr(model_cs, summary = FALSE)[["residual__"]][["sd"]][ , 1]^2

# Now, we can compute the total variance
var_tot_cs <- var_plas_cs + var_gen_cs + var_res_cs
mean(var_tot_cs)
var(tbl_data[["Phen"]])

## Comparison with the polynomial model and computation of R² ---

# Comparison of V_plas from the polynomial and character-state models
mean(var_plas_p2)
mean(var_plas_cs)
# We can see that V_plas measured from the p2 model is lower than that of the CS model,
# this is due to imperfected modelling of the reaction norm shape by the quadratic function.

# In fact, we can measure how well (or badly) we fitted the "true" reaction norm by
# comparing both variance:
R2_mod_p2 <- var_plas_p2 / var_plas_cs
mean(R2_mod_p2)
HPDinterval(as.mcmc(R2_mod_p2))
# Here the quadratic model p2 retrieves 71% of the variance arising from the reaction norm.
# It's not incredibly bad, but it's not great either!

# Actually, we could do the same kind of comparison for the other components
mean(var_gen_p2)
mean(var_gen_cs)
# Again, the genetic variance retrieved by the p2 model is lower, because it does not
# manages to fit the true shape of the reaction norm well enough.

# Of course, this is offset by an increase in the "unexplained" residual variance in p2
mean(var_res_p2)
mean(var_res_cs)
# Note the huge difference in residual variance here!

# Graphical comparison of the models
# Polynomial is in blue, character-state is in red
p_rn_comp <-
    p_rn +
    geom_ribbon(data = tbl_p2_mod,
                mapping = aes(x = Env, ymin = Predict_Low, ymax = Predict_Up),
                alpha = 0.2,
                fill = "#0055FF") +
    geom_line(data = tbl_p2_mod,
              mapping = aes(x = Env, y = Predict),
              linewidth = 1,
              colour = "#0055FF") +
geom_pointrange(data = tbl_cs_mod,
                mapping = aes(x = Env, y = Predict, ymin = Predict_Low, ymax = Predict_Up),
                colour = "#C71414")

plot(p_rn_comp)
# We can see how the quadratic function "misses" the height of the peak, while the character-state
# closely follows the true average shape of the reaction norms

## ------------------ A study of reaction norm using a non-linear model ----

## Running a non-linear model ----
# Setting up the formula beforehand (e.g. to use sparse arg to spare memory)
# form_nl <- brmsformula(Phen ~ cmax * exp(
#                                     - exp(rho * (Env - xopt) - 6) -       # Gompertz part
#                                         sigmagaus * (Env - xopt)^2        # Gaussian part
#                                 ),
#                        cmax + xopt ~ 1 + (1 | ID1 | Individual),
#                        rho + sigmagaus ~ 1,
#                        nl = TRUE,
#                        sparse = TRUE)
# prior_nl <- prior(uniform(0, 100), nlpar = "rho", lb = 0, ub = 100) +
#             prior(uniform(0, 10), nlpar = "sigmagaus", lb = 0, ub = 10)
# model_nl <-
#     brm(formula     = form_nl,
#         data        = tbl_data,
#         save_pars = save_pars(group = FALSE),
#         prior       = prior_nl,
#         chains      = n_chains,
#         cores       = n_chains,
#         seed        = seed,
#         init        = rep(list(list(b_cmax      = array(data = 1),
#                                     b_xopt      = array(data = 0.9),
#                                     b_rho       = array(data = 8),
#                                     b_sigmagaus = array(data = 0.4))), 4),
#         iter        = n_iter,
#         warmup      = n_warm,
#         thin        = n_thin)
# # Now saving the model output
# saveRDS(model_nl, file = here("model_nl_ds.rds"))
# The previous commands can be commented out and the model just loaded with this command
model_nl <- readRDS(here("model_nl_ds.rds"))

## Checking the model ----
# Everything should look OK
summary(model_nl)
# plot(model_nl)

# We can add the predictions of the model on top of the reaction norms
tbl_nl_mod <-
    tbl_data |>
    mutate(Predict = predict(model_nl, re_formula = NA) |>
                     as_tibble()) |>
    unpack(Predict) |>
    select(Env,
           Predict = Estimate,
           Predict_Low = Q2.5,
           Predict_Up  = Q97.5) |>
    summarise(across(starts_with("Predict"), mean),
              .by = Env)

p_rn_nl <-
    p_rn +
    geom_ribbon(data = tbl_nl_mod,
                mapping = aes(x = Env, ymin = Predict_Low, ymax = Predict_Up),
                alpha = 0.2) +
    geom_line(data = tbl_nl_mod,
              mapping = aes(x = Env, y = Predict),
              linewidth = 1)

plot(p_rn_nl)
# The fit is much better that for the quadratic function, isn't it?

## Getting estimates ----

# Getting estimates
theta_nl <- fixef(model_nl, summary = FALSE)
colnames(theta_nl) <- str_remove(colnames(theta_nl), "_Intercept")
# Setting up the posterior distribution as a list, it'll be convenient later
theta_nl <- apply(theta_nl, 1, \(vec) { vec }, simplify = FALSE)

# Getting the G-matrix
G_Theta_nl <-
    VarCorr(model_nl, summary = FALSE)[["Individual"]][["cov"]] |>
    apply(1, \(mat) { mat }, simplify = FALSE)
G_Theta_nl <- map(G_Theta_nl,
                  \(mat) { colnames(mat) <- rownames(mat) <- str_remove(colnames(mat), "_Intercept"); mat })

## Computing V_Plas and related components ----

# Generating a function for the Gaussian-Gompertz shape
expr <- expression(
    cmax * exp(
        - exp(rho * (x - xopt) - 6) -
            sigmagaus * (x - xopt)^2
    ))
gausgompz <- QGrn_generate_shape(expr, c("cmax", "xopt", "rho", "sigmagaus"))

# Environment vector
env <- unique(tbl_data[["Env"]])

# Now we can compute the posterior distribution of V_plas
# We will use the QGrn_vplas function from the QGglmm package
# Note the fixed = c(3, 4) argument, that is informing QGrn_vplas that
# the parameters rho (third) and sigmagauss (fourth) are missing from G_Theta (no genetic variation)
var_plas_nl <-
    map2_dbl(theta_nl, G_Theta_nl,
             \(th, g) { QGrn_vplas(env      = env,
                                   shape    = gausgompz,
                                   theta    = th,
                                   G_theta  = g,
                                   fixed    = c(3, 4)) },
             .progress = TRUE)
mean(var_plas_nl)
HPDinterval(as.mcmc(var_plas_nl))
# Note how close our value using the CS model was to that one:
mean(var_plas_cs)
# In a way, the character-state model really is very robust, because of its lack of
# parametric assumptions regarding the shape of the reaction norm. However, the other side
# of the coin is that it is more limited in its predictive ability when evolution occurs
# across environments. Here, we have access to a parametric curve for the reaction norm, which is
# more actionable.

# Of course, the fit of this non-linear model is just as good as the CS model
R2_mod_nl <- var_plas_nl / var_plas_cs
mean(R2_mod_nl)
# An impressive >99%, of course due to the fact that the curve used here is exactly the same
# as the curve used to generate the data!

# We can have a look at the difference between the average for each environment and
# the function of the average parameters, using the internal rn_avg_e() function of QGlmm
avg_e_nl <-
    map2(theta_nl, G_Theta_nl,
         \(th, g) { map_dbl(env,
                            \(e) QGglmm:::rn_avg_e(e        = e,
                                                   shape    = gausgompz,
                                                   theta    = th,
                                                   G_theta  = g,
                                                   fixed    = c(3, 4))) },
         .progress = TRUE) |>
    map(enframe, name = "Env", value = "Value") |>
    list_rbind() |>
    mutate(Env = env[Env])

# Computing estimates and 95% CI for the graph
tbl_avg_e <-
    avg_e_nl |>
    summarise(Estimate = mean(Value),
              CI_Low = HPDinterval(as.mcmc(Value))[1],
              CI_Up = HPDinterval(as.mcmc(Value))[2],
              .by = Env)

p_rnavg_nl <-
    p_rn_nl +
    geom_ribbon(data = tbl_avg_e,
                mapping = aes(x = Env, ymin = CI_Low, ymax = CI_Up),
                colour = "#FF2E2E",
                alpha = 0.2) +
    geom_line(data = tbl_avg_e,
              mapping = aes(x = Env, y = Estimate),
              colour = "#FF2E2E",
              linewidth = 1)
plot(p_rnavg_nl)
# This time, it is a bit closer than the example in the article

## Computing V_Gen and related components ----

# Computing the V_gen by environment
# We will use the QGrn_vgen function from the QGglmm package
# (Note the average = FALSE, so that we can obtain the variance at each environment)
var_gen_e_nl <-
    map2(theta_nl, G_Theta_nl,
             \(th, g) { QGrn_vgen(env      = env,
                                  shape    = gausgompz,
                                  theta    = th,
                                  G_theta  = g,
                                  fixed    = c(3, 4),
                                  average  = FALSE) },
             .progress = TRUE)
# Computing the total (average) V_gen
var_gen_nl <- map_dbl(var_gen_e_nl, mean)
mean(var_gen_nl)
HPDinterval(as.mcmc(var_gen_nl))
# Again, very close to the one of the character-state model
mean(var_gen_cs)

# Computing the additive genetic variance V_A
# First, we need to compute the gradient of our function for the reaction norm.
# Fortunately, QGglmm provides a helper function to do that
d_gausgompz <- QGrn_generate_gradient(expr, c("cmax", "xopt"), c("cmax", "xopt", "rho", "sigmagaus"))

# Computing V_A by environment with gamma-decomposition
# We will use the QGrn_vgen function from the QGglmm package
# Note that we provide the gradient of our function, as computed above
# (Note the average = FALSE, so that we can obtain the variance at each environment)
var_a_e_nl <-
    map2(theta_nl, G_Theta_nl,
             \(th, g) { QGrn_va(env      = env,
                                d_shape  = d_gausgompz,
                                theta    = th,
                                G_theta  = g,
                                fixed    = c(3, 4),
                                average  = FALSE) },
             .progress = TRUE)
var_a_nl <- map_dfr(var_a_e_nl, colMeans)
mean(var_a_nl[["V_A"]])
HPDinterval(as.mcmc(var_a_nl[["V_A"]]))

# Note that the character-state cannot provide this value for V_A, because it lacks the information
# about the true pathway between the trait and parameters. Nevertheless, we could have obtained an estimate # of V_A is we had access to the relatedness between individuals in our setting. However, such a model
# would prove quite large to fit, and thus slow to execute.

# We can then look at the gamma-decomposition, to see which parameters are important in
# generating some additive genetic variance
mean(var_a_nl[["Gamma:cmax"]])                      # 0.64
HPDinterval(as.mcmc(var_a_nl[["Gamma:cmax"]]))
mean(var_a_nl[["Gamma:xopt"]])                      # 0.49
HPDinterval(as.mcmc(var_a_nl[["Gamma:xopt"]]))
mean(var_a_nl[["Gamma:cmax-xopt"]])                 # -0.13
HPDinterval(as.mcmc(var_a_nl[["Gamma:cmax-xopt"]]))

# In our case, most of the additive genetic variation comes from cmax, then xopt, with a negative
# co-contribution of cmax and xopt. This might seem counter-intuitive, since the covariance between
# cmax and xopt is positive in the G_Theta matrix
map_dbl(G_Theta_nl, \(g) g[1, 2]) |> mean()
# However here, the contribution of gradient of the reaction norm according to xopt is often negative,
# while cmax is always positive. As a result of this difference in sign of the gradients, the
# co-contribution of cmax and xopt is, on average, negative. This means that (counter-intuitively?)
# increasing the covariance between cmax and xopt would result in a decrease of the trait V_A.


## Plotting the (additive) genetic variances across environments ----

# Formatting V_Gen_e and V_A_e for a plot
tbl_vg_e <-
    bind_rows(
        var_gen_e_nl |> map(enframe, name = "Env", value = "Value") |> list_rbind(),
        var_a_e_nl |> map("V_A") |> map(enframe, name = "Env", value = "Value") |> list_rbind(),
        .id = "Variance"
    ) |>
    mutate(Variance = recode(Variance, `1` = "V_Gen", `2` = "V_A"),
           Env = env[Env],
           Variance = as_factor(Variance))

p_vg_e <-
    ggplot(tbl_vg_e) +
    geom_violin(aes(x = as_factor(Env), y = Value, fill = Variance),
                scale = "width") +
    ylab("Variance") + xlab("Environment")
plot(p_vg_e)
# Note that a few values are quite divergent for V_gen, but recall that we have 12000 iterations,
# so a few weird samples are to be expected.

## Computing V_Res and V_Tot ----

# Getting the residual variance
var_res_nl <- VarCorr(model_nl, summary = FALSE)[["residual__"]][["sd"]][ , 1]^2

# Computing the total variance
var_tot_nl <- var_plas_nl + var_gen_nl + var_res_nl
mean(var_tot_nl)
# Again, very close to the actual phenotypic variance
var(tbl_data[["Phen"]])
