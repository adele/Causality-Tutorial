rm(list=ls())

library(FCI.Utils)

########################################################
# Defining true underlying causal diagram and true PAG #
########################################################

allvars <- c("Age", "Sex", "Wealth", "BedNet", "Malaria")
p <- length(allvars)
amat <- matrix(0, p, p)
colnames(amat) <- rownames(amat) <- allvars
amat["Age","Wealth"] <- 0; amat["Wealth","Age"] <- 1; # Age -> Wealth
amat["Sex","Wealth"] <- 0; amat["Wealth","Sex"] <- 1; # Sex -> Wealth

amat["Wealth","BedNet"] <- 0; amat["BedNet","Wealth"] <- 1; # Wealth -> BedNet
amat["Age","Malaria"] <- 0; amat["Malaria","Age"] <- 1; # Age -> Malaria
amat["BedNet","Malaria"] <- 0; amat["Malaria","BedNet"] <- 1; # BedNet -> Malaria

# Adjacency Matrix
print(amat)

# Generating DAG image
renderDAG(amat, fileid = "trueDAG", add_index = FALSE)

lat <- c()
trueDAG <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
dagitty::latents(trueDAG) <- lat
truePAG <- getTruePAG(trueDAG, verbose = FALSE)
true.amat.pag <- truePAG@amat[colnames(amat), colnames(amat)]
true.amat.pag

renderAG(true.amat.pag, add_index = FALSE)

trueMAG <- dagitty::toMAG(trueDAG)
trueImpliedCI <- dagitty::impliedConditionalIndependencies(trueMAG, type = "missing.edge")
trueImpliedCI



########################
# Simulating Variables #
########################

seed <-30868460
set.seed(seed)

check_fit = FALSE

n_individuals <- 50000 # in Mancio Lima Study: 1620

saved_objs_folder <- paste0("./saved_objs/example_malaria_indep/", seed,
                            "_n", n_individuals, "/")
if (!file.exists(saved_objs_folder)) {
  dir.create(saved_objs_folder, recursive=TRUE)
}

#######
# Age #
#######

malaria_dataset <- data.frame(person_id = 1:n_individuals)
malaria_dataset[, "Age"] <- rnorm(n_individuals, 45, 10)
summary(malaria_dataset[, "Age"])

#######
# Sex #
#######
malaria_dataset <- cbind.data.frame(malaria_dataset,
                                    Sex = rbinom(n_individuals, 1, prob = c(0.5)))
table(malaria_dataset$Sex) # 0: female; 1: male

##########
# Wealth #
##########

malaria_dataset[, "Wealth"] <-
  0.3 * malaria_dataset$Sex +
  + 0.15 * malaria_dataset$Age +
  rnorm(n_individuals, 0, 1)

summary(malaria_dataset[, "Wealth"])

if (check_fit) {
  fit_wealth <- lm(Wealth ~ Age + Sex, data=malaria_dataset)
  summary(fit_wealth)
}

##########
# BedNet #
##########

lp_yes = -5 + 0.75 * as.numeric(malaria_dataset$Wealth)

# Simulate binary outcome
p_yes <- 1 / (1 + exp(-lp_yes))
summary(p_yes)

bed_net_use <- rbinom(n_individuals, 1, p_yes) # 0: no; 1: yes
malaria_dataset[, "BedNet"] <- bed_net_use
table(bed_net_use)

if (check_fit) {
  fit_bednet <- glm(BedNet ~ Wealth + Age + Sex,
                  family="binomial", data=malaria_dataset)
  summary(fit_bednet)
}

###########
# Malaria #
###########

eta <- 2.5 +
  - 0.15 * malaria_dataset$Age + # age effect
  - 0.5 * malaria_dataset$BedNet # + # reduced effect for using bed net
mu <- exp(eta)

is_zero <- rbinom(n_individuals, size = 1, prob = 0.25)

nb_values <- rnbinom(n_individuals, size = 2, # lower = more over dispersion
                     mu = mu)

malaria_dataset[, "Malaria"] <- ifelse(is_zero == 1, 0, nb_values)
table(malaria_dataset[, "Malaria"] )

if (check_fit) {
  fit_malaria <- pscl::zeroinfl(Malaria ~ Age +  BedNet + Wealth + Sex | 1,
                                dist = "negbin", data = malaria_dataset)
  summary(fit_malaria)
}

##############################
# Part III: Causal Discovery #
##############################

############################################
# Computing Conditional Independence Tests #
############################################

library(FCI.Utils)

# In case you want to run in parallel
library(doFuture)
library(future.apply)
n_cores <- 12 # change to the number of cores available
#plan("multisession", workers = n_cores)
plan("multicore", workers = n_cores) # forking


dat <- malaria_dataset[,colnames(true.amat.pag)]
head(dat)

run_dcfci = TRUE

vars_names <- colnames(dat)
covs_names = c()
indepTest <- mixedCITest

suffStat <- getMixedCISuffStat(dat, vars_names, covs_names)
vars_df <- dat[,vars_names, drop=FALSE]

suffStat$types["Malaria"] <- "count"
suffStat$count_regr <-  "simplzeroinfl" # for zero-inflated negative binomial distr.
suffStat$verbose <- TRUE
suffStat$comb_p_method <- "min"
suffStat$citestResults <- NULL

saved_suffstat <- paste0(saved_objs_folder, "suffStat.RData")
if (!file.exists(saved_suffstat)) {
  if (run_dcfci) {
    citestResults <- getAllCITestResults(vars_df, indepTest, suffStat,
                                         computeProbs = TRUE,
                                         eff_size = 0.001)
  } else {
    citestResults <- getAllCITestResults(vars_df, indepTest, suffStat)
  }
  suffStat$citestResults <- citestResults
  save(suffStat, file=saved_suffstat)
} else {
  load(saved_suffstat)
}

# Checking faithfulness degree
faithf_degree <- getFaithfulnessDegree(true.amat.pag, suffStat$citestResults,
                                       alpha = 0.05, bayesian = run_dcfci)
faithf_degree$f_citestResults
faithf_degree$faithful_pprop
faithf_degree$faithful_bprop

subset(faithf_degree$f_citestResults, bf == FALSE)

###########################
# Running the FCI from    #
# pcalg R package         #
###########################

library(pcalg)

# Under faithfulness, FCI recovers the true PAG
labels <- colnames(suffStat$dataset)
alpha <- 0.05
fci_fit <- fci(suffStat, indepTest, alpha, labels)
renderAG(fci_fit@amat, add_index = FALSE)

shd_fci <- shd_PAG(amat.trueP = true.amat.pag, amat.estP = fci_fit@amat)
shd_fci

faithf_fci <- getFaithfulnessDegree(fci_fit@amat, suffStat$citestResults,
                                    bayesian = FALSE)
faithf_fci


#######################
# Running dcFCI from  #
# dcFCI R package     #
#######################

if (run_dcfci) {
  library(dcFCI)
  library(dplyr)

  saved_fit_dcfci <-  paste0(saved_objs_folder, "fit_dcfci.RData")
  if (!file.exists(saved_fit_dcfci)) {
  # dcFCI is more robust under unfaithfulness
    fit_dcfci <- dcFCI(suffStat, indepTest,
                     labels=labels, alpha=0.01,
                     sel_top = 1,
                     prob_sel_top = TRUE,
                     verbose = TRUE,
                     run_parallel = TRUE)
    save(fit_dcfci, file=saved_fit_dcfci)
  } else {
    load(saved_fit_dcfci)
  }

  top_dcpag <- fit_dcfci$allPAGList[[1]]
  renderAG(top_dcpag$amat.pag, add_index = FALSE)
  head(fit_dcfci$top_scoresDF)

  shd_dcfci <- shd_PAG(amat.trueP = true.amat.pag, amat.estP = top_dcpag$amat.pag)
  shd_dcfci

  faithf_dcfci <- getFaithfulnessDegree(top_dcpag$amat.pag, suffStat$citestResults,
                                      bayesian = FALSE)
  faithf_dcfci
}

####################################################################
# Part IV: Effect Identification from the Markov Equivalence Class #
####################################################################

####################################################################
# Effect Identification via Generalized Adjustment Criterion (GAC) #
####################################################################


# Checking identifiability of P(y|do(x)) through generalized adjustment criterion (GAC).

# estPAG <- fci_fit@amat # FCI PAG
estPAG <- top_dcpag$amat.pag # dcFCI PAG
renderAG(estPAG, add_index = TRUE) # showing variable indices

labels <- colnames(estPAG)
x = which(labels == "BedNet")
y =  which(labels == "Malaria")

adj <- adjustment(amat = estPAG, amat.type = "pag",
                  x = x, y = y, set.type = "all")

if (length(adj) > 1) {
  print(paste0("P(", labels[y], "|do(", labels[x], ")) is identifiable via adjustment over:"))
  sapply(adj, function(x) { paste(labels[x], collapse=",") } ) # List of all sets admissible for adjustment
} else {
  print(paste0("P(", labels[y], "|do(", labels[x], ")) is not identifiable via adjustment."))
}


# Checking whether a particular set is admissible for adjustment.

zlabels <- c("Age", "Sex")
z <- which(labels %in% zlabels)
print(paste0("Is the effect of ", labels[x], " on ", labels[y], " identifiable through adjustment over {", paste(zlabels, collapse=","), "}? "))

gac_out <- gac(estPAG, x, y, z, type = "pag")
is_admissible <- gac_out$gac

print(paste0("Answer: ", is_admissible))

############################################
# Effect Identification via CIDP Algorithm #
############################################

if (!require(PAGId)) {
  library(devtools, warn.conflicts=F, quietly=T)
  devtools::install_github("adele/PAGId", dependencies=TRUE)
}

library(PAGId)

y = "Malaria"
x = "BedNet"
z = c()

retPAG <- CIDP(estPAG, x, y, z, verbose = FALSE)

# This shows the steps taken by the CIDP algorithm
# By substitution and simplication, we will get the same adjustment formula
print(paste0("Is P(", y, "|do(", x, ")) identifiable? ", retPAG$id))
if (retPAG$id) {
  print(paste0("P(", y, "|do(", x, ")) = "))
  print(retPAG$Qexpr)
}



