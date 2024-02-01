####################
# Helper Functions #
####################

# A <- f(Ua) = Ua
# D <- f(Ucd, Ud) = beta_D.Ucd * Ucd + Ud
# B <- f(A, D, Ub) = beta_B.A * A - beta_B.D * D + Ub
# C <- f(Ucd, B, Uc) = beta_C.B * B - beta_C.Ucd * Ucd + Uc
# betaList is a list with entries for:
# beta_D.Ucd, beta_B.A, beta_B.D, beta_C.Ucd, and beta_C.B
# if betaList is NULL, then coefficients are randomly selected
getIVSEM <- function(betaList=NULL) {
  allvars <- c("A", "B", "C", "D", "Ucd")
  p <- length(allvars)
  beta <- matrix(0, p, p)
  colnames(beta) <- rownames(beta) <- allvars
  topolOrd <- c( "Ucd", "A", "D", "B", "C")

  if (is.null(betaList)) {
    betaList <- as.list(sample(c(-1, 1), p, replace = T) * runif(p, 0.3, 0.7))
    names(betaList) <- c("beta_D.Ucd", "beta_B.A", "beta_B.D", "beta_C.Ucd", "beta_C.B")
  }

  beta["Ucd","D"] <- betaList$beta_D.Ucd; # coeff of Ucd in the function for d
  beta["A", "B"] <- betaList$beta_B.A; # coeff of a in the function for b
  beta["D", "B"] <- betaList$beta_B.D; # coeff of d in the function for b
  beta["Ucd", "C"] <- betaList$beta_C.Ucd; # coeff of Ucd in the function for c
  beta["B", "C"] <- betaList$beta_C.B; # coeff of b in the function for c

  # Note that beta is triangular, implying that the SEM is recursive (or acyclic)
  beta <- beta[topolOrd, topolOrd]
  lat <- c("Ucd")

  return(list(beta=beta, lat=lat))
}

# The SEM is constructed as
# Y = beta Y + eps, where
# beta is a matrix of the coefficients for all variables (V and U)
generateDatasetFromSEM <- function(beta, lat, N) {
  p <- ncol(beta)
  ident <- diag(1, p, p)
  colnames(ident) <- rownames(ident) <- colnames(beta)
  # Or, similarly,
  # Y = (I - Beta)^{-1} eps
  # For Gaussian Y and errors eps independent from each other (we have all the Us),
  # we just need to simulate Y from a multivariate distribution
  # of mean zero and cov Sigma=I
  IminusBinv <- ginv(ident - beta)
  Sigma = t(IminusBinv) %*% ident %*% IminusBinv

  valR <- matrixcalc::is.symmetric.matrix(Sigma) &&
    matrixcalc::is.positive.definite(Sigma, tol=1e-8)
  if (!valR) {
    stop("This SEM generates a non-positive definite covariance matrix. Try another SEM.")
  }

  dat <-  MASS::mvrnorm(N, rep(0, p), Sigma, empirical = FALSE)
  dat <- as.data.frame(dat)
  colnames(dat) <- colnames(beta)
  head(dat)

  uvars <- which(colnames(dat) %in% lat)
  dat <- dat[,-uvars]

  return(dat)
}

# Randomly generate a linear SEM following a dagitty DAG, adag,
# and then draw samples from it.
generateDatasetFromDAG <- function(adag, N, ntries=30) {
  done <- FALSE
  tries <- 0
  obs.dat <- NULL
  while (!done && tries <= ntries) {
    done <- tryCatch(
      {
        obs.dat <- dagitty::simulateSEM(adag, b.lower = -0.6, b.upper = 0.6, N=N)
        R <- cor(obs.dat)
        valR <- matrixcalc::is.symmetric.matrix(R) &&
          matrixcalc::is.positive.definite(R, tol=1e-8)
        valR
      }, error=function(cond) {
        message(cond)
        FALSE
      })
    tries <- tries + 1
  }
  return(obs.dat)
}

dagittyOracleCI <- function(x, y, S, suffStat) {
  g <- suffStat$g
  labels <- names(g)
  if (dagitty::dseparated(g, labels[x], labels[y], labels[S])) {
    return(1)
  } else {
    return(0)
  }
}

# receives a dagitty g of type "dag" or "mag" and returns
# the true PAG as an pcalg fci object
getTruePAG <- function(g, verbose = FALSE) {
  indepTest <- dagittyOracleCI
  if (graphType(g) == "dag") {
    g <- dagitty::toMAG(g)
  }
  suffStat <- list(g=g)
  truePag <- pcalg::fci(suffStat,
                        indepTest = indepTest,
                        labels= names(suffStat$g), alpha = 0.9999,
                        verbose = verbose)
  return(truePag)
}



# A -> B -> C; B <- D <- Ucd -> C
getIVGraph <- function() {
  allvars <- c("A", "B", "C", "D", "Uab", "Ucd")
  p <- length(allvars)
  amat <- matrix(0, p, p)
  colnames(amat) <- rownames(amat) <- allvars
  amat["Uab","A"] <- 0; amat["A","Uab"] <- 1; # Uab -> a
  amat["Uab","B"] <- 0; amat["B","Uab"] <- 1; # Uab -> b
  amat["B","C"] <- 0; amat["C","B"] <- 1; # b -> c
  amat["D","B"] <- 0; amat["B","D"] <- 1; # d -> b
  amat["Ucd","C"] <- 0; amat["C","Ucd"] <- 1; # Ucd -> c
  amat["Ucd","D"] <- 0; amat["D","Ucd"] <- 1; # Ucd -> c

  lat <- c("Uab", "Ucd")
  adag <- pcalg::pcalg2dagitty(amat, colnames(amat), type="dag")
  dagitty::latents(adag) <- lat
  coordinates(adag) <- list( x=c(A=0.5, Uab=0.65, B=0.80, C=1, D=0.88, Ucd=1),
                        y=c(A=0.15, Uab=0.05, B=0.15, C=0.15, D=0, Ucd=0) )


  idag <- igraph_from_graphNel(as(t(amat), "graphNEL"), lat)

  return(list(amat.dagitty = adag, # adjacency matrix as a dagitty object
              amat.igraph = idag, # adjacency matrix as an igraph object
              amat.pcalg = amat  # adjacency matrix as a matrix using pcalg notation
              ))
}


# n: number of nodes
# e_prob: probability of having an edge between any two nodes
# l_prob: proportion of the nodes to be considered latent (unmeasured variables)
getRandomGraph <- function(n, e_prob, l_prob, force_nedges=TRUE)  {
  if (force_nedges) {
    min_nedges <- floor((n*(n-1)/2)*(e_prob * 1))
  } else {
    min_nedges <- 0
  }
  nedges <- 0
  while (nedges <= min_nedges) {
    dag <- pcalg::randomDAG(n, e_prob, V = paste0("V", sprintf(
      paste0("%0", nchar(n), "d"), seq(from = 1, to = n))))

    amat.dag <- (as(dag, "matrix") > 0) * 1
    adag <- pcalg::pcalg2dagitty(amat.dag, colnames(amat.dag), type="dag")
    lat <- sample(colnames(amat.dag), l_prob * n)

    dagitty::latents(adag) <- lat

    magg <- dagitty::toMAG(adag)
    nedges <- nrow(edges(magg))
  }

  return(adag)
}

# types can be: png, pdf, or svg
#' @importFrom rsvg rsvg_png
#' @importFrom DOT dot
#' @export renderAG
renderAG <- function(amat, output_folder=NULL, fileid=NULL, type="png",
                     width=NULL, height=NULL, labels=NULL, add_index=TRUE) {
  if (is.null(labels)) {
    labels <- colnames(amat)
  }

  if (is.null(output_folder)) {
    output_folder = "./tmp/"
  }

  if (!file.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }

  if (is.null(fileid)) {
    fileid <- "ag"
  }

  if (type == "png") {
    png_filename <- paste0(output_folder, fileid, ".png")
    if (is.null(width)) {
      width = 1024
    }
    if (is.null(height)) {
      height = 1024
    }
  } else if (type == "pdf") {
    pdf_filename <- paste0(output_folder, fileid, ".pdf")
    if (is.null(width)) {
      width = 5
    }
    if (is.null(height)) {
      height = 5
    }
  }

  dot_filename <- paste0(output_folder, fileid, ".dot")
  svg_filename <- paste0(output_folder, fileid, ".svg")

  graphFile <- file(dot_filename, "w")
  cat('digraph graphname {', file=graphFile)
  cat('node [shape = oval];\n', file=graphFile)

  formatted_labels <- labels
  if (add_index) {
    formatted_labels <- paste0(formatted_labels, "_", 1:length(labels))
  }

  for (v in formatted_labels) {
    cat(v, "[label=", v, "]\n", file=graphFile)
  }

  # For each row:
  for (i in 1:(nrow(amat)-1)) {
    # For each column:
    for (j in (i+1):ncol(amat)) {
      label_i <- formatted_labels[i]
      label_j <- formatted_labels[j]

      if (amat[i,j] > 0) {
        cat(label_i, "->", label_j, "[color=black, dir=both,",  file=graphFile)
        if (amat[i,j] == 1) {
          cat("arrowhead=odot, ", file=graphFile)
        } else if (amat[i,j] == 2) {
          cat("arrowhead=normal, ", file=graphFile)
        } else if (amat[i,j] == 3) {
          cat("arrowhead=none, ", file=graphFile)
        }

        if (amat[j,i] == 1) {
          cat("arrowtail=odot", file=graphFile)
        } else if (amat[j,i] == 2) {
          cat("arrowtail=normal", file=graphFile)
        } else if (amat[j,i] == 3) {
          cat("arrowtail=none", file=graphFile)
        }
        cat("];\n", file=graphFile)
      }
    }
  }
  cat("}\n", file=graphFile)
  close(graphFile)

  DOT::dot(paste(readLines(dot_filename), collapse=" "), file=svg_filename)

  if (type == "png") {
    rsvg::rsvg_png(svg_filename, png_filename, width = 1024, height = 1024)
  } else if (type == "pdf") {
    rsvg::rsvg_pdf(svg_filename, pdf_filename, width = 5, height = 5)
  }
}

igraph_from_graphNel <- function(graphN, latNodes){
  igraph_dag <- igraph::igraph.from.graphNEL(graphN, weight = FALSE)
  for (n in latNodes) {
    adj_list <- graphN@edgeL[[n]]$edges
    if (length(adj_list) == 2) {
      igraph_dag <- igraph::add_edges(igraph_dag, c(adj_list[1], adj_list[2], adj_list[2], adj_list[1]))
      igraph_dag <- igraph::set.edge.attribute(graph = igraph_dag,
                                               name ="description",
                                               index = c(length(igraph::E(igraph_dag))-1, length(igraph::E(igraph_dag))), value = "U")
    }
  }
  for (n in latNodes){
    igraph_dag <- igraph::delete_vertices(igraph_dag, n)
  }
  return(igraph_dag)
}

# returns an pcalg amat (adjacency matrix) of
# type amat.pag (same type for MAGs), where:
# 0: No edge
# 1: Circle
# 2: Arrowhead
# 3: Tail
dagitty2amat <- function(adagg, type="mag") {
  edg <- dagitty:::edges(adagg)
  node_names <- dagitty:::names.dagitty(adagg)
  ans_mat <- matrix(
    data = 0, nrow = length(node_names),
    ncol = length(node_names),
    dimnames = list(node_names, node_names)
  )

  diredg <- subset(edg, e == "->")

  ans_mat[as.matrix(diredg[c("w", "v")])] <- 3
  ans_mat[as.matrix(diredg[c("v", "w")])] <- 2

  bidiredg <-  subset(edg, e == "<->")
  ans_mat[as.matrix(bidiredg[c("w", "v")])] <- 2
  ans_mat[as.matrix(bidiredg[c("v", "w")])] <- 2

  return(ans_mat)
}

