library(gplots)

p.R.z = function(R, z, prior.beta) {
  # Likelhood of data (R) given cluster assignment z. Evaluated by integrating
  # out possible values of the Bernoulli probabilities of the relations (which
  # is governed by prior.beta)
  z.uniq = unique(z)
  z.len = length(z.uniq)
  m.ab = matrix(0, nrow = z.len, ncol = z.len)
  m.ab.bar = matrix(0, nrow = z.len, ncol = z.len)
  # Loop through matrix. length(z.uniq))
  for (i in 1:nrow(R)) {
    for (j in 1:ncol(R)) {
      a = z[i]; b = z[j]
      R.ij = R[i, j]
      m.ab[a, b] = m.ab[a, b] + R.ij
      m.ab.bar[a, b] = m.ab.bar[a, b] + (1 - R.ij)
    }
  }
  ll = 0
  for (a in 1:nrow(m.ab)) {
    for (b in 1:nrow(m.ab)) {
      ll.ab = beta(m.ab[a, b] + prior.beta, m.ab.bar[a, b] + prior.beta) /
        beta(prior.beta, prior.beta)
      ll = ll + log(ll.ab)
    }
  }
  ll
}

p.z = function(z, prior.gamma) {
  # Likelihood of this z being drawn from a CRP given GAMMA.
  # TODO: Sorting is probably really slow...
  z = sort(z)

  ll = 0
  # Keeping track of members already part of the clusters.
  # Requires that we know how many unique clusters there are
  n = rep(0, length(unique(z)))
  stopifnot(z[1] == 1)
  for (i in 1:length(z)) {
    a = z[i]
    n.a = n[a]
    if (n.a == 0) { # New cluster
      ll = ll + log(prior.gamma / (i - 1 + prior.gamma))
    } else { # Existing cluster
      ll = ll + log(n.a / (i - 1 + prior.gamma))
    }
    # Regardless, increment n.a
    n[a] = n[a] + 1
  }
  ll
}

renumber = function(arr) {
  # First, collapse extras
  arr = match(arr, unique(sort(arr)))
  # Second, how many clusters are there now?
  n.clus = length(unique(arr))
  a.map = rep(NA, n.clus)
  # Now loop from beinning to end of arr. When seeing a new number
  next.a = 1
  for (i in 1:length(arr)) {
    a = arr[i]
    if (is.na(a)) {
      next
    }
    if (is.na(a.map[a])) {
      a.map[a] = next.a
      arr[i] = next.a
      next.a = next.a + 1
    } else {
      arr[i] = a.map[a]
    }
  }
  arr
}

samp.z.probs = function(obs, prior.gamma) {
  # The probabilities of cluster assignments given already observed values.
  # The ith value is the length of obs + 1
  i = length(obs) + 1
  n.clus = length(unique(obs))
  n = rep(0, n.clus)
  for (j in 1:length(obs)) {
    a = obs[j]
    n[a] = n[a] + 1
  }
  # Now consider sampling i into any of the known as, OR a new a.
  p = rep(0, n.clus + 1)
  for (a in 1:n.clus) {
    # What is the probability that i belongs to a?
    n.a = n[a]
    stopifnot(n.a != 0)
    p[a] = (n.a / (i - 1 + prior.gamma))
  }
  # Sample from n.clus + 1 because you could sample a new cluster.
  p[n.clus + 1] = 1 - sum(p)
  p
}

samp.z = function(obs, prior.gamma) {
  # Samples from samp.z.probs.
  p = samp.z.probs(obs, prior.gamma)
  sample(length(p), 1, prob = p)
}

j = function(z.i, i, obs, prior.gamma) {
  # What was the loglik of sampling z.i as the ith position of obs?
  probs = samp.z.probs(renumber(obs[-i]), prior.gamma)
  log(probs[z.i])
}

# Simple 2-dimensional IRM. For each assignment "sample", I sweep through the
# cluster assignments of each individual point, and use a metropolis-hastings
# update.
irm = function(R, sweeps = 1000, prior.beta = 1, prior.gamma = 1) {

  Z = matrix(NA, nrow = sweeps, ncol = nrow(R))

  z = 1:nrow(R)
  n = length(z)

  # Sweeps.
  for (s in 1:sweeps) {
    for (i in 1:n) {
      # Propose changing this value by sampling from the existing CRP.
      # Specifically, since this is exchangeable, we imagine that z.i is the last
      # value to be sampled, and we sample from the probabilities.
      z.i = z[i]
      # Setup new vector. ith value to be filled in
      z.star = rep(z)
      z.star[i] = NA
      # Renumber (leaves NAs alone)
      z.star = renumber(z.star)

      # Get vector without the ith value and sample a new z.i
      z.i.star = samp.z(z.star[-i], prior.gamma)
      z.star[i] = z.i.star

      # Now compare likelihoods of 1) original z, and 2) new z.samp
      # Initial likelihood of original z and new z.samp
      log.ll = (p.R.z(R, z.star, prior.beta) +
                  p.z(z.star, prior.gamma)) -
        (p.R.z(R, z, prior.beta) + p.z(z, prior.gamma))
      # Next, relative likelihood of proposal distribution
      # likelihood of sampling z.i given z.i.star vs likelihood of sampling
      # z.i.star given z.i
      log.prop = j(z.i, i, z.star, prior.gamma) - j(z.i.star, i, z, prior.gamma)

      # Final acceptance ratio
      log.r = log.ll + log.prop

      if (log(runif(1)) < log.r) {
        z = renumber(z.star)
      }
    }
    # We ignore storing results of sweeps
    Z[s, ] = z
  }

  Z
}

# Collapse Z to strings
top.n = function(Z, n = 10) {
  Z.str = rep("", nrow(Z))
  for (i in 1:nrow(Z)) {
    Z.str[i] = paste(Z[i, ], collapse = "")
  }

  Z.counts = table(Z.str)
  # Only keep top 10 counts
  Z.top = head(sort(Z.counts, decreasing = TRUE), n = n)
  c(Z.top)
}

mode.irm = function(Z) {
  top.1 = names(top.n(Z, n = 1))
  as.numeric(strsplit(top.1, "")[[1]])
}

plot.R = function(R, assn = NULL) {
  if (is.null(assn)) {
    heatmap.2(R,
              trace = 'none',
              dendrogram = 'none',
              Colv = 'none', Rowv = 'none',
              key = FALSE,
              col = c('#eeeeee', '#000000')
              )
  } else {
    assn.sort = sort(assn)
    # Where does assn change?
    assn.indices = c(1, 1 + which(diff(assn.sort) != 0)) - 1
    R.assn = R[order(assn), order(assn)]
    heatmap.2(R.assn,
              trace = 'none',
              dendrogram = 'none',
              Colv = 'none', Rowv = 'none',
              colsep = assn.indices,
              rowsep = assn.indices,
              key = FALSE,
              col = c('#eeeeee', '#000000')
              )
  }
}
