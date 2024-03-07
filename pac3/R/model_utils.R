order_check_quiet = function(tree) {
  # The following is a precondition for the C++ code.
  E = nrow(tree$edge)
  stopifnot("expect even number of edges"=E %% 2 == 0)
  for (e in seq(E, 1, -2)) {
    i1 = tree$edge[e  , 1]
    i2 = tree$edge[e-1, 1]
    stopifnot("expect edges to come in pairs"=i1==i2)
  }
}

run_model = function(tree, alignment, trait_tbl=NULL, transition_matrix, stationary_distro,
                    n_iter=10, inner_iter=500,
                    prior_zr_1=0.5,
                    Phi_prior=c(1.0, 5.0, 1.0, 5.0, 1.0, 1.0),
                    prior_r2=c(5.0, 0.04),
                    prior_r3=c(10.0, 0.2),
                    prior_root=c(0.0, 1.0),
                    prior_lv=c(0.0, 2.0),
                    prior_lb2=c(0.0, 1.0),
                    prior_lb3=c(0.0, 1.0),
                    uniform_prior_lv=FALSE,
                    verbose=FALSE) {
  do_Y = !is.null(trait_tbl)
  if (!do_Y) { warning("no trait supplied") }
  stopifnot(0 <= prior_zr_1 && prior_zr_1 <= 1)
  tree = reorder(tree, "postorder")
  order_check_quiet(tree) # Want to fail immediately in event reorder does not return paired edges in future.
  ape_names = with(tree, c(tip.label, node.label))
  ape_nodes = sort(unique(c(tree$edge)))
  N = length(ape_names)
  L = length(tree$tip.label)
  R = ape_nodes[!(ape_nodes %in% tree$edge[,2])] # the root is not a destination
  E = nrow(tree$edge)
  
  # Alignments are handled by APE.
  alignment = alignment[match(ape_names[1:L], rownames(alignment)),] # put alignment in APE order
  stopifnot(ape_names[1:L]==rownames(alignment))
  convert_1234 = function(aln) {
    aln_ = match(aln, c("a","c","g","t","-"))
    stopifnot(all(complete.cases(aln_)))
    dim(aln_) = dim(aln)
    aln_
  }
  alignment = convert_1234(alignment)
  # Trait codes from data frame with columns "NAME" and "Y". Note, this step pads Y with NA for internal nodes.
  if (do_Y) {
    stopifnot(nrow(trait_tbl)==L)
    stopifnot("Y" %in% colnames(trait_tbl))
    stopifnot(all(trait_tbl$NAME %in% ape_names[1:L]))
    trait <- merge(x=data.frame(ID=ape_nodes, NAME=ape_names), y=trait_tbl, by.x="NAME", by.y="NAME", all.x=T)
    trait <- trait[order(trait$ID),] # ensure rows are in APE order after join
    stopifnot(nrow(trait)==N)
  } else {
    trait <- data.frame(Y=rnorm(N))
  }
  # These changes accommodate 0 based indexing.
  prep_aln = alignment - 1
  prep_aln = rbind(prep_aln, matrix(0, nrow(prep_aln) - 1, ncol(prep_aln))) # pad matrix to include space for latent variables
  prep_edge = as.matrix(tree$edge-1)
  prep_R = R - 1
  # Run C++ code.
  prior_gc_a = Phi_prior[1]
  prior_gc_b = Phi_prior[2]
  prior_lc_a = Phi_prior[3]
  prior_lc_b = Phi_prior[4]
  prior_ac_a = Phi_prior[5]
  prior_ac_b = Phi_prior[6]
  mcmc_trace=Sample_n(n_iter=n_iter, L=L, R=prep_R, edge=prep_edge, distance=tree$edge.length,
                      X=prep_aln, Y=trait$Y, transition_matrix=transition_matrix, pi_distro=stationary_distro,
                      prior_r2_a=prior_r2[1],
                      prior_r2_b=prior_r2[2],
                      prior_r3_a=prior_r3[1],
                      prior_r3_b=prior_r3[2],
                      prior_zr_1=prior_zr_1,
                      prior_gc_a=prior_gc_a,
                      prior_gc_b=prior_gc_b,
                      prior_lc_a=prior_lc_a,
                      prior_lc_b=prior_lc_b,
                      prior_ac_a=prior_ac_a,
                      prior_ac_b=prior_ac_b,
                      prior_root_u =prior_root[1],
                      prior_root_sd=prior_root[2],
                      prior_lv_u=prior_lv[1],
                      prior_lv_sd=prior_lv[2],
                      prior_lb2_u=prior_lb2[1],
                      prior_lb2_sd=prior_lb2[2],
                      prior_lb3_u=prior_lb3[1],
                      prior_lb3_sd=prior_lb3[2],
                      do_Y=do_Y, uniform_prior_lv=uniform_prior_lv, inner_iter=inner_iter, verbose=verbose)
  mcmc_context=list(tree=tree, ape_names=ape_names, ape_nodes=ape_nodes, N=N, L=L, R=R, alignment=alignment)
  colnames(mcmc_trace$z_trace) <- ape_names
  colnames(mcmc_trace$r_trace) <- c("r1","r2","r3")
  colnames(mcmc_trace$Phi_trace) <- c("phi.nc","phi.ca","phi.ac")
  colnames(mcmc_trace$v_trace) <- c("log.v","log.beta2","log.beta3")
  colnames(mcmc_trace$y_trace) <- ape_names
  list(mcmc_trace=mcmc_trace, mcmc_context=mcmc_context)
}

plot_mixing = function(mcmc_out, ith=1) {
  R = mcmc_out$mcmc_context$R
  n_samples = nrow(mcmc_out$mcmc_trace[[1]])
  thin = seq(1, n_samples, by=ith)

  par(mfrow=c(3,3))
#  zs = mcmc_out$mcmc_trace$z_trace[thin,node_index]
#  pz1 = colMeans(zs==0)
#  pz2 = colMeans(zs==1)
#  pz3 = 1 - (pz1 + pz2)
#  plot(pz1, xaxt="n", main="post. Z", yaxt="n", ylab="", xlab="", ylim=c(-1,1), type="n", col="black")
#  points(-pz2, xaxt="n", main="post. Z", yaxt="n", ylab="", xlab="", ylim=c(-1,1), type="h", col="blue",lwd=2)
#  points( pz3, xaxt="n", main="post. Z", yaxt="n", ylab="", xlab="", ylim=c(-1,1), type="h", col="red",lwd=2)
#  axis(side=1, at=1:N, labels = FALSE)
#  plot(colMeans(zs), xaxt="n", main="post. Z", yaxt="n", ylab="", xlab="", ylim=c(0,2), type="p")
#  axis(side=1, at=1:N, labels = FALSE)
#  text(x=1:N, par("usr")[3], labels=ape_names[node_index], srt = 45, pos = 1, xpd = TRUE)
#  axis(side=2, at=0:2, labels = FALSE)
#  text(y=0:2, par("usr")[3], labels = c("neu", "con", "acc"), pos = 2, xpd = TRUE)
#  plot(c(1,1),c(1,1),type="n",xlab="",ylab="")  
  lvb = mcmc_out$mcmc_trace$v_trace[thin,] # also used for 3rd row
  plot(lvb[,3]-lvb[,2], ylab="log b3:b2", type="l", main="log b3:b2")

  
  rs = mcmc_out$mcmc_trace$r_trace[thin,]
  plot(rs[,2], type="l", main="r_2 con")
  plot(rs[,3], type="l", main="r_3 acc")
  
  phis = mcmc_out$mcmc_trace$Phi_trace[thin,]
  plot(phis[,1], type="l",ylim=c(0,1), main="P(neu->con)",ylab="Phi_12 = c")
  plot(phis[,2], type="l",ylim=c(0,1), main="P(con->acc)",ylab="Phi_23 = a")
  plot(phis[,3], type="l",ylim=c(0,1), main="P(acc->con)",ylab="Phi_32 = b")
  
#  ys = mcmc_out$mcmc_trace$y_trace[thin,]
#  plot(ys[,R], type="l", main="Y_root")
#  hist(1:5)
  
  #lvb = mcmc_out$mcmc_trace$v_trace[thin,]
  plot(lvb[,1], type="l", ylab="lv", main="neu")
  plot(lvb[,2], ylab="log b2", type="l", main="con")
  plot(lvb[,3], ylab="log b3", type="l", main="acc")
  par(mfrow=c(1,1))
}

sim_y = function(tree, neu_branches, con_branches, acc_branches, v=1, b2=0.2, b3=2, yR=0) {
  tree = reorder(tree, "postorder")
  ape_names = with(tree, c(tip.label, node.label))
  ape_nodes = sort(unique(c(tree$edge)))
  N = length(ape_names)
  L = length(tree$tip.label)
  R = ape_nodes[!(ape_nodes %in% tree$edge[,2])] # the root is not a destination
  E = nrow(tree$edge)

  Y = rep(NA, N)
  Y[R] = yR
  for (e in seq(E,1,-1)) {
    i = tree$edge[e,1]
    j = tree$edge[e,2]
    nn = ape_names[j]
    if (nn %in% acc_branches) {
      sim_v = v * b3 * tree$edge.length[e]
    } else if (nn %in% con_branches) {
      sim_v = v * b2 * tree$edge.length[e]
    } else {
      stopifnot(nn %in% neu_branches)
      sim_v = v *      tree$edge.length[e]
    }
    Y[j] = rnorm(1, mean=Y[i], sd=sqrt(sim_v));
  }
  trait_tbl=data.frame(NAME=tree$tip.label, Y=Y[1:L])
}

sim_X = function(
  PI,
  Q,
  tree,
  neu_branches, # provide names branches falling into each Z category
  con_branches, # 
  acc_branches,
  con_m=0.1, # relative rates of evolution
  acc_m=0.7,
  S=200 # number of sites in alignment
) {
  tree = reorder(tree, "postorder")
  ape_names = with(tree, c(tip.label, node.label))
  ape_nodes = sort(unique(c(tree$edge)))
  N = length(ape_names)
  L = length(tree$tip.label)
  R = ape_nodes[!(ape_nodes %in% tree$edge[,2])] # the root is not a destination
  E = nrow(tree$edge)
  X = matrix(NA, nrow=N, ncol=S)
  X[R,] = sample(1:4, size=S, replace=T, prob=PI)
  for (e in seq(E,1,-1)) {
    i = tree$edge[e,1]
    j = tree$edge[e,2]
    nn = ape_names[j]
    if (nn %in% acc_branches) {
      rt = tree$edge.length[e] * acc_m
    } else if (nn %in% con_branches) {
      rt = tree$edge.length[e] * con_m
    } else {
      stopifnot(nn %in% neu_branches)
      rt = tree$edge.length[e]
    }
    Qrt = expm::expm(Q*rt)
    #print(Qrt)
    for (s in 1:S) {
      X[j,s] = sample(1:4, size=1, replace=T, prob=Qrt[X[i,s],])
    }
  }
  X_char = matrix(c("a","c","g","t")[X], nrow=N, ncol=S)
  rownames(X_char) = ape_names
  X_char
}

save_fasta = function(X_char, file) {
  as_strings = apply(X_char, 1, paste, collapse="", simplify=T)
  as_fasta_strings = sprintf(">%s\n%s\n", row.names(X_char), as_strings)
  cat(paste(as_fasta_strings, sep="\n"), file=file, append=F)
}

sim_z = function(tree,
                 Z_R_prior=c(0.5,0.5,0), # prior on Z at root of tree
                 P_nc, # to gain conservation
                 P_ca, # to lose conservation
                 P_ac, # to regain conservation
                 plot_it=F)
{
  tree = reorder(tree, "postorder")
  ape_names = with(tree, c(tip.label, node.label))
  ape_nodes = sort(unique(c(tree$edge)))
  N = length(ape_names)
  L = length(tree$tip.label)
  R = ape_nodes[!(ape_nodes %in% tree$edge[,2])] # the root is not a destination
  E = nrow(tree$edge)
  
  Z = rep(NA, N) # store rate category here
  Z[R] = sample(1:3, size=1, prob=Z_R_prior) # N.B. treating Z in {1,2,3} in R code
  p_change = c(P_nc, P_ca, P_ac) # a.k.a. Phi entries
  # evolve everything together
  for (e in seq(E,1,-1)) {
    i = tree$edge[e,1]
    j = tree$edge[e,2]
    stopifnot(j!=R) # Z_R is set before loop
    stopifnot(!is.na(Z[i])) # from edge must be set
    stopifnot( is.na(Z[j])) # to edge must not be set
    if (Z[i] < 3) {
      Z[j] = Z[i] + rbinom(1, 1, prob=p_change[Z[i]])
    } else {
      Z[j] = Z[i] - rbinom(1, 1, prob=p_change[Z[i]]) # acc -> con is 3 -> 2
    }
    stopifnot(Z[j] %in% 1:3)
  }
  if (plot_it) {
    plot(tree, edge.color=c("black","blue","red")[Z[tree$edge[,2]]], edge.width=2, show.node.label=F)
  }
  names(Z) = ape_names
  Z
}

plot_z = function(mcmc_out, drop_n=1000, thin_n=1) {
  R = mcmc_out$mcmc_context$R
  N = mcmc_out$mcmc_context$N
  n_samples = nrow(mcmc_out$mcmc_trace[[1]])
  thin = seq(1, n_samples, by=thin_n)
  acc_prob = colMeans(2 == mcmc_out$mcmc_trace$z_trace[-(1:drop_n),])
  con_prob = colMeans(1 == mcmc_out$mcmc_trace$z_trace[-(1:drop_n),])
  neu_prob = 1 - acc_prob - con_prob
  edge_color = rgb(acc_prob, 0, con_prob)[mcmc_out$mcmc_context$tree$edge[,2]]
  plot(mcmc_out$mcmc_context$tree, show.tip.label=T, use.edge.length=T, edge.color=edge_color, main="posterior Z", edge.width=2)
}

plot_v = function(mcmc_out, drop_n=1000, thin_n=1) {
  n_samples = nrow(mcmc_out$mcmc_trace[[1]])
  thin = seq(1, n_samples, by=thin_n)
  log_b2 = mcmc_out$mcmc_trace$v_trace[-(1:drop_n),2]
  log_b3 = mcmc_out$mcmc_trace$v_trace[-(1:drop_n),3]
  hist(log_b3-log_b2, main="posterior log b3:b2")
  abline(v=mean(log_b3-log_b2), lty=2)
}

#plot_v = function(mcmc_out, drop_n=1000, thin_n=1) {
#  R = mcmc_out$mcmc_context$R
#  N = mcmc_out$mcmc_context$N
#  n_samples = nrow(mcmc_out$mcmc_trace[[1]])
#  thin = seq(1, n_samples, by=thin_n)
#  log_b23 = mcmc_out$mcmc_trace$v_trace[-(1:drop_n),2:3]
#  sm::sm.density(log_b23, model="Normal", display="slice", xlab="log[b2]", ylab="log[b3]", xlim=c(-3,3), ylim=c(-3,3), props=c(10,50,90), lwd=2)
#  abline(h=0, lty=3)
#  abline(v=0, lty=3)
#  null_data = MASS::mvrnorm(n=9000, mu=c(0,0), Sigma=matrix(c(1,0,0,1), ncol=2))
#  sm::sm.density(null_data, display="slice", model="Normal", lty=2, add=T,  props=c(10,50,90))
#  text(-2.5, 2.5, "establishing")
#  text( 2.5,-2.5, "enabling")
#  title(main="posterior beta")
#}

#plot_r = function(mcmc_out) {
#	stopifnot(FALSE)
#}
