##===== Prepare workspace =====#

## Load required libraries; install if not found
load.libs <- function() {

  libs <- c('ggplot2', 'ggforce', 'dplyr', 'reshape2','RColorBrewer', 'scales',
            'gridExtra', 'grid', 'magrittr', 'outbreaker', 'outbreaker2',
            'EpiEstim')
  
  for(i in libs) {
    if(!require(i, character.only = TRUE)) {
      install.packages(i)
    }
  }

  if(!require('distcrete')) {
    devtools::install_github('reconhub/distcrete')
  }
  
}


##===== Main analysis functions =====#

## Simulate outbreaks and reconstruct them using T, TC, TG and TCG
analyse <- function(param, w.dens, f.dens, mu.transi, seq.length, offsp) {

  ## Defining global properties of run
  n.hosts <- 200
  min.cases <- 60
  max.cases <- min.cases + 1
  rate.import.case <- 0
  prior.eps <- c(0,1)
  R0 = seq(0, 20, 0.1)

  runs <- param$runs
  eps <- param$eps
  lambda <- psi.to.lambda(1, param$psi, min.cases, 1)
  
  analysis <- data.frame(variable = rep(c("accuracy", "entropy"), 4*runs),
                         model = rep(c("t", "tc", "tg", "tcg"), runs, each = 2),
                         value = NA)

  store.outbreak <- store.CTD <- store.t.result <- store.tc.result <-
    store.tg.result <- store.tcg.result <- list()

  ind  <- 0

  while((ind + 1) <= runs) {

    ## Simulate outbreak
    outbreak <- outbreaker::simOutbreak(R0 = R0,
                                        infec.curve = w.dens,
                                        n.hosts = n.hosts,
                                        duration = 60,
                                        mu.transi = mu.transi,
                                        group.freq = offsp,
                                        rate.import.case = rate.import.case,
                                        seq.length = seq.length)
    if(outbreak$n < min.cases) next
    if(outbreak$n > min.cases) {
      outbreak <- subset.outbreak(outbreak, 1:min.cases)
    }

    ## Draw sampling time from incubation period
    sample_delay  <- sample(0:(length(f.dens) - 1), outbreak$n,
                            prob = f.dens, replace = TRUE)
    outbreak$onset <- outbreak$onset + sample_delay

    ## Simulate contact data
    tTree <- data.frame(i = outbreak$ances, j = outbreak$id)
    CTD <- sim_ctd(tTree, eps = eps, lambda = lambda)[,1:2]
    if(nrow(CTD) == 0) CTD <- NULL

    ## Reconstruct using T
    data <- list(dates = outbreak$onset, w_dens = w.dens, f_dens = f.dens)
    config <- list(n_iter = 10000, sample_every = 50, find_import = FALSE,
                   move_kappa = FALSE, move_pi = FALSE, init_pi = 1,
                   init_tree = 'star')

    t.result <- outbreaker2::outbreaker(data = data, config = config)
    analysis$value[8*ind + 1]  <- get.acc(t.result, outbreak)
    analysis$value[8*ind + 2]  <- get.ent(t.result)

    ## Reconstruct using TC
    data$ctd <- CTD

    tc.result <- outbreaker2::outbreaker(data = data, config = config)
    analysis$value[8*ind + 3]  <- get.acc(tc.result, outbreak)
    analysis$value[8*ind + 4]  <- get.ent(tc.result)

    ## Reconstruct using TG
    data$ctd <- NULL
    data$dna <- outbreak$dna
    config$init_tree  <- 'seqTrack'

    tg.result <- outbreaker2::outbreaker(data = data, config = config)
    analysis$value[8*ind + 5]  <- get.acc(tg.result, outbreak)
    analysis$value[8*ind + 6]  <- get.ent(tg.result)

    ## Reconstruct using TCG
    data$ctd <- CTD

    tcg.result <- outbreaker2::outbreaker(data = data, config = config)
    analysis$value[8*ind + 7]  <- get.acc(tcg.result, outbreak)
    analysis$value[8*ind + 8]  <- get.ent(tcg.result)

    ## Remove genetic data to free up space
    outbreak$dna <- NULL

    ## store results
    store.outbreak[[ind + 1]] <- outbreak
    store.CTD[[ind + 1]] <- CTD
    store.t.result[[ind + 1]] <- t.result
    store.tc.result[[ind + 1]] <- tc.result
    store.tg.result[[ind + 1]] <- tg.result
    store.tcg.result[[ind + 1]] <- tcg.result

    ind <- ind + 1

  }

  param <- list(n.hosts = n.hosts, min.cases = min.cases, w.dens = w.dens,
                n.hosts = n.hosts, max.cases = max.cases, eps = eps,
                lambda = lambda, runs = runs)

  out <- list()
  out$param <- param
  out$outbreak <- store.outbreak
  out$CTD <- store.CTD
  out$t.result <- store.t.result
  out$tc.result <- store.tc.result
  out$tg.result <- store.tg.result
  out$tcg.result <- store.tcg.result
  out$analysis <- analysis

  return(out)
}

## Run analysis across SARS- and Ebola-like outbreaks
run.analysis <- function(param) {

  ## Mutation rates are scaled by two thirds, as mu.transv is added by default
  
  ebola <- analyse(param = param,
                   w.dens = discr.gamma(mean = 15.3, sd = 9.3),
                   f.dens = discr.gamma(mean = 9.1, sd = 7.3),
                   offsp = create.offsp(R0 = 1.8, k = 0.18),
                   mu.transi = 1.24e-3/365*(2/3),
                   seq.length = 18958)

  sars <- analyse(param = param,
                  w.dens = discr.gamma(mean = 8.7, sd = 3.6),
                  f.dens = discr.gamma(mean = 6.37, sd = 4.09),
                  offsp = create.offsp(R0 = 3.5, k = 0.16),
                  mu.transi = 1.14e-5*(2/3),
                  seq.length = 29750)

  out <- list(ebola = ebola, sars = sars)

  return(out)

}

## Extract summary results from raw results
create.store <- function(dir) {
  
  store = list()
  
  create.est <- function(disease,temp.runs) {

    calc.est <- function(model, var){
      sapply(seq_len(temp.runs),
             function(i) mean(r[[disease]][[paste0(model,".result")]][[i]][[var]][22:201]))
    }

    data.frame(Model=rep(c("tcg","tc","tg","t"),2,each=temp.runs),
               Variable=rep(c("eps","lambda"),each=4*temp.runs),
               Value=c(calc.est("tcg","eps"),
                       calc.est("tc","eps"),
                       calc.est("tg","eps"),
                       calc.est("t","eps"),
                       calc.est("tcg","lambda"),
                       calc.est("tc","lambda"),
                       calc.est("tg","lambda"),
                       calc.est("t","lambda")))

  }

  create.offsp <- function(r, disease) {
    as.vector(sapply(seq_len(r[[disease]]$param$runs),
                     function(i) get.offsp(r[[disease]]$outbreak[[i]])))
  }

  add.r <- function(store,r) {

    ## Adding Ebola
    temp.runs <- r$ebola$param$runs
    
    if(is.null(store[[run]]$ebola$param)) {
      store[[run]]$ebola$param  <- r$ebola$param
    } else {
      store[[run]]$ebola$param$runs  <- store[[run]]$ebola$param$runs + temp.runs
    }

    store[[run]]$ebola$analysis <- rbind(store[[run]]$ebola$analysis,
                                         r$ebola$analysis)

    ## Adding est; the mean parameter estimate of the posterior distribution
    store[[run]]$ebola$est <- rbind(store[[run]]$ebola$est,
                                    create.est("ebola",temp.runs=temp.runs))

    ## Adding the offspring distribution
    store[[run]]$ebola$offsp <- c(store[[run]]$ebola$offsp, create.offsp(r, 'ebola'))


    ## Adding SARS
    temp.runs <- r$sars$param$runs
    
    if(is.null(store[[run]]$sars$param)) {
      store[[run]]$sars$param  <- r$sars$param
    } else {
      store[[run]]$sars$param$runs  <- store[[run]]$sars$param$runs + temp.runs
    }
    
    store[[run]]$sars$analysis <- rbind(store[[run]]$sars$analysis,
                                        r$sars$analysis)

    ## Adding est
    store[[run]]$sars$est <- rbind(store[[run]]$sars$est,
                                   create.est("sars",temp.runs=temp.runs))

    ## Adding the offspring distribution
    store[[run]]$sars$offsp <- c(store[[run]]$sars$offsp, create.offsp(r, 'sars'))

    return(store)
    
  }

  ## Load previously downloaded results into store
  files <- list.files(dir)
  files <- files[-grep('store', files)]
  if(length(files) == 0) files <- list.files(dir)
  pb <- txtProgressBar(min = 1, max = length(files),style = 3)
  for(file in files){
    setTxtProgressBar(pb, which(files == file))
    load(paste0(dir, file))
    run <- paste0("eps", r$ebola$param$eps, "lambda", r$ebola$param$lambda)
    store <- add.r(store, r)
  }

  return(store)
}


##===== Plotting functions ======##

## Plots the accuracy of transmission tree inference
vis.acc <- function(store, mod = c("t", "tc", "tg", "tcg")) {

  psi <- lambda.to.psi(1, 1, store[[1]]$ebola$param$min.cases, 1)

  df <- mk.dwide(store, mod, 'accuracy')

  y_axs <- 1:length(store) %>%
    sapply(function(i) store[[i]]$ebola$param$lambda) %>%
    unique %>%
    sort %>%
    "*"(psi) %>%
    round(1)

  x_axs <- sort(unique(sapply(1:length(store), function(i) store[[i]]$ebola$param$eps)))

  mk.lab <- function(variable) {

    tmp <- c(sars = 'SARS-CoV', ebola = 'EBOV')
    tmp[variable]

  }


  p <- ggplot(df, aes(x = Eps, y = Lambda, fill = Accuracy)) +
    geom_tile(color = "white") +
    facet_grid(Disease ~ Data, labeller = labeller(Disease = mk.lab)) +
    scale_fill_gradientn(colours=c("white","darkgreen"),
                         limits=c(min(df$Accuracy),1)) +
    theme_gray(base_size = 10) +
    xlab(expression(paste("Contact reporting coverage (",epsilon,")"))) +
    ylab(expression(paste("Number of non-infectious contacts per person (",psi,")"))) +
    scale_x_continuous(expand=c(0,0), breaks = x_axs) +
    scale_y_continuous(expand=c(0,0), breaks = y_axs) +
    theme(axis.ticks=element_blank(),
          panel.border = element_rect(fill=NA,colour="black",size=1)) +
    guides(fill = guide_colorbar(title.vjust = 0,title.hjust=0.5,reverse=TRUE,title="Accuracy"))

  p
}

## Plot change in accuracy (absolute) between dat1 and dat2
vis.rel <- function(store, dat1 = 'tc', dat2 = 'tg', var = 'accuracy') {

  psi <- lambda.to.psi(1, 1, store[[1]]$ebola$param$min.cases, 1)

  dwide <- mk.dwide(store, var = var)
  dwide$Disease <- mk.name(dwide$Disease)

  df1 <- filter(dwide, Data == toupper(dat1))
  df2 <- filter(dwide, Data == toupper(dat2))

  df1$Accuracy <- df1$Accuracy - df2$Accuracy

  if(var == 'entropy') {
    grad <- scale_fill_gradient2(low="darkgreen",mid="white",high="darkred")
  } else {
    grad <- scale_fill_gradient2(low="darkred",mid="white",high="darkgreen")
  }

  p <- ggplot(df1, aes(x = Eps, y = Lambda, fill = Accuracy)) +
    geom_tile(color = "white") + facet_grid(Disease ~ ., labeller = label_parsed) +
    grad +
                                        #scale_fill_gradientn(colours=c("white","darkgreen"),limits=c(min(dwide$Accuracy),1)) +
    theme_gray(base_size = 10) +
    xlab(expression(paste("Contact reporting coverage (",epsilon,")"))) +
    ylab(expression(paste("Number of non-infectious contacts per person (",psi,")"))) +
                                        #ylab(expression(paste("Number of non-infectious contacts (",lambda,")"))) +
    scale_x_continuous(expand=c(0,0),
                       breaks = round(sort(unique(df1$Eps)), 2)) +
    scale_y_continuous(expand=c(0,0),
                       breaks = round(sort(unique(df1$Lambda)), 1)) +
    theme(axis.ticks   = element_blank(),
          panel.border = element_rect(fill=NA,colour="black",size=1)) +
    guides(fill = guide_colorbar(title.vjust = 0,
                                 title.hjust=0.5,
                                 reverse=TRUE,
                                 title=paste0("Change in\n", var)))

  p

}

## Plots the posterior entropy of ancestry assignments
vis.ent <- function(store, mod = c("t", "tc", "tg", "tcg")) {

  psi <- lambda.to.psi(1, 1, store[[1]]$ebola$param$min.cases, 1)

  df <- mk.dwide(store, mod, var = 'entropy')

  y_axs <- 1:length(store) %>%
    sapply(function(i) store[[i]]$ebola$param$lambda) %>%
    unique %>%
    sort %>%
    "*"(psi) %>%
    round(1)

  x_axs <- sort(unique(sapply(1:length(store), function(i) store[[i]]$ebola$param$eps)))

  mk.lab <- function(variable) {

    tmp <- c(sars = 'SARS-CoV', ebola = 'EBOV')
    tmp[variable]

  }

  p <- ggplot(df, aes(x = Eps, y = Lambda, fill = Accuracy)) +
    geom_tile(color = "white") +
    facet_grid(Disease ~ Data, labeller = labeller(Disease = mk.lab)) +
    scale_fill_gradientn(colours=c("darkgreen","white"),
                         limits=c(min(df$Accuracy),max(df$Accuracy))) +
    theme_gray(base_size = 10) +
    xlab(expression(paste("Contact reporting coverage (",epsilon,")"))) +
    ylab(expression(paste("Number of non-infectious contacts per person (",psi,")"))) +
    scale_x_continuous(expand=c(0,0), breaks = x_axs) +
    scale_y_continuous(expand=c(0,0), breaks = y_axs) +
    theme(axis.ticks=element_blank(),
          panel.border = element_rect(fill=NA,colour="black",size=1)) +
    guides(fill = guide_colorbar(title.vjust = 0,title.hjust=0.5,reverse=TRUE,title="Entropy"))

  p
}

## Plot the parameter estimates for eps and lambda
vis.est <- function(store, dis = c("Ebola", "SARS")) {
  
  df <- NULL

  for(i in seq_along(store)) {
    for(disease in c("ebola","sars")) {
      for(mod in c("tc","tcg")) {
        for(var in c("eps","lambda")){
          temp <- subset(store[[i]][[disease]]$est, Model == mod & Variable == var)
          temp$Eps <- store[[i]][[disease]]$param$eps
          temp$Lambda <- store[[i]][[disease]]$param$lambda
          if(disease == "ebola") temp$Disease <- "EBOV"
          if(disease == "sars") temp$Disease <- "SARS"
          df <- rbind(df, temp)
        }
      }
    }
  }

  df %<>% filter(Eps != 0)

  df$Model <- toupper(df$Model)
  names(df)[which(names(df) == 'Model')] <- 'Data'

  df1 <- subset(df, Variable == "eps" & Disease %in% dis)
  df2 <- subset(df, Variable == "lambda" & Disease %in% dis)

  df1 <- df1[rev(order(df1$Lambda)),]
  df2 <- df2[rev(order(df2$Lambda)),]

  df1$lab.Eps <- factor(paste0("epsilon==",df1$Eps))
  df1$lab.Lambda <- factor(paste0("lambda==", round(df1$Lambda, 2)))

  df2$lab.Eps <- factor(paste0("epsilon==",df2$Eps))
  df2$lab.Lambda <- factor(paste0("lambda==", round(df2$Lambda, 2)))

  df1$lab.Lambda <- factor(df1$lab.Lambda, levels = rev(levels(df1$lab.Lambda)))
  df2$lab.Lambda <- factor(df2$lab.Lambda, levels = rev(levels(df2$lab.Lambda)))

  create.histgrid <- function(df, var, main, lab.x, lim, brk) {

    p <- ggplot(df, aes(Value)) + geom_density(aes(y = ..scaled.., fill = Data),
                                               alpha = 0.5, size = 0.2, adjust = 1.5) +
      scale_x_continuous(breaks = brk, limits = lim) +
      facet_grid(lab.Lambda ~ lab.Eps, labeller = label_parsed, scales = 'free') +
      scale_fill_manual(values = match.col("data")) +
      ggtitle(main) + xlab(lab.x) + ylab("Density") +
      theme_gray(base_size = 10) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.text.y      = element_blank(),
            axis.ticks.y     = element_blank(),
            axis.text.x      = element_text(size = 4),
            strip.text       = element_text(size = 6),
            plot.title       = element_text(size = 10),
            legend.key.size  = unit(0.5, 'line'),
            legend.position  = "bottom")

    dat <- ggplot_build(p)
    yend <- max(dat$data[[1]]$density)
    p <- p + geom_segment(aes_string(x = var, y = 0, xend = var, yend = 1),
                          linetype = "dotted", size = 0.5)

    return(p)
  }

  p1 <- create.histgrid(df1, "Eps",
                        expression(paste("Estimating the contact reporting coverage (", epsilon, ")")),
                        expression(paste(epsilon)),
                        lim = c(0, 1), brk = c(0, 0.5, 1))
  p2 <- create.histgrid(df2, "Lambda",
                        expression(paste("Estimating the non-infectious contact probability (", lambda, ")")),
                        expression(paste(lambda)),
                        lim = c(0, 0.6), brk = c(0, 0.3, 0.3))

  p <- joint.legend(p1, p2, 'bottom')

  return(p)
}


##===== Plot saving functions =====##

## Save plots to a given directory
fig.save <- function(p, dir, name, ext = 'svg', ...) {

  ggsave(p, file=paste0(dir, name, ".", ext), ...)
  return(NULL)
  
}

## Create the manuscript figures
create.figs <- function(store, dir) {

  p <- vis.acc(store, mod = c("t", "tc", "tg", "tcg"))
  fig.save(p, dir, "Fig1", width = 7.5, height = 3.9, dpi = 600, ext = 'tiff')

  p <- vis.rel(store)
  fig.save(p, dir, "Fig2", width = 3.9, height = 3.9, dpi = 600, ext = 'tiff')

  p <- vis.ent(store, mod = c("t", "tc", "tg", "tcg"))
  fig.save(p, dir, "S1Fig", width = 7.5, height = 3.9, dpi = 600, ext = 'tiff')

  p <- vis.est(store, 'EBOV')
  fig.save(p, dir, "S2Fig", width = 7.5, height = 3.9, dpi = 600, ext = 'tiff')

  p <- vis.est(store, 'SARS')
  fig.save(p, dir, "S3Fig", width = 7.5, height = 3.9, dpi = 600, ext = 'tiff')

}


##===== Auxilliary functions =====##

## Simulate contact data from a transmission tree using eps and lambda
sim_ctd <- function(tTree, eps, lambda) {

  if(any(c(eps, lambda) < 0) | any(c(eps, lambda) > 1)) {
    stop('eps and lambda must be probabilities')
  }

  id <- unique(c(tTree[,1], tTree[,2]))
  id <- id[!is.na(id)]

  ## Sort tTree by value or alphabetically, This ensures A:B and B:A are both
  ## recognised when querying the contacts dataframe for transmission pairs
  tTree <- tTree %>%
    stats::na.omit() %>%
    apply(1, sort, FALSE) %>%
    t() %>%
    as.data.frame(stringsAsFactors = FALSE)

  tTree <- tTree[order(tTree[1]),]
  names(tTree) <- c('V1', 'V2')
  if(nrow(tTree) == 0) stop("No transmission observed")

  ## Create a dataframe of all potential contacts
  contacts <- as.data.frame(t(utils::combn(id, 2)))

  ## Create a column of logicals indicating whether a pair represents a
  ## transmission pair or not
  tTree$tPair <- TRUE

  ## Mark non-transmission pairs in the contacts dataframe. The merge function
  ## will mark pairs found in contacts but not in tTree as 'NA'. These are
  ## then converted to FALSE
  contacts <- merge(contacts, tTree, by = c('V1', 'V2'), all.x = TRUE)
  contacts$tPair[is.na(contacts$tPair)] <- FALSE

  ## Sample a number of rows given by a binomial distribution
  sampler <- function(x, prob) {
    x[sample(1:nrow(x), stats::rbinom(1, nrow(x), prob)), 1:3]
  }

  ## Sample transmission pairs with probability eps
  ## Sample non-transmission pairs with probability eps*lambda
  ctd <- rbind(sampler(contacts[contacts$tPair,], eps),
               sampler(contacts[!contacts$tPair,], eps*lambda))

  rownames(ctd) <- NULL

  return(ctd)
}

## Returns a discretized gamma distribution
discr.gamma <- function(mean, sd) {
  w <- sapply(0:100, EpiEstim::DiscrSI, mean, sd)
  return(w)
}

## Creates an offspring distribution from mean R0 and k
create.offsp <- function(R0, k) {

  ## v is the distribution of individual reproduction numbers
  v <- distcrete::distcrete('gamma', 1, shape = k, rate = k/R0)
  v$d(x = seq(0, 10, 1))

}

## Returns the accuracy of transmission tree inference
get.acc <- function(res, sim, burnin = 1000) {

  inferred <- summary(res, burnin = burnin)$tree$from
  true <- sim$ances

  inferred[is.na(inferred)] <- 0
  true[is.na(true)] <- 0

  acc <- mean(inferred == true)

  return(acc)

}

## Returns the mean entropy of ancestry assignments
get.ent <- function(res) {

  i <- grep('alpha', names(res))

  calc.ent <- function(x) {
    x <- as.character(x)
    fk <- table(x)/sum(table(x))
    return(-sum(log(fk)*fk))
  }

  ent <- sapply(res[i], calc.ent)

  return(round(mean(ent), 2))

}

## Returns the offspring distribution of a simOutbreak objects
get.offsp <- function(sim) {

  offsp <- table(sim$ances)

  ## add zeroes which aren't listed by table
  offsp <- as.vector(c(rep(0, sim$n - length(offsp)), offsp))

  return(offsp)

}

## Calculate lambda from psi
psi.to.lambda <- function(eps, psi, n, I) {
  psi/((n-3 + 2*I/n)*eps)
}

## Calculates the theoretical psi from CTD parameters
lambda.to.psi <- function(eps, lambda, n, I) {
  (n-3 + 2*I/n)*eps*lambda
}

## Subset an outbreak object to keep only the first n.cases
subset.outbreak <- function(outbreak, to.keep) {

  outbreak$n <- length(to.keep)
  outbreak$dna <- outbreak$dna[to.keep,]
  outbreak$onset <- outbreak$onset[to.keep]
  outbreak$id <- outbreak$id[to.keep]
  outbreak$ances <- outbreak$ances[to.keep]
  outbreak$group <- outbreak$group[to.keep]
  outbreak$nmut <- outbreak$nmut[to.keep]
  outbreak$ngen <- outbreak$ngen[to.keep]

  return(outbreak)
}

## Make a dataframe with mean accuracy values per grid point
mk.dwide <- function(store, mod = c("t", "tc", "tg", "tcg"), var = 'accuracy') {

  psi <- lambda.to.psi(1, 1, store[[1]]$ebola$param$min.cases, 1)

  dwide <- data.frame(matrix(nrow=length(store)*2*4,ncol=5,
                             dimnames=list(c(),c("Disease","Data","Eps","Lambda","Accuracy"))))
  counter <- 1

  for(run in names(store)){
    for(disease in c("ebola","sars")){
      for(dat in c("t","tc","tg","tcg")){
        dwide$Disease[counter] <- disease
        dwide$Data[counter] <- toupper(dat)
        dwide$Eps[counter] <- store[[run]][[disease]]$param$eps
        dwide$Lambda[counter] <- store[[run]][[disease]]$param$lambda*psi
        dwide$Accuracy[counter] <- mean(subset(store[[run]][[disease]]$analysis,
                                               variable==var & model==dat)$value)
        counter <- counter + 1
      }
    }
  }

  df <- subset(dwide, Data %in% toupper(mod))
  df$Data <- factor(df$Data, levels = toupper(mod))

  return(df)

}

## Extracts a ggplot legend
extract.legend <- function(gplot) {

  tmp <- ggplot_gtable(ggplot_build(gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)

}

## Arranges two plots horizontally with a join legend at the bottom (from p1)
joint.legend <- function(p1, p2, pos = 'bottom') {

  leg <- extract.legend(p1)

  if(pos == 'bottom') {
    p <- grid.arrange(arrangeGrob(p1 + theme(legend.position = 'none'),
                                  p2 + theme(legend.position = 'none'),
                                  ncol = 2),
                      leg, nrow = 2, ncol = 1, heights = c(10, 0.5))
  }

  if(pos == 'right') {
    p <- grid.arrange(arrangeGrob(p1 + theme(legend.position = 'none'),
                                  p2 + theme(legend.position = 'none'),
                                  nrow = 2),
                      leg, nrow = 1, ncol = 2, widths = c(10, 0.5))
  }

  return(p)

}

## Returns a named vector for colour/fill_manual; for consistent colour scheme
match.col <- function(var) {

  if(var == "data") {
    cols = rep(brewer.pal(4,"Set1"), 2)
    dats <- c("t", "tc", "tg", "tcg")
    names(cols) <- c(dats, toupper(dats))
    return(cols)
  }

  if(var == "disease") {
    pal <- wesanderson::wes_palette("Royal1",4)[c(3,2)]
    cols <- rep(pal, 2)
    names(cols) <- c("Ebola", "SARS", "ebola", "sars")
    return(cols)
  }
}

## Convert 'ebola' and 'sars' to 'EBOV' and 'SARS-CoV' respectively
mk.name <- function(x) {

  x[x == 'ebola'] <- 'EBOV'
  x[x == 'sars'] <- 'SARS-CoV'

  x

}
