# show distinct color for positive and negative value
plot.matrix <- function (phi, p, name = NULL) {
  B <- phi
  if (nrow(B) == 1) {
    B <- matrix(B[, 1:ncol(B)], nrow = 1)
  }
  else {
    B <- B[, 1:ncol(B)]
  }
  k <- nrow(B)
  s1 <- 0
  m <- 0
  s <- 0
  s <- s + s1
  text <- c()
  for (i in 1:p) {
    text1 <- as.expression(bquote(bold(Phi)^(.(i))))
    text <- append(text, text1)
  }
  if (m > 0) {
    for (i in (p + 1):(p + s + 1)) {
      text1 <- as.expression(bquote(bold(beta)^(.(i - p - 
                                                    s1))))
      text <- append(text, text1)
    }
  }
  f <- function(m) t(m)[, nrow(m):1]
  rgb.palette <- colorRampPalette(c("blue", "white", "red"), space = "Lab")
  at <- seq(k/2 + 0.5, p * (k) + 0.5, by = k)
  if (m > 0) {
    at2 <- seq(p * k + m/2 + 0.5, p * k + s * m + 0.5, by = m)
  }
  else {
    at2 = c()
  }
  at <- c(at, at2)
  se2 = seq(1.75, by = k, length = k)
  L2 <- levelplot(as.matrix(f(B)), at =seq( -1, 1, length=101),col.regions = rgb.palette(100), 
                  colorkey = NULL, xlab = NULL, ylab = NULL, main = list(label = name, 
                                                                         cex = 1), panel = function(...) {
                                                                           panel.levelplot(...)
                                                                           panel.abline(a = NULL, b = 1, h = seq(1.5, m * s + 
                                                                                                                   p * k + 0.5, by = 1), v = seq(1.5, by = 1, length = p * 
                                                                                                                                                   k + m * s), lwd = 0.5)
                                                                           bl1 <- seq(k + 0.5, p * k + 0.5, by = k)
                                                                           b23 <- seq(p * k + 0.5, p * k + 0.5 + s * m, by = m)
                                                                           b1 <- c(bl1, b23)
                                                                           panel.abline(a = NULL, b = 1, v = p * k + 0.5, lwd = 3)
                                                                           panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
                                                                         }, scales = list(x = list(alternating = 1, labels = text, 
                                                                                                   cex = 1, at = at, tck = c(0, 0)), y = list(alternating = 0, 
                                                                                                                                              tck = c(0, 0))))
  return(L2)
}



#' Generating non-stationary ARMA data.
#' 
#' @param nobs number of time points
#' @param arlags the true AR order
#' @param malags the true MA order
#' @param cnst the constant
#' @param phi parameter matrix of the AR model
#' @param theta parameter matrix of the MA model
#' @param sigma covariance matrix of the white noise
#' @return Matrices of time series data and white noise data 
var.sim.break.tdist <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL,
                                 theta = NULL, skip = 200, sigma, brk = nobs+1, df = 1 ) {
  if (!is.matrix(sigma)) 
    sigma = as.matrix(sigma)
  k <- nrow(sigma); m <- length(brk); nT <- nobs + skip
  
  #generate multivariate normal distributed data as the white noise data
  # at <- rmvnorm(nT, rep(0, k), sigma)
  at <- rmvt(nT, sigma = sigma, df = df)
  
  #generate the ARMA time series data
  nar <- length(arlags); p <- 0
  if (nar > 0) {
    arlags <- sort(arlags)
    p <- arlags[nar]
  }
  
  nma <- length(malags); q <- 0
  if (nma > 0) {
    malags <- sort(malags)
    q <- malags[nma]
  }
  
  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  if (length(cnst) == 0) 
    cnst = rep(0, k)
  if (m == 1){
    for (it in ist:nT) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
  }
  
  #if there are some break points
  if (m > 1){
    for (it in ist:(skip+brk[1]-1)) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
    for ( mm in 1:(m-1)){
      for (it in (skip+brk[mm]):(skip+brk[mm+1]-1) ) {
        tmp = matrix(at[it, ], 1, k)
        if (nma > 0) {
          for (j in 1:nma) {
            jdx = (j - 1) * k
            thej = theta[, (jdx + 1):(jdx + k)]
            atm = matrix(at[it - malags[j], ], 1, k)
            tmp = tmp - atm %*% t(thej)
          }
        }
        if (nar > 0) {
          for (i in 1:nar) {
            idx = (i - 1) * k
            phj = phi[, ((mm)*p*k+idx + 1):((mm)*p*k+idx + k)]
            ztm = matrix(zt[it - arlags[i], ], 1, k)
            tmp = tmp + ztm %*% t(phj)
          }
        }
        zt[it, ] = cnst + tmp
      }
    }
  }
  
  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsim <- list(series = zt, noises = at)
}



#' Generating non-stationary ARMA data.
#' 
#' @param nobs number of time points
#' @param arlags the true AR order
#' @param malags the true MA order
#' @param cnst the constant
#' @param phi parameter matrix of the AR model
#' @param theta parameter matrix of the MA model
#' @param sigma covariance matrix of the white noise
#' @return Matrices of time series data and white noise data 
var.sim.break <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL,theta = NULL, skip = 200, sigma, brk = nobs+1) {
  if (!is.matrix(sigma)) 
    sigma = as.matrix(sigma)
  k <- nrow(sigma); m <- length(brk); nT <- nobs + skip
  
  #generate multivariate normal distributed data as the white noise data
  at <- rmvnorm(nT, rep(0, k), sigma)
  
  #generate the ARMA time series data
  nar <- length(arlags); p <- 0
  if (nar > 0) {
    arlags <- sort(arlags)
    p <- arlags[nar]
  }
  
  nma <- length(malags); q <- 0
  if (nma > 0) {
    malags <- sort(malags)
    q <- malags[nma]
  }
  
  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  if (length(cnst) == 0) 
    cnst = rep(0, k)
  if (m == 1){
    for (it in ist:nT) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
  }
  
  #if there are some break points
  if (m > 1){
    for (it in ist:(skip+brk[1]-1)) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
    for ( mm in 1:(m-1)){
      for (it in (skip+brk[mm]):(skip+brk[mm+1]-1) ) {
        tmp = matrix(at[it, ], 1, k)
        if (nma > 0) {
          for (j in 1:nma) {
            jdx = (j - 1) * k
            thej = theta[, (jdx + 1):(jdx + k)]
            atm = matrix(at[it - malags[j], ], 1, k)
            tmp = tmp - atm %*% t(thej)
          }
        }
        if (nar > 0) {
          for (i in 1:nar) {
            idx = (i - 1) * k
            phj = phi[, ((mm)*p*k+idx + 1):((mm)*p*k+idx + k)]
            ztm = matrix(zt[it - arlags[i], ], 1, k)
            tmp = tmp + ztm %*% t(phj)
          }
        }
        zt[it, ] = cnst + tmp
      }
    }
  }
  
  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsim <- list(series = zt, noises = at)
}


#' block segmentation scheme (BSS).
#' 
#' @description Perform the block segmentation scheme (BSS) algorithm to detect the structural breaks 
#' in large scale high-dimensional non-stationary VAR models.
#' 
#' @param data input data matrix, with each column representing the time series component 
#' @param lambda.1.cv tuning parmaeter lambda_1 for fused lasso
#' @param lambda.2.cv tuning parmaeter lambda_2 for fused lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param block.size the block size
#' @param blocks the blocks
#' @param refine logical; if TRUE, use local refinement in the exhaustive search step. Default is TRUE.
#' @return A list object, which contains the followings
#' \describe{
#'   \item{pts.1}{a set of selected break point after the first block fused lasso step}
#'   \item{pts.2}{a set of selected break point after the second local screening step}
#'   \item{pts.3}{a set of selected break point after the thrid exhaustive search step}
#'   \item{an}{the selected neighborhood size a_n after the grid search}
#' }
bss <- function(data, lambda.1.cv = NULL, lambda.2.cv = NULL, q = 1, 
                max.iteration = 100, tol = 10^(-2), block.size = NULL, blocks = NULL,
                refine = FALSE, use.BIC= FALSE, an.grid = NULL){
  T <- length(data[,1]); p <- length(data[1,]); 
  second.brk.points <- c(); pts.final <- c();
  
  ############# block size and blocks ###########
  if(is.null(block.size) && is.null(blocks) ){
    block.size = floor(sqrt(T));
    blocks <- seq(0, T, block.size);
  }else if( !is.null(block.size) && is.null(blocks)){
    blocks <- seq(0, T, block.size);
  }else if(!is.null(block.size) && !is.null(blocks)){
    #check if the block.size and blocks match
    n.new <- length(blocks) - 1;
    blocks.size.check <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
    if( sum(blocks.size.check[1: (length(blocks.size.check)-1 )] != block.size ) > 0 ){
      stop("Error: The block.size and blocks can't match!")
    }
  }
  
  if(blocks[length(blocks)] < T){
    blocks <- c(blocks[-length(blocks)], T)
  }
  
  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  
  
  #sample the cv index for cross-validation
  bbb <- floor(n.new/5);
  # aaa <- sample(1:5, 1);
  aaa <- 4
  cv.index <- seq(aaa, n.new, floor(n.new/bbb));
  
  ############# Tuning parameter ################
  if(is.null(lambda.1.cv)){
    lambda.1.max <- lambda_warm_up(data, q, blocks, cv.index)$lambda_1_max
    if(blocks[2] <= 2*p ){
      epsilon <-  10^(-3)
    }
    if(blocks[2] >= 2*p ){
      epsilon <-  10^(-4)
    }
    nlam <- 10 
    lambda.1.min <-  lambda.1.max*epsilon
    delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
    lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
  }
  
  if(is.null(lambda.2.cv)){
    lambda.2.cv <-  c(10*sqrt(log(p)/T),1*sqrt(log(p)/T),0.10*sqrt(log(p)/T))
  }
  
  print("lambda.1.cv:")
  print(lambda.1.cv)
  print("lambda.2.cv:")
  print(lambda.2.cv)
  ######################################################################
  ######## First Step: Initial Break Points Selection ##################
  ######################################################################
  
  #run the first block fused lasso step 
  temp.first <- first.step.blocks( data, lambda.1.cv, lambda.2.cv, q, max.iteration = max.iteration, tol = tol, cv.index, blocks=blocks)
  first.brk.points <- temp.first$brk.points;
  print(first.brk.points)
  phi.est.full <- temp.first$phi.full
  print(temp.first$cv)
  print(temp.first$cv1.final)
  print(temp.first$cv2.final)

  #construct the grid values of neighborhood size a_n
  n <- T - q;
  an.lb <- max(floor(mean(blocks.size)),floor( (log(n)*log(p))^1 ));
  an.ub <-  min(10*an.lb,0.95*(min(first.brk.points)-1-q),0.95*(n - max(first.brk.points)-1))
  if(is.null(an.grid)){
    an.grid <- seq(an.lb,  an.ub, length.out = 5);
  }
 
  an.idx.final <- length(an.grid)
  an.grid <- floor(an.grid);
  final.pts.res <- vector("list",length(an.grid));
  final.phi.hat.list.res <- vector("list",length(an.grid));
  flag <- c("FALSE");
  an.idx <- 0;
  phi.local.1.full <- vector("list", length(an.grid));
  phi.local.2.full <- vector("list", length(an.grid));
  
  #for each a_n, run the second and thrid step
  while(an.idx < length(an.grid)  ){
    an.idx <- an.idx + 1;
    an <- an.grid[an.idx]
    print("an:")
    print(an)
    
    #remove the boundary points
    remove.ind <- c();
    if(length(first.brk.points) != 0){
      for(i in 1:length(first.brk.points)){
        if ( first.brk.points[i] < (an-1-q)   ){remove.ind <- c(remove.ind,i);}
        if ( (T-first.brk.points[i]) < (an-1-q)   ){remove.ind <- c(remove.ind,i);}
      }
    }
    if( length(remove.ind) > 0  ){first.brk.points <- first.brk.points[-remove.ind];}
    
    #if there are selected break points after the first step
    if( length(first.brk.points) != 0){
      
      #####################################################
      ######## Second Step: Local Screening      ##########
      #####################################################
      eta <- (1/1)*(log(2*an)*log(p))/(2*an); # the tuning parameter for second and third steps.
      #run the second local screening step
      temp <- second.step.local(data, eta = eta, q, max.iteration = 1000, tol = tol, first.brk.points, an, 
                                phi.est.full = phi.est.full, blocks, use.BIC)
      
      #record the selected break points after local screening step
      second.brk.points <- temp$pts;
      phi.local.1.full[[an.idx]] <- temp$phi.local.1
      phi.local.2.full[[an.idx]] <- temp$phi.local.2
      
      ######################################################
      ######## Thrid Step: Exhaustive Search      ##########
      ######################################################
      print("second.brk.points:")
      print(second.brk.points)
      pts.final <- second.brk.points;
      #keep running until none of the selected break points close to any other selected break points 
      while( min(abs(diff(pts.final)), 3*an ) <=  2*an  ){
        if( length(pts.final) != 0){
          #cluster the selected break points by size 2a_n
          pts.list <- block.finder(pts.final, 2*an)
          # run the third exhaustive search step for each cluster
          print(pts.list)
          temp.third <- third.step.exhaustive(data, q, max.iteration = 1000, tol = tol, pts.list, an, eta)
          pts.final <- temp.third
          # print("pts.final:")
          # print(pts.final)
         
        }
      }

      #record the final selected break points for each given a_n
      final.pts.res[[an.idx]] <- pts.final
      # final.phi.hat.list.res[[an.idx]] <- phi.hat.list
      
      #terminate the grid search of an if the number of final selected break points is stable
      if(an.idx > 2){
        if( length(final.pts.res[[an.idx]]) == length(final.pts.res[[an.idx-1]]) && length(final.pts.res[[an.idx-1]]) == length(final.pts.res[[an.idx-2]]) ){
          flag <- c("TRUE");
          an.idx.final <- an.idx;
          an.sel <- an.grid[an.idx];
          break;
        }
      }
      
    }
  }
  
  #if the stable criterion hasn't been met
  #find the length that happen the most
  #if there are multiple lengths with same occurrence, find the longest one
  if(flag == FALSE){
    loc.final <- rep(0,length(an.grid));
    for(i in 1:length(an.grid)){
      loc.final[i] <- length(final.pts.res[[i]]);
    }
    loc.table <- table(loc.final)
    counts.final <- sort(loc.table,decreasing=TRUE)[1]
    if(counts.final >=3){
      len.final <- max(as.integer(names(loc.table)[loc.table == counts.final]))
     
      
    }else{
      # choose the longest one instead
      len.final <- max(loc.final)
    }
    an.idx.final <- max(c(1:length(loc.final))[loc.final == len.final])
    an.sel <- an.grid[an.idx.final];
    
  }
  
  if(refine == TRUE){
    pts.final = final.pts.res[[an.idx.final]]
    local.idx = match(pts.final, first.brk.points)
    print("pts.final:")
    print(pts.final)
    # print("local.idx:")
    # print(local.idx)
    phi.local.1 <- phi.local.1.full[[an.idx.final]]
    phi.local.2 <- phi.local.2.full[[an.idx.final]]
    print(length(phi.local.1))
    if(length(pts.final) > 0){
      pts.list <- block.finder(pts.final, 2*an.sel)
      # temp.forth <- forth.step.refine(data, q, max.iteration = 1000, tol = tol, pts.list, an.sel, phi.est.full= phi.est.full, blocks)
      temp.forth <- forth.step.refine(data, q, max.iteration = 1000, tol = tol, pts.list,
                                      an.sel, phi.est.full= phi.est.full,
                                      phi.local.1 = phi.local.1[local.idx],
                                      phi.local.2 = phi.local.2[local.idx],
                                      blocks)
      phi.hat.list <- temp.forth$phi.hat.list
      pts.final <- temp.forth$pts
      
    }else{
      idx <- floor(n.new/2)
      phi.hat.list <- phi.est.full[[idx]]
    }
    
  }else{
    pts.final = final.pts.res[[an.idx.final]]
    print("an.sel")
    print(an.sel)
    print(pts.final)
    if(length(pts.final) > 0){
      pts.list <- block.finder(pts.final, 2*an.sel)
      pts.list.full <- pts.list
      pts.list.full <- c(1, pts.list.full , T)
      phi.hat.list <- vector("list", length(pts.final) + 1)
      cp.index.list <- vector("list", length(pts.final) + 2);
      cp.index.list[[1]] <- c(1);
      cp.index.list[[length(pts.final)+2]] <- c(n.new+1);
      for(i in 1:length(pts.final)){
        pts.temp <- pts.list.full[[i+1]];
        cp.index.list[[i+1]] <- match(pts.temp, blocks)
      }
      print("n.new")
      print(n.new)
      # print(cp.index.list)
      
      #construct the interval for performing the lasso and computing the loss function
      for(i in 1:length(pts.final)){
        idx <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2);
        phi.hat.list[[i]] <- phi.est.full[[idx]]
        print(idx)
      }
      idx <- floor((min(cp.index.list[[length(pts.final)+2]]) + max(cp.index.list[[length(pts.final)+1]]))/2);
      print(idx)
      phi.hat.list [[length(pts.final)+1]] <- phi.est.full[[idx]]
      
    }else{
      idx <- floor(n.new/2)
      phi.hat.list <- phi.est.full[[idx]]
    }
    
    
      
    
  }
  
  
  
  return(list(first.selected.points = first.brk.points, second.selected.points = second.brk.points, 
              final.selected.points = pts.final, 
              final.selected.points.grid = final.pts.res, 
              an = an.sel, phi.est.full = phi.est.full, final.phi.hat.list = phi.hat.list )) 
}


bss.timing <- function(data, lambda.1.cv = NULL, lambda.2.cv = NULL, q = 1, 
                max.iteration = 100, tol = 10^(-2), block.size = NULL, blocks = NULL,
                refine = FALSE, use.BIC= FALSE, an.grid = NULL){
  T <- length(data[,1]); p <- length(data[1,]); 
  second.brk.points <- c(); pts.final <- c();
  
  ############# block size and blocks ###########
  if(is.null(block.size) && is.null(blocks) ){
    block.size = floor(sqrt(T));
    blocks <- seq(0, T, block.size);
  }else if( !is.null(block.size) && is.null(blocks)){
    blocks <- seq(0, T, block.size);
  }else if(!is.null(block.size) && !is.null(blocks)){
    #check if the block.size and blocks match
    n.new <- length(blocks) - 1;
    blocks.size.check <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
    if( sum(blocks.size.check[1: (length(blocks.size.check)-1 )] != block.size ) > 0 ){
      stop("Error: The block.size and blocks can't match!")
    }
  }
  
  if(blocks[length(blocks)] < T){
    blocks <- c(blocks[-length(blocks)], T)
  }
  
  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  
  
  #sample the cv index for cross-validation
  bbb <- floor(n.new/5);
  # aaa <- sample(1:5, 1);
  aaa <- 4
  cv.index <- seq(aaa, n.new, floor(n.new/bbb));
  
  ############# Tuning parameter ################
  if(is.null(lambda.1.cv)){
    lambda.1.max <- lambda_warm_up(data, q, blocks, cv.index)$lambda_1_max
    if(blocks[2] <= 2*p ){
      epsilon <-  10^(-3)
    }
    if(blocks[2] >= 2*p ){
      epsilon <-  10^(-4)
    }
    nlam <- 10 
    lambda.1.min <-  lambda.1.max*epsilon
    delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
    lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
  }
  
  if(is.null(lambda.2.cv)){
    lambda.2.cv <-  c(10*sqrt(log(p)/T),1*sqrt(log(p)/T),0.10*sqrt(log(p)/T))
  }
  
  print("lambda.1.cv:")
  print(lambda.1.cv)
  print("lambda.2.cv:")
  print(lambda.2.cv)
  ######################################################################
  ######## First Step: Initial Break Points Selection ##################
  ######################################################################
  time.comparison <- rep(0, 4)
  #run the first step
  ptm.temp <- proc.time()
  
  
  #run the first block fused lasso step 
  temp.first <- first.step.blocks( data, lambda.1.cv, lambda.2.cv, q, max.iteration = max.iteration, tol = tol, cv.index, blocks=blocks)
  
  time.temp <- proc.time() - ptm.temp;
  time.comparison[1] <- c(time.temp[3])
  
  
  first.brk.points <- temp.first$brk.points;
  print(first.brk.points)
  phi.est.full <- temp.first$phi.full
  print(temp.first$cv)
  print(temp.first$cv1.final)
  print(temp.first$cv2.final)
  
  #construct the grid values of neighborhood size a_n
  n <- T - q;
  an.lb <- max(floor(mean(blocks.size)),floor( (log(n)*log(p))^1 ));
  an.ub <-  min(10*an.lb,0.95*(min(first.brk.points)-1-q),0.95*(n - max(first.brk.points)-1))
  if(is.null(an.grid)){
    an.grid <- seq(an.lb,  an.ub, length.out = 5);
  }
  
  an.idx.final <- length(an.grid)
  an.grid <- floor(an.grid);
  final.pts.res <- vector("list",length(an.grid));
  final.phi.hat.list.res <- vector("list",length(an.grid));
  flag <- c("FALSE");
  an.idx <- 0;
  phi.local.1.full <- vector("list", length(an.grid));
  phi.local.2.full <- vector("list", length(an.grid));
  
  
  
  
  #for each a_n, run the second and thrid step
  while(an.idx < length(an.grid)  ){
    an.idx <- an.idx + 1;
    an <- an.grid[an.idx]
    print("an:")
    print(an)
    
    #remove the boundary points
    remove.ind <- c();
    if(length(first.brk.points) != 0){
      for(i in 1:length(first.brk.points)){
        if ( first.brk.points[i] < (an-1-q)   ){remove.ind <- c(remove.ind,i);}
        if ( (T-first.brk.points[i]) < (an-1-q)   ){remove.ind <- c(remove.ind,i);}
      }
    }
    if( length(remove.ind) > 0  ){first.brk.points <- first.brk.points[-remove.ind];}
    
    #if there are selected break points after the first step
    if( length(first.brk.points) != 0){
      
      #####################################################
      ######## Second Step: Local Screening      ##########
      #####################################################
      eta <- (1/1)*(log(2*an)*log(p))/(2*an); # the tuning parameter for second and third steps.
      #run the second local screening step
      
      ptm.temp <- proc.time()
      
      temp <- second.step.local(data, eta = eta, q, max.iteration = 1000, tol = tol, first.brk.points, an, 
                                phi.est.full = phi.est.full, blocks, use.BIC)
      
      time.temp <- proc.time() - ptm.temp;
      time.comparison[2] <- time.comparison[2] + c(time.temp[3])
      
      #record the selected break points after local screening step
      second.brk.points <- temp$pts;
      phi.local.1.full[[an.idx]] <- temp$phi.local.1
      phi.local.2.full[[an.idx]] <- temp$phi.local.2
      
      ######################################################
      ######## Thrid Step: Exhaustive Search      ##########
      ######################################################
      print("second.brk.points:")
      print(second.brk.points)
      pts.final <- second.brk.points;
      
      
      ptm.temp <- proc.time()
      
      #keep running until none of the selected break points close to any other selected break points 
      while( min(abs(diff(pts.final)), 3*an ) <=  2*an  ){
        if( length(pts.final) != 0){
          #cluster the selected break points by size 2a_n
          pts.list <- block.finder(pts.final, 2*an)
          # run the third exhaustive search step for each cluster
          print(pts.list)
          temp.third <- third.step.exhaustive(data, q, max.iteration = 1000, tol = tol, pts.list, an, eta)
          pts.final <- temp.third
          # print("pts.final:")
          # print(pts.final)
          
        }
      }
      
      time.temp <- proc.time() - ptm.temp;
      time.comparison[3] <- time.comparison[3] + c(time.temp[3])
      
      
      #record the final selected break points for each given a_n
      final.pts.res[[an.idx]] <- pts.final
      # final.phi.hat.list.res[[an.idx]] <- phi.hat.list
      
      #terminate the grid search of an if the number of final selected break points is stable
      if(an.idx > 2){
        if( length(final.pts.res[[an.idx]]) == length(final.pts.res[[an.idx-1]]) && length(final.pts.res[[an.idx-1]]) == length(final.pts.res[[an.idx-2]]) ){
          flag <- c("TRUE");
          an.idx.final <- an.idx;
          an.sel <- an.grid[an.idx];
          break;
        }
      }
      
    }
  }
  
  #if the stable criterion hasn't been met
  #find the length that happen the most
  #if there are multiple lengths with same occurrence, find the longest one
  if(flag == FALSE){
    loc.final <- rep(0,length(an.grid));
    for(i in 1:length(an.grid)){
      loc.final[i] <- length(final.pts.res[[i]]);
    }
    loc.table <- table(loc.final)
    counts.final <- sort(loc.table,decreasing=TRUE)[1]
    if(counts.final >=3){
      len.final <- max(as.integer(names(loc.table)[loc.table == counts.final]))
      
      
    }else{
      # choose the longest one instead
      len.final <- max(loc.final)
    }
    an.idx.final <- max(c(1:length(loc.final))[loc.final == len.final])
    an.sel <- an.grid[an.idx.final];
    
  }
  
  if(refine == TRUE){
    pts.final = final.pts.res[[an.idx.final]]
    local.idx = match(pts.final, first.brk.points)
    print("pts.final:")
    print(pts.final)
    # print("local.idx:")
    # print(local.idx)
    phi.local.1 <- phi.local.1.full[[an.idx.final]]
    phi.local.2 <- phi.local.2.full[[an.idx.final]]
    print(length(phi.local.1))
    if(length(pts.final) > 0){
      pts.list <- block.finder(pts.final, 2*an.sel)
      # temp.forth <- forth.step.refine(data, q, max.iteration = 1000, tol = tol, pts.list, an.sel, phi.est.full= phi.est.full, blocks)
      temp.forth <- forth.step.refine(data, q, max.iteration = 1000, tol = tol, pts.list,
                                      an.sel, phi.est.full= phi.est.full,
                                      phi.local.1 = phi.local.1[local.idx],
                                      phi.local.2 = phi.local.2[local.idx],
                                      blocks)
      phi.hat.list <- temp.forth$phi.hat.list
      pts.final <- temp.forth$pts
      
    }else{
      idx <- floor(n.new/2)
      phi.hat.list <- phi.est.full[[idx]]
    }
    
  }else{
    pts.final = final.pts.res[[an.idx.final]]
    print("an.sel")
    print(an.sel)
    print(pts.final)
    if(length(pts.final) > 0){
      pts.list <- block.finder(pts.final, 2*an.sel)
      pts.list.full <- pts.list
      pts.list.full <- c(1, pts.list.full , T)
      phi.hat.list <- vector("list", length(pts.final) + 1)
      cp.index.list <- vector("list", length(pts.final) + 2);
      cp.index.list[[1]] <- c(1);
      cp.index.list[[length(pts.final)+2]] <- c(n.new+1);
      for(i in 1:length(pts.final)){
        pts.temp <- pts.list.full[[i+1]];
        cp.index.list[[i+1]] <- match(pts.temp, blocks)
      }
      print("n.new")
      print(n.new)
      # print(cp.index.list)
      
      #construct the interval for performing the lasso and computing the loss function
      for(i in 1:length(pts.final)){
        idx <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2);
        phi.hat.list[[i]] <- phi.est.full[[idx]]
        print(idx)
      }
      idx <- floor((min(cp.index.list[[length(pts.final)+2]]) + max(cp.index.list[[length(pts.final)+1]]))/2);
      print(idx)
      phi.hat.list [[length(pts.final)+1]] <- phi.est.full[[idx]]
      
    }else{
      idx <- floor(n.new/2)
      phi.hat.list <- phi.est.full[[idx]]
    }
    
    
    
    
  }
  
  
  
  return(list(first.selected.points = first.brk.points, second.selected.points = second.brk.points, 
              final.selected.points = pts.final, 
              final.selected.points.grid = final.pts.res, 
              an = an.sel, phi.est.full = phi.est.full, final.phi.hat.list = phi.hat.list,
              timing = time.comparison)) 
}


#' block fused lasso step (first step).
#' 
#' @description Perform the block fused lasso to detect candidate break points.
#' 
#' @param data.temp input data matrix, with each column representing the time series component 
#' @param lambda.1.cv tuning parmaeter lambda_1 for fused lasso
#' @param lambda.2.cv tuning parmaeter lambda_2 for fused lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param cv.index the index of time points for cross-validation
#' @param blocks the blocks
#' @return A list object, which contains the followings
#' \describe{
#'   \item{brk.points}{a set of selected break point after the first block fused lasso step}
#'   \item{cv}{the cross validation values for tuning parmeter selection}
#'   \item{cv1.final}{the selected lambda_1}
#'   \item{cv2.final}{the selected lambda_2}
#' }

first.step.blocks <- function(data.temp, lambda.1.cv, lambda.2.cv, q, max.iteration = max.iteration, tol = tol,cv.index, blocks){
  
  cv.l <- length(cv.index); data.org <- data.temp; T.org <- length(data.temp[,1]); p <- length(data.temp[1,]); n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  
  #create the tuning parmaeter combination of lambda1 and lambda2
  lambda.full <- expand.grid(lambda.1.cv,lambda.2.cv)
  kk <- length(lambda.full[,1]);
  
  cv <- rep(NA,kk); 
  phi.final <- vector("list",kk);
  T <- length(data.temp[,1]); p <- length(data.temp[1,]);
  brk.points.final <- vector("list",kk);
  flag.full <- rep(0,kk);
  
  #cross-validation for each values of lambda1 and lambda2
  nlam1 <- length(lambda.1.cv)
  nlam2 <- length(lambda.2.cv)
  kk <- nlam1*nlam2
  i = 1
  while(i <= kk) {
    print(i)
    i.lam1 <- i%% nlam1
    if(i.lam1 == 0 ){
      i.lam1 = nlam1
    }
    i.lam2 <- floor((i-1)/nlam1)+1
    if ( i == 1){
      test <- var_break_fit_block_cpp(data.temp, lambda.full[i,1],lambda.full[i,2], q, max.iteration, tol = tol, initial_phi = 0.0+matrix(0.0,p,p*q*n.new), blocks, cv.index)
      flag.full[i] <- test$flag;
    }else if(is.na(cv[i-1]) ){
      test <- var_break_fit_block_cpp(data.temp, lambda.full[i,1],lambda.full[i,2], q, max.iteration, tol = tol, initial_phi = 0.0+matrix(0.0,p,p*q*n.new), blocks, cv.index)
      flag.full[i] <- test$flag;
    }else{
      initial.phi <- phi.final[[(i-1)]]
      if(max(abs(phi.final[[(i-1)]])) > 10^3  ){initial.phi <- 0*phi.final[[(i-1)]];}
      test <- var_break_fit_block_cpp(data.temp, lambda.full[i,1],lambda.full[i,2], q, max.iteration, tol = tol, initial_phi = initial.phi, blocks, cv.index)
      flag.full[i] <- test$flag;
    }
    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;
    
    ll <- c(0);
    brk.points.list <- vector("list",length(ll));
    
    for(j in 1:length(ll)){
      phi.hat <- phi.hat.full;
      n <- T - q;
      m.hat <- 0; brk.points <- rep(0,n.new);
      
      for (iii in 1:(n.new-1))
      {
        if ( sum((phi.hat[,((iii-1)*p*q+1):(iii*p*q)] )^2 ) > tol   ){
          m.hat <- m.hat + 1; brk.points[m.hat] <- blocks[iii+1];
        }
      }
      
      
      loc <- rep(0,m.hat);
      brk.points <- brk.points[1:m.hat];
      
      #remove the bounary points and clean up 
      brk.points <- brk.points[which(brk.points > 3*mean(blocks.size))]; brk.points <- brk.points[which(brk.points < (n-3*mean(blocks.size)))];
      m.hat <- length(brk.points);
      del <- 0;
      if(m.hat >= 2){
        while(del < m.hat){
          if(length(brk.points) <= 1){break;}
          del <- del + 1;
          deleted <- 0;
          for (i.3 in 2:length(brk.points)) {
            if(deleted == 0 &&  abs(brk.points[i.3] - brk.points[i.3-1]) <= (1)*max(q,min(blocks.size))  ){
              brk.points <- brk.points[-i.3]; deleted <- 1;
            }
          }
        }
      }
      
      brk.points.list[[j]] <- brk.points;
    }
    
    brk.points.final[[i]] <- brk.points;
    m.hat <- length(brk.points);
    
    #forecast the time series based on the estimated matrix Phi
    #and compute the forecast error
    phi.full.all <- vector("list",n.new);
    forecast <- matrix(0,p,T);
    phi.full.all[[1]] <- phi.hat[,(1):(p*q)];
    for(i.1 in 2:n.new){
      phi.full.all[[i.1]] <- phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*p*q+1):(i.1*p*q)];
      forecast[,(blocks[i.1]+1):(blocks[i.1+1])] <- pred.block(t(data.org),phi.full.all[[i.1-1]],q,blocks[i.1],p,blocks[i.1+1]-blocks[i.1]);
    }
    forecast.new <- matrix(0,p,cv.l);
    for(j in (1):cv.l){
      forecast.new[,j] <- pred(t(data.org),phi.full.all[[(cv.index[j])]],q,blocks[cv.index[j]+1]-1,p,1)
    }
    temp.index <- rep(0,cv.l);
    for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1];}
    cv[i] <- (1/(p*cv.l))*sum( (forecast.new - t(data.org[temp.index,])  )^2 );
    
    #break condition 
    if(nlam1 >= 2){
      if( !(i %in% c(seq(1, kk, nlam1), seq(2, kk, nlam1)))  && cv[i] > cv[i-1] && cv[i-1] > cv[i-2]){
        i.lam2 <- i.lam2 + 1
        i <- (i.lam2-1)*nlam1 + 1
      }else{
        i <- i + 1
      }
      
    }else{
      i <- i+1
    }
  }
  
  #select the tuning parmaete that has the small cross-validation value
  lll <- min(which(cv == min(cv, na.rm = TRUE)));
  phi.hat.full <- phi.final[[lll]];
  #compute the estimated phi
  phi.par.sum <- vector("list", n.new);
  phi.par.sum[[1]] <- phi.hat.full[, 1:(p*q)];
  for(i in 2:n.new){
    phi.par.sum[[i]] <- phi.par.sum[[i-1]] + phi.hat.full[,((i-1)*p*q+1):(i*p*q)];
  }
  

  return(list(brk.points = brk.points.final[[lll]], cv = cv, 
              cv1.final = lambda.full[lll,1], cv2.final = lambda.full[lll,2],
              phi.full = phi.par.sum
              ))
}


#' BIC  and HBIC function
#' @param residual residual matrix
#' @param phi estimated coefficient matrix of the model
#' @param gamma.val hyperparameter for HBIC, if HBIC == TRUE.
BIC <- function(residual, phi, gamma.val = 1){
  p<- length(phi[, 1]);
  q <- length(phi[1, ])/p;
  T.new <- length(residual[1, ]);
  count <- 0;
  # for (i in 1:p){
  #   for (j in 1:(p*q)){
  #     if(phi[i,j] != 0){
  #       count <- count + 1;
  #     }
  #   }
  # }
  count = sum(phi !=0)
  print("nonzero count"); print(count)
  print("p:")
  print(p)
  
  sigma.hat <- 0*diag(p);
  for(i in 1:T.new){sigma.hat <- sigma.hat +  residual[, i]%*%t(residual[, i]);  }
  sigma.hat <- (1/(T.new))*sigma.hat;
  ee.temp <- min(eigen(sigma.hat)$values);
  if(ee.temp <= 10^(-8)){
    # print("nonpositive eigen values!")
    sigma.hat <- sigma.hat + (2.0)*(abs(ee.temp) + 10^(-3))*diag(p);
  }
  
  log.det <- log(det(sigma.hat));
  count <- count
  print("log.det:")
  print(log.det)
  print("BIC:")
  print(log.det + log(T.new)*count/T.new)
  print("HBIC:")
  print(log.det + 2*gamma.val*log(p*q*p)*count/T.new)
  return(list(BIC = log.det + log(T.new)*count/T.new , HBIC = log.det + 2*gamma.val*log(p*q*p)*count/T.new))
}



#' local sreening step (second step).
#' 
#' @description Perform the local screening to "thin out" redundant break points. 
#' 
#' @param data input data matrix, with each column representing the time series component 
#' @param eta tuning parameter eta for lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param pts the selected break points after the first step
#' @param an the neighborhood size a_n
#' @return A list object, which contains the followings
#' \describe{
#'   \item{pts}{a set of selected break point after the second local screening step}
#'   \item{omega}{the selected Omega value}
#' }
second.step.local <- function(data, eta, q, max.iteration = 1000, tol = 10^(-4), pts, an, 
                              phi.est.full = NULL, blocks = NULL, use.BIC = FALSE){
  m <- length(pts); if( m == 0){break;}
  p = length(data[1,])
  
  #compute the local loss functions for each selected break points 
  try <- break.var.local.new(data, eta, q, max.iteration = 1000, tol = tol, pts, an);
  # record the local loss function that include or exclude some break point 
  L.n.1 = try$L.n.1; L.n.2 = try$L.n.2; 
  #record the local estimate left (1) and right (2)
  phi.local.1 = try$phi.local.1
  phi.local.2 = try$phi.local.2
  
  #OMEGA is selected by data-driven method
  #first, compute the V value as the difference of loss functions that include and exclude some break point 
  V = rep(0, m)
  for(i in 1:m ){
    V[i] = L.n.2[i] - (L.n.1[2*i-1] +  L.n.1[2*i])
  }
  
  #add two bounary points as reference points (by assumption, their V values shoule be extremly small)
  T <- length(data[,1]);
  pts.redundant <- c(an+q, T-an)
  try.redundant <- break.var.local.new(data, eta, q, max.iteration = 1000, tol = tol,  pts.redundant, an);
  L.n.1.redundant = try.redundant$L.n.1; L.n.2.redundant = try.redundant$L.n.2;
  V.redundant <- rep(0, 2)
  for(i in 1:2 ){
    V.redundant[i] = L.n.2.redundant[i] - (L.n.1.redundant[2*i-1] +  L.n.1.redundant[2*i])
  }
  
  #use the maximum value of V.redundant as the reference V value
  # V <- c(V,rep(max(V.redundant),floor(2*length(V))))
  V <- c(V, rep(max(V.redundant), 2))
  print("V:")
  print(V)
  
  if(use.BIC == FALSE){
    if( length(unique(V)) <= 2 ){
      omega <- max(V) + 10^(-6);
    }
    if( length(unique(V)) > 2 ){
      #use kmeans to cluster the V 
      clus.2 <- kmeans(V, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss; 
      if(fit.2 < 0.20){
        omega <- max(V) + 10^(-6);
      }
      if( fit.2 >= 0.20 ){
        #if the reference point is in the subset with larger center, this means no change points: set omeage = max(V)
        #otherwise, set omega = min(V) -1 
        loc <- clus.2$cluster;
        if( clus.2$centers[1] > clus.2$centers[2]  ){
          omega <- min(V[which(loc==1)]) - 10^(-6) ;
          if(loc[length(loc)] == 1){
            omega <- max(V) + 10^(-6);
          }
        }
        if( clus.2$centers[1] < clus.2$centers[2]  ){
          omega <- min(V[which(loc==2)]) - 10^(-6) ;
          if(loc[length(loc)] == 2){
            omega <- max(V) + 10^(-6);
          }
        }
      }
    }
    
  }
  
  
  
  if(use.BIC == TRUE){
    n.new <- length(blocks) - 1
    ###### use BIC to determine the k-means
    BIC.diff <- 1
    BIC.old <- 10^8
    pts.sel <- c()
    omega <- 0
    loc.block.full <- c()
    while(BIC.diff > 0 & length(unique(V)) > 1 ){
      pts.sel.old <- pts.sel
      omega.old <- omega
      #use kmeans to cluster the V 
      clus.2 <- kmeans(V, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss; 
      print(fit.2)
      if(fit.2 < 0.20){
        #no change points: set omeage = max(V)
        omega <- max(V) + 10^(-6);
        pts.sel <- c(pts.sel);
        break
      }
      if( fit.2 >= 0.20 ){
        #if the reference point is in the subset with larger center, 
        #this means no change points: set omeage = max(V)
        #otherwise, set omega = min(V) -1 
        loc <- clus.2$cluster;
        if( clus.2$centers[1] > clus.2$centers[2]  ){
          omega <- min(V[which(loc==1)]) - 10^(-6) ;
          if(loc[length(loc)] == 1){
            print("large reference! break")
            omega <- max(V)  + 10^(-6);
            pts.sel <- c(pts.sel);
            break
          }else{
            loc.idx <- which(loc==1);
          }
        }
        if( clus.2$centers[1] < clus.2$centers[2]  ){
          omega <- min(V[which(loc==2)]) - 10^(-6) ;
          if(loc[length(loc)] == 2){
            print("large reference! break")
            omega <- max(V) + 10^(-6);
            pts.sel <- c(pts.sel);
            break
          }else{
            loc.idx <- which(loc==2);
          }
        }
        pts.sel <- sort(c(pts.sel, pts[loc.idx]))
        V[loc.idx] <- V[length(V)]
        loc.block.full <- match(pts.sel, blocks)
        print(pts.sel)
      }
      print("pts.sel:")
      print(pts.sel)
      
      m.temp <- length(pts.sel)
      phi.est.new <- vector("list", m.temp + 1);
      cp.index.list <- vector("list", m.temp + 2);
      cp.index.list[[1]] <- c(1);
      cp.index.list[[m.temp+2]] <- c(n.new+1);
      for(i.1 in 1:m.temp){
        pts.temp <- pts.sel[i.1]
        cp.index.list[[i.1+1]] <- match(pts.temp, blocks)
      }
      # print(cp.index.list)
      phi.full.all <- vector("list", n.new);
      for(i.1 in 1:(m.temp + 1)){
        idx <- floor( (cp.index.list[[i.1+1]] +cp.index.list[[i.1]] )/2 )
        phi.est.new[[i.1]] <- matrix(phi.est.full[[idx]], ncol = p*q);
        for(i.2 in cp.index.list[[i.1]]: (cp.index.list[[i.1+1]]-1)){
          phi.full.all[[i.2]] <-   phi.est.new[[i.1]]
          
        }
        
      }
     
     
      forecast.all.new <- matrix(0, p, T);
      forecast.all.new[, (blocks[1]+1+q):(blocks[1+1])] <-
        sapply(c((blocks[1]+1+q):(blocks[1+1])), 
               function(jjj) pred(t(data), phi.full.all[[1]], q, jjj-1 , p, 1) )
      for(i.1 in 2:n.new){
        forecast.all.new[, (blocks[i.1]+1):(blocks[i.1+1])] <-
          sapply(c((blocks[i.1]+1):(blocks[i.1+1])), function(jjj) pred(t(data), phi.full.all[[i.1]], q, jjj-1 , p, 1) )
      }
      residual <- t(data[( (1+q) :T), ]) - forecast.all.new[, (1+q) :T];
      print("Use BIC!")
      BIC.new <- BIC(residual, phi = do.call(cbind, phi.est.new))$BIC
      
      print("BIC.new:"); print(BIC.new)
      BIC.diff <- BIC.old - BIC.new
      print("BIC.diff:");print(BIC.diff)
      BIC.old <- BIC.new
      if(BIC.diff <= 0){
        pts.sel <- sort(pts.sel.old)
        omega <- omega.old
        break
      }
    }
    
    print("BIC stop")
  }
  
  
  
  
  #select the break points by localized information criterion (LIC)
  L.n.1.temp <- L.n.1; L.n.2.temp <- L.n.2; L.n.plot <- rep(0,m+1); L.n.plot[1] <- sum(L.n.1) + m*omega; 
  mm <- 0; ic <- 0; add.temp <- 0; pts.full <- vector("list",m+1); pts.full[[1]] <- pts; ind.pts <- rep(0,m);
  while(mm < m){
    mm <- mm + 1;
    L.n.temp <- rep(0,length(pts));
    for(i in 1:length(pts)){
      L.n.temp[i] <- sum(L.n.1.temp) - L.n.1.temp[(2*i-1)] - L.n.1.temp[(2*i)] + L.n.2.temp[i] + 1*add.temp;
    }
    ll <- min(which.min(L.n.temp)); ind.pts[mm] <- ll;
    pts <- pts[-ll]; 
    L.n.1.temp <- L.n.1.temp[-c(2*ll-1,2*ll)]; add.temp <- add.temp + 1*L.n.2.temp[ll]; 
    L.n.2.temp <- L.n.2.temp[-ll]; 
    L.n.plot[mm+1] <- L.n.temp[ll] + (m - mm)*omega;
    pts.full[[mm+1]] <- pts;
  }
  
  ind <- 0;
  ind <- min(which.min(L.n.plot))
  
  return(list(pts = pts.full[[ind]], omega = omega, 
              phi.local.1= phi.local.1, phi.local.2 = phi.local.2 ))
}

# second.step.local <- function(data, eta, q, max.iteration = 1000, tol = 10^(-4), pts, an){
#   m <- length(pts); if( m == 0){break;}
#   
#   #compute the local loss functions for each selected break points 
#   try <- break.var.local.new(data, eta, q, max.iteration = 1000, tol = tol, pts, an);
#   # record the local loss function that include or exclude some break point 
#   L.n.1 = try$L.n.1; L.n.2 = try$L.n.2; 
#   
#   #OMEGA is selected by data-driven method
#   #first, compute the V value as the difference of loss functions that include and exclude some break point 
#   V = rep(0, m)
#   for(i in 1:m ){
#     V[i] = L.n.2[i] - (L.n.1[2*i-1] +  L.n.1[2*i])
#   }
#   
#   #add two bounary points as reference points (by assumption, their V values shoule be extremly small)
#   T <- length(data[,1]);
#   pts.redundant <- c(an+q, T-an)
#   try.redundant <- break.var.local.new(data, eta, q, max.iteration = 1000, tol = tol,  pts.redundant, an);
#   L.n.1.redundant = try.redundant$L.n.1; L.n.2.redundant = try.redundant$L.n.2;
#   V.redundant <- rep(0, 2)
#   for(i in 1:2 ){
#     V.redundant[i] = L.n.2.redundant[i] - (L.n.1.redundant[2*i-1] +  L.n.1.redundant[2*i])
#   }
# 
#   #use the maximum value of V.redundant as the reference V value
#   V <- c(V,rep(max(V.redundant),floor(2*length(V))))
#   if( length(unique(V)) <= 2 ){
#     omega <- max(V);
#   }
#   if( length(unique(V)) > 2 ){
#     #use kmeans to cluster the V 
#     clus.2 <- kmeans(V, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss; 
#     if(fit.2 < 0.20){
#       omega <- max(V);
#     }
#     if( fit.2 >= 0.20 ){
#       #if the reference point is in the subset with larger center, this means no change points: set omeage = max(V)
#       #otherwise, set omega = min(V) -1 
#       loc <- clus.2$cluster;
#       if( clus.2$centers[1] > clus.2$centers[2]  ){
#         omega <- min(V[which(loc==1)]) -1 ;
#         if(loc[length(loc)] == 1){
#           omega <- max(V);
#         }
#       }
#       if( clus.2$centers[1] < clus.2$centers[2]  ){
#         omega <- min(V[which(loc==2)]) -1 ;
#         if(loc[length(loc)] == 2){
#           omega <- max(V);
#         }
#       }
#     }
#   }
#   
#   #select the break points by localized information criterion (LIC)
#   L.n.1.temp <- L.n.1; L.n.2.temp <- L.n.2; L.n.plot <- rep(0,m+1); L.n.plot[1] <- sum(L.n.1) + m*omega; 
#   mm <- 0; ic <- 0; add.temp <- 0; pts.full <- vector("list",m+1); pts.full[[1]] <- pts; ind.pts <- rep(0,m);
#   while(mm < m){
#     mm <- mm + 1;
#     L.n.temp <- rep(0,length(pts));
#     for(i in 1:length(pts)){
#       L.n.temp[i] <- sum(L.n.1.temp) - L.n.1.temp[(2*i-1)] - L.n.1.temp[(2*i)] + L.n.2.temp[i] + 1*add.temp;
#     }
#     ll <- min(which.min(L.n.temp)); ind.pts[mm] <- ll;
#     pts <- pts[-ll]; 
#     L.n.1.temp <- L.n.1.temp[-c(2*ll-1,2*ll)]; add.temp <- add.temp + 1*L.n.2.temp[ll]; 
#     L.n.2.temp <- L.n.2.temp[-ll]; 
#     L.n.plot[mm+1] <- L.n.temp[ll] + (m - mm)*omega;
#     pts.full[[mm+1]] <- pts;
#   }
#   
#   ind <- 0;
#   ind <- min(which.min(L.n.plot))
#   
#   return(list(pts = pts.full[[ind]], omega = omega  ))
# }



#' Compute local loss function.
#' 
#' @param data input data matrix, with each column representing the time series component 
#' @param eta tuning parameter eta for lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param pts the selected break points after the first step
#' @param an the neighborhood size a_n
#' @return A list oject, which contains the followings
#' \describe{
#'   \item{L.n.1}{A vector of loss functions that include some break point}
#'   \item{L.n.2}{A vector of loss functions that exclude some break point}
#' }

break.var.local.new <- function(data, eta, q, max.iteration = 1000, tol = 10^(-4),  pts, an){
  p <- length(data[1,]); T <- length(data[,1]); m <- length(pts);
  
  #construct the local interval for computing the loss function
  bounds.1 <- vector("list",2*m); bounds.2 <- vector("list",m);
  for(i in 1:m){
    bounds.1[[(2*i-1)]] <- c(pts[i] - an, pts[i] - 1 );
    bounds.1[[(2*i)]] <- c(pts[i], pts[i] + an );
    bounds.2[[(i)]] <- c(pts[i] - an, pts[i] + an );
  }
  
  #compute the local loss function that include the given break point
  L.n.1 <- c()
  # add a hashtable to store the local estimate
  phi.local.1 <- vector("list", m); 
  phi.local.2 <- vector("list", m); 
  
  for(mm in 1:(2*m)){
    data.temp <- data[(bounds.1[[mm]][1]):(bounds.1[[mm]][2]),];
    try <- var_lasso_brk(data = data.temp, eta, q, 1000, tol = tol)
    L.n.1 <- c(L.n.1 , try$pred.error)
    key <- ceiling((mm)/2)
    if(mm %% 2 == 1){
      phi.local.1[[key]] <- try$phi.hat
    }else{
      phi.local.2[[key]] <- try$phi.hat
    }
    
  }
  
  #compute the local loss function that include the given break point
  L.n.2 <- c()
  for(mm in 1:m){
    data.temp <- data[(bounds.2[[mm]][1]):(bounds.2[[mm]][2]),];
    try <- var_lasso_brk(data = data.temp, eta, q, 1000, tol = tol)
    L.n.2 <- c(L.n.2 , try$pred.error)
  }
  
  return(list(L.n.1 = L.n.1, L.n.2 = L.n.2,
              phi.local.1 = phi.local.1, phi.local.2 = phi.local.2))
}


#' exhaustive search step (third step).
#' 
#' @description Perform the exhaustive search to select the break point for each cluster. 
#' 
#' @param data input data matrix, with each column representing the time series component 
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param pts.list the selected break points clustered by a_n after the second step
#' @param an the neighborhood size a_n
#' @param eta tuning parameter eta for lasso
#' @return A list object, which contains the followings
#' \describe{
#'   \item{pts}{a set of final selected break point after the third exhaustive search step}
#' }
third.step.exhaustive <- function(data, q, max.iteration = 1000, tol = tol,  pts.list, an, eta, refine = TRUE ){
  N <- length(data[,1]); p <- length(data[1,]);
  n <- length(pts.list);  #number of cluster
  final.pts <- rep(0,n);
  pts.list.full <- pts.list
  pts.list.full <- c(1, pts.list.full , N)

  #construct the interval for performing the lasso and computing the loss function
  for(i in 1:n){
    pts.temp <- pts.list.full[[i+1]];
    m <- length(pts.temp);
    if( m <= 1  ) {
      final.pts[i] <- pts.temp;
    }
    if( m > 1  ){
      bounds.1 <- vector("list",2*m);
      lb <- max(pts.list.full[[i]]) + an;
      ub <- min(pts.list.full[[i+2]]) - an;
      for(ii in 1:m){
        bounds.1[[(2*ii-1)]] <- c(pts.temp[ii] - an , pts.temp[ii] - 1 );
        bounds.1[[(2*ii)]] <- c(pts.temp[ii], pts.temp[ii] + an );
      }

      L.n <- c()
      for(mm in 1:(2*m)){
        data.temp <- data[(bounds.1[[mm]][1]):(bounds.1[[mm]][2]),];
        try <- var_lasso_brk(data = data.temp, eta, q, 1000, tol = tol)
        L.n <- c(L.n , try$pred.error)
      }

      sse.full <- rep(0,m)
      for(ii in 1:m ){
        sse.full[ii] = abs(L.n[2*ii-1]  +  L.n[2*ii])
      }

      #select the point that has the smallest SSE among the cluster
      final.pts[i] <- pts.list.full[[i+1]][min(which(sse.full == min(sse.full)))];
    }
  }

  return(pts = final.pts)
}

 
#' local refinement step (forth step).
#' 
#' @description Perform the exhaustive search to select the break point for each cluster. 
#' 
#' @param data input data matrix, with each column representing the time series component 
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso 
#' @param pts.list the selected break points clustered by a_n after the second step
#' @param an the neighborhood size a_n
#' @param phi.est.full list of local parameter estimator
#' @param blocks the blocks
#' @return A list object, which contains the followings
#' \describe{
#'   \item{pts}{a set of final selected break point after the third exhaustive search step}
#' }
forth.step.refine <- function(data, q, max.iteration = 1000, tol = tol, pts.list, 
                              an, phi.est.full = NULL, phi.local.1 = NULL, phi.local.2 = NULL,
                              blocks = NULL ){
  N <- length(data[,1]); p <- length(data[1,]);
  n <- length(pts.list);  #number of cluster
  n.new <- length(phi.est.full)
  final.pts <- rep(0, n);
  pts.list.full <- pts.list
  pts.list.full <- c(1, pts.list.full , N)
  phi.hat.list <- vector("list", n + 1)
  cp.index.list <- vector("list", n + 2);
  cp.index.list[[1]] <- c(1);
  cp.index.list[[n+2]] <- c(n.new+1);

    
  cp.list.full <- vector("list", n+2);
  cp.list.full[[1]] <- c(1);
  cp.list.full[[n+2]] <- c(N+1);
  
  
  for(i in 1:n){
    pts.temp <- pts.list.full[[i+1]];
    m <- length(pts.temp);
    if( m <= 1  ) {
      cp.list.full[[i+1]] <- c((pts.temp-an + 1 ):(pts.temp + an -1) )
      cp.index.list[[i+1]] <- match(pts.temp, blocks)
      
    }
    if( m > 1  ){
      cp.list.full[[i+1]] <- c((pts.temp[1] ):(pts.temp[length(pts.temp)]) )
      cp.index.list[[i+1]] <- match(pts.temp, blocks)
    }
  }
  
  
  #construct the interval for performing the lasso and computing the loss function
  for(i in 1:n){
    idx <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2);
    phi.hat.list[[i]] <- phi.est.full[[idx]]

    
    pts.temp <- pts.list.full[[i+1]];
    m <- length(pts.temp);
    
    #compare the SSE of first num and last num
    num  = cp.list.full[[i+1]][1]
    lb.1 <- min(pts.temp) - an;
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    # idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
    # phi.hat <- phi.est.full[[idx.1]]
    phi.hat <- phi.local.1[[i]]
    forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data), matrix(phi.hat, ncol = p*q), q, jjj-1 , p, 1) )
    if(len.1 == 1){
      temp.1 <- sum( (data[lb.1:ub.1,]-forecast)^2 );
    }else{
      temp.1 <- sum( (t(data[lb.1:ub.1,])-forecast)^2 );
    }
    
    
    lb.2 <- num ;
    ub.2 <- max(pts.temp) +  an -1;
    len.2 <- ub.2 - lb.2 + 1;
    # idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    phi.hat <- phi.local.2[[i]]
    # phi.hat <- phi.est.full[[idx.2]]
    forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data), matrix(phi.hat, ncol = p*q), q, jjj-1 , p, 1) )
    if(len.2 == 1){
      temp.2 <- sum( ( data[lb.2:ub.2,]-forecast)^2 );
    }else{
      temp.2 <- sum( (t(data[lb.2:ub.2,])-forecast)^2 );
    }
    
    sse1 <- temp.1 + temp.2;
    num  <- cp.list.full[[i+1]][length(cp.list.full[[i+1]])]
    lb.1 <- min(pts.temp) - an;
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
    # phi.hat <- phi.est.full[[idx.1]]
    phi.hat <- phi.local.1[[i]]
    forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data), matrix(phi.hat, ncol = p*q), q, jjj-1 , p, 1) )
    if(len.1 == 1){
      temp.1 <- sum( (data[lb.1:ub.1,]-forecast)^2 );
    }else{
      temp.1 <- sum( (t(data[lb.1:ub.1,])-forecast)^2 );
    }
    
    lb.2 <- num ;
    ub.2 <- max(pts.temp) + an -1;
    len.2 <- ub.2 - lb.2 + 1;
    idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    # phi.hat <- phi.est.full[[idx.2]]
    phi.hat <- phi.local.2[[i]]
    forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data), matrix(phi.hat, ncol = p*q), q, jjj-1 , p, 1) )
    if(len.2 == 1){
      temp.2 <- sum( ( data[lb.2:ub.2,]-forecast)^2 );
    }else{
      temp.2 <- sum( (t(data[lb.2:ub.2,])-forecast)^2 );
    }
    sse2 <- temp.1 + temp.2;
    
    
    if(sse1 <= sse2){
      sse.full <- 0;
      ii <- 0
      for(num in cp.list.full[[i+1]]  ){
        ii <- ii + 1
        lb.1 <- min(pts.temp) - an;
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        # phi.hat <- phi.est.full[[idx.1]]
        phi.hat <- phi.local.1[[i]]
        forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data), matrix(phi.hat, ncol = p*q), q, jjj-1 , p, 1) )
        if(len.1 == 1){
          temp.1 <- sum( (data[lb.1:ub.1,]-forecast)^2 );
        }else{
          temp.1 <- sum( (t(data[lb.1:ub.1,])-forecast)^2 );
        }
        
        lb.2 <- num ;
        ub.2 <- max(pts.temp) +  an - 1;
        len.2 <- ub.2 - lb.2 + 1;
        idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
        # phi.hat <- phi.est.full[[idx.2]]
        phi.hat <- phi.local.2[[i]]
        forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data), matrix(phi.hat, ncol = p*q), q, jjj-1 , p, 1) )
        if(len.2 == 1){
          temp.2 <- sum( ( data[lb.2:ub.2,]-forecast)^2 );
        }else{
          temp.2 <- sum( (t(data[lb.2:ub.2,])-forecast)^2 );
        }
        sse.full[ii] <- temp.1 + temp.2;
        # print(ii)
        # print(sse.full[ii])
        if(ii >= min(round(3/2*an), length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full, 0.20) ){
          break
        }
      }
      #select the point that has the smallest SSE among the cluster
      final.pts[i] <- cp.list.full[[i+1]][min(which(sse.full == min(sse.full)))];
    }
    if(sse1 > sse2){
      sse.full <- 0;
      ii <- 0
      for(num in rev(cp.list.full[[i+1]])  ){
        ii <- ii + 1
        lb.1 <- min(pts.temp) - an;
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        # phi.hat <- phi.est.full[[idx.1]]
        phi.hat <- phi.local.1[[i]]
        forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data), matrix(phi.hat, ncol = p*q), q, jjj-1 , p, 1) )
        if(len.1 == 1){
          temp.1 <- sum( (data[lb.1:ub.1,]-forecast)^2 );
          
        }else{
          temp.1 <- sum( (t(data[lb.1:ub.1,])-forecast)^2 );
        }
        
        
        lb.2 <- num ;
        ub.2 <- max(pts.temp) +  an -1;
        len.2 <- ub.2 - lb.2 + 1;
        idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
        # phi.hat <- phi.est.full[[idx.2]]
        phi.hat <- phi.local.2[[i]]
        forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data), matrix(phi.hat, ncol = p*q), q, jjj-1 , p, 1) )
        if(len.2 == 1){
          temp.2 <- sum( (data[lb.2:ub.2,]-forecast)^2 );
        }else{
          temp.2 <- sum( (t(data[lb.2:ub.2,])-forecast)^2 );
        }
        sse.full[ii] <- temp.1 + temp.2;
        if(ii >= min(round(3/2*an), length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full,0.20) ){
          break
        }
      }
      #select the point that has the smallest SSE among the cluster
      final.pts[i] <- cp.list.full[[i+1]][length(cp.list.full[[i+1]]) + 1 - min(which(sse.full == min(sse.full)))];
      
    }
    
    
    
  }
  

  idx <- floor((min(cp.index.list[[n+2]]) + max(cp.index.list[[n+1]]))/2);
  phi.hat.list [[n+1]] <- phi.est.full[[idx]]
  return( list(pts = final.pts, phi.hat.list = phi.hat.list ))
}




#' cluster the points by neighborhood size a_n
block.finder <- function(pts,an){
  nn <- length(pts);
  if( nn == 1){b <- pts;}
  if( nn > 1){
    b <- vector("list",nn);
    i.ind <- 1;
    jj <- 0;
    while (i.ind < nn) {
      ct <- 1;
      jj <- jj + 1;
      for (j in (i.ind+1):nn) {
        if( abs(pts[i.ind] - pts[j]  ) <= an   ){ct <- ct + 1;}
      }
      b[[jj]] <- pts[(i.ind):(i.ind+ct-1)];
      i.ind <- i.ind + ct;
    }
    l <- length(b[[jj]]);
    if(b[[jj]][l] != pts[nn]  ){
      jj <- jj + 1;
      b[[(jj)]] <- c(pts[nn])   
    }
    b <- b[(1):(jj)];
  }
  
  return(b = b)
}





#' Prediction function 1
pred.block <- function(Y,phi,q,T,p,h){
  concat.Y <- matrix(0,p,q+h); concat.Y[,1:q] <- Y[,(T-q+1):T];
  for ( j in 1:h){
    temp <- matrix(0,p,1);
    for (i in 1:q){temp <- temp +  phi[,((i-1)*p+1):(i*p)]%*%concat.Y[,q+j-i];}
    concat.Y[,q+j] <- temp; 
  }
  return(as.matrix(concat.Y[,(q+1):(q+h)]))
}

#' Prediction function 2
pred <- function(Y,phi,q,T,p,h){
  concat.Y <- matrix(0,p,q+h); concat.Y[,1:q] <- Y[,(T-q+1):T];
  for ( j in 1:h){
    temp <- matrix(0,p,1);
    for (i in 1:q){temp <- temp +  phi[,((i-1)*p+1):(i*p)]%*%concat.Y[,q+j-i];}
    concat.Y[,q+j] <- temp; 
  }
  return(as.matrix(concat.Y[,q+h]))
}

#' Plot the AR coefficient matrix
plot.matrix <- function (phi,p,name = NULL) {
  B <- phi
  if (nrow(B) == 1) {
    B <- matrix(B[, 1:ncol(B)], nrow = 1)
  }
  else {
    B <- B[, 1:ncol(B)]
  }
  k <- nrow(B)
  s1 <- 0
  m <- 0
  s <- 0
  s <- s + s1
  text <- c()
  for (i in 1:p) {
    text1 <- as.expression(bquote(bold(Phi)^(.(i))))
    text <- append(text, text1)
  }
  if (m > 0) {
    for (i in (p + 1):(p + s + 1)) {
      text1 <- as.expression(bquote(bold(beta)^(.(i - p - 
                                                    s1))))
      text <- append(text, text1)
    }
  }
  f <- function(m) t(m)[, nrow(m):1]
  rgb.palette <- colorRampPalette(c("blue", "white", "red"), space = "Lab")
  at <- seq(k/2 + 0.5, p * (k) + 0.5, by = k)
  if (m > 0) {
    at2 <- seq(p * k + m/2 + 0.5, p * k + s * m + 0.5, by = m)
  }
  else {
    at2 = c()
  }
  at <- c(at, at2)
  se2 = seq(1.75, by = k, length = k)
  L2 <- levelplot(as.matrix(f(B)), at =seq( -1, 1, length=101),col.regions = rgb.palette(100), 
                  colorkey = NULL, xlab = NULL, ylab = NULL, main = list(label = name, 
                                                                         cex = 1), panel = function(...) {
                                                                           panel.levelplot(...)
                                                                           panel.abline(a = NULL, b = 1, h = seq(1.5, m * s + 
                                                                                                                   p * k + 0.5, by = 1), v = seq(1.5, by = 1, length = p * 
                                                                                                                                                   k + m * s), lwd = 0.5)
                                                                           bl1 <- seq(k + 0.5, p * k + 0.5, by = k)
                                                                           b23 <- seq(p * k + 0.5, p * k + 0.5 + s * m, by = m)
                                                                           b1 <- c(bl1, b23)
                                                                           panel.abline(a = NULL, b = 1, v = p * k + 0.5, lwd = 3)
                                                                           panel.abline(a = NULL, b = 1, v = b1, lwd = 2)
                                                                         }, scales = list(x = list(alternating = 1, labels = text, 
                                                                                                   cex = 2, at = at, tck = c(0, 0)), y = list(alternating = 0, 
                                                                                                                                              tck = c(0, 0))))
  return(L2)
}




remove.extra.pts <- function(pts, brk){
  m.hat <- length(brk)-1;
  if(length(pts) <= m.hat){break;}
  pts.temp <- rep(0, m.hat);
  for(i in 1:m.hat){
    origin <- brk[i];
    dis <- rep(0, length(pts));
    for(j in 1:length(pts)){
      dis[j] <- abs(origin - pts[j]);
    }
    ll <- min(which.min(dis));
    pts.temp[i] <- pts[ll];
  }
  
  pts <- pts.temp;
  return(pts)
}



detection.check <- function(pts.final, brk, N){
  m <- length(brk); len <- rep(0, N);
  for(i in 1:N){len[i] <- length(pts.final[[i]]);}
  freq <- as.matrix(table(len)/N);
  pts.final.full.1 <- vector("list", N); pts.final.full.2 <- vector("list", N);
  for(i in 1:N){
    if ( length(pts.final[[i]]) > (m-1)   ){
      pts.final[[i]] <- remove.extra.pts(pts.final[[i]], brk);
      print("extra points exists!"); 
      print(i)
    }
    if ( length(pts.final[[i]]) == 0 ) {  
      pts.final[[i]] <- rep(0, m-1);
      # pts.final.full.1 only counts those points within 1/5 intervals in order to compute the selection rate
      # pts.final.full.2 counts all the points in order to calculate the mean and sd for the location error
      pts.final.full.1[[i]] <- rep(0, m-1); 
      pts.final.full.2[[i]] <- rep(0, m-1);  
    }
    if ( length(pts.final[[i]]) > 0 && length(pts.final[[i]]) <= (m-1) ){
      ll <- length(pts.final[[i]]); 
      pts.final.full.1[[i]] <- rep(0, m-1); pts.final.full.2[[i]] <- rep(0, m-1);
      for(j in 1:ll){
        if(m == 2){
          if (  pts.final[[i]][j] < (brk[1] + (1/3)*(brk[2] - brk[1]) ) && pts.final[[i]][j] >= (brk[1] - (1/3)*(brk[1]) )   ) {
            pts.final.full.1[[i]][1] <- pts.final[[i]][j];
          }
          if (  pts.final[[i]][j] < (brk[1] + (1/2)*(brk[2] - brk[1]) )   ) {
            pts.final.full.2[[i]][1] <- pts.final[[i]][j];
          }
        }else{
          if (  pts.final[[i]][j] < (brk[1] + (1/5)*(brk[2] - brk[1]) ) && pts.final[[i]][j] >= (brk[1] - (1/5)*(brk[1]) )   ) {
            pts.final.full.1[[i]][1] <- pts.final[[i]][j];
          }
          if (  pts.final[[i]][j] < (brk[1] + (1/2)*(brk[2] - brk[1]) )   ) {
            pts.final.full.2[[i]][1] <- pts.final[[i]][j];
          }
          
          for(kk in 2:(m-1)){
            if (  pts.final[[i]][j] >=  (brk[(kk)] - (1/5)*(brk[kk] - brk[(kk-1)]) ) && pts.final[[i]][j] <  (brk[(kk)] + (1/5)*(brk[kk+1] - brk[(kk)]) )){
              pts.final.full.1[[i]][kk] <- pts.final[[i]][j];
            }
            if(kk == m-1){
              if (  pts.final[[i]][j] >=  (brk[(kk)] - (1/2)*(brk[kk] - brk[(kk-1)]) ) ){
                pts.final.full.2[[i]][kk] <- pts.final[[i]][j];
              }
              
            }else{
              if (  pts.final[[i]][j] >=  (brk[(kk)] - (1/2)*(brk[kk] - brk[(kk-1)]) ) && pts.final[[i]][j] <  (brk[(kk)] + (1/2)*(brk[kk+1] - brk[(kk)]) )){
                pts.final.full.2[[i]][kk] <- pts.final[[i]][j];
              }
            }
            
          }
          
        }
        
        
      }
    }
  }
  
  detection <- matrix(0, m, 7)
  # print(pts.final)
  # print(pts.final.full.1)
  # print(pts.final.full.2)
  
  detection[1, 1] <- c("break points"); detection[1, 2] <- c("truth"); 
  detection[1, 3] <- c("mean(errors)"); detection[1, 4] <- c("std(errors)"); 
  detection[1, 5] <- c("selection rate");
  detection[1, 6] <- c("mean(location)"); detection[1,7] <- c("std(location)");
  for(i in 1:(m-1)){
    detection[(i+1),1] <- c(i); detection[(i+1), 2] <- c(brk[i]); 
    # loc <- rep(0, N); 
    loc.1 <- rep(0, N); loc.2 <- rep(0, N);
    for(j in 1:N){
      # temp <- pts.final[[j]]; l <- length(temp); loc[j] <- temp[i];
      temp.1 <- pts.final.full.1[[j]];loc.1[j] <- temp.1[i];
      temp.2 <- pts.final.full.2[[j]];loc.2[j] <- temp.2[i];
    }
    # loc <-  loc[which(loc!=0)]
    loc.1 <- loc.1[which(loc.1!=0)]; loc.2 <- loc.2[which(loc.2!=0)];
    T.new.1 <- length(loc.1); 
    detection[(i+1),3] <- mean(abs(loc.2-brk[i])); detection[(i+1),4] <- sd(abs(loc.2-brk[i]));
    detection[(i+1),5] <- T.new.1/N;
    detection[(i+1),6] <- mean(loc.2/T); detection[(i+1),7] <- sd(loc.2/T);
  }
  
  for(i in 2:(m)){
    for(j in 2:7){
      detection[i,j] <- round(as.numeric(detection[i,j]), digits = 4)
    }
  }
  return(list(pts.final.full.1 = pts.final.full.1, pts.final.full.2 = pts.final.full.2,  detection = detection, freq = freq))
}





hausdorff.check <- function(pts.final,  brk, N){
  m <- length(brk); len <- rep(0,N);
  for(i in 1:N){len[i] <- length(pts.final[[i]]);}
  
  #hausdorff_true_est d(A_n, \tilde{A}_n^f)
  #for each estimated cp, find the closest true CP and compute the distance,
  # then take the maximum of distances
  pts.final.full.1 <- vector("list",N);
  for(i in 1:N){
    ll <- length(pts.final[[i]]); pts.final.full.1[[i]] <- rep(NA,ll);
    if(ll>0){
      for(j in 1:ll){
        if (   pts.final[[i]][j] < (brk[1] + (1/2)*(brk[2] - brk[1]) )    ) {pts.final.full.1[[i]][j] <- brk[1];}
        if (m -2 > 0){
          if (   pts.final[[i]][j] >= (brk[m-2] + (1/2)*(brk[m-1] - brk[m-2]) )    ) {pts.final.full.1[[i]][j] <- brk[m-1];}
        }
        if(m-2 >=2){
          for(kk in 2:(m-2)){
            if (  pts.final[[i]][j] >=  (brk[(kk-1)] + (1/2)*(brk[kk] - brk[(kk-1)])) && pts.final[[i]][j] <  (brk[(kk)] + (1/2)*(brk[kk+1] - brk[(kk)]) )){
              pts.final.full.1[[i]][j] <- brk[kk];
            }
          }
        }
      }
    }else{
      #print("zero estimated CP!")
      #print(i)
    }
  }
  
  
  #hausdorff_est_true d(\tilde{A}_n^f, A_n)
  #for each true cp, find the closest estimate CP and compute the distance,
  # then take the maximum of distances
  pts.final.full.2<- vector("list",N);
  for(i in 1:N){
    ll <- length(pts.final[[i]]); pts.final.full.2[[i]] <- rep(NA, m -1 );
    if(ll > 0){
      if(ll > 1){
        for(j in 1: (m-1)){
          if (   brk[j] < (pts.final[[i]][1] + (1/2)*(pts.final[[i]][2] - pts.final[[i]][1]) )    ) {pts.final.full.2[[i]][j] <- pts.final[[i]][1];}
          if (   brk[j] >= (pts.final[[i]][ll-1] + (1/2)*(pts.final[[i]][ll] - pts.final[[i]][ll-1]) )    ) {pts.final.full.2[[i]][j] <- pts.final[[i]][ll];}
          if(ll >= 3){
            for(kk in 2:(ll-1)){
              if (  brk[j] >=  (pts.final[[i]][(kk-1)] + (1/2)*(pts.final[[i]][kk] - pts.final[[i]][(kk-1)])) 
                    && brk[j] <  (pts.final[[i]][(kk)] + (1/2)*(pts.final[[i]][kk+1] - pts.final[[i]][(kk)]) )){
                pts.final.full.2[[i]][j] <- pts.final[[i]][kk];
              }
            }
          }
        }
        
      }else{
        pts.final.full.2[[i]] <- rep(pts.final[[i]], m -1 );
      }
      
    }
  }
  
  
  detection <- matrix(0, 2, 6)
  
  detection[1,1] <- c("mean(hausdorff_true_est)"); detection[1,2] <- c("std(hausdorff_true_est)"); 
  detection[1,3] <- c("mean(hausdorff_est_true)"); detection[1,4] <- c("std(hausdorff_est_true)");
  detection[1,5] <- c("median(hausdorff_true_est)");
  detection[1,6] <- c("median(hausdorff_est_true)");
  
  distance.1 <- rep(NA, N); distance.2 <- rep(NA, N)
  for(j in 1:N){
    temp.1 <- pts.final.full.1[[j]];
    if(length(temp.1) > 0){
      distance.1[j] <- max(abs(pts.final.full.1[[j]] - pts.final[[j]] ))
    }
    temp.2 <- pts.final.full.2[[j]];
    distance.2[j] <- max(abs(pts.final.full.2[[j]] - brk[1:(m-1)] ))
  }
  detection[2,1] <- mean(distance.1,na.rm= TRUE); detection[2,2] <- sd(distance.1,na.rm= TRUE);
  detection[2,3] <- mean(distance.2,na.rm= TRUE); detection[2,4] <- sd(distance.2,na.rm= TRUE);
  detection[2,5] <- median(distance.1,na.rm= TRUE); 
  detection[2,6] <- median(distance.2,na.rm= TRUE); 
  
  for(j in 1:4){
    detection[2,j] <- round(as.numeric(detection[2,j]),digits = 4)
  }
  return(list(pts.final.full.1 = pts.final.full.1, pts.final.full.2 = pts.final.full.2, detection = detection))
}




