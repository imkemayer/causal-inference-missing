
pscores.fun <- function(treat=Z, outs=Y, covs=X){
#
N <- nrow(covs)
nouts <- 1 
ncovs <- ncol(covs)
#
# first set up places to store results
res <- matrix(0,nouts,2)
bal <- matrix(0,ncovs,2)
#
# estimate p-scores
	dat <- cbind.data.frame(treat=treat,covs)
	mod <- glm(dat,fam=binomial)
	qx <<- mod$linear	
#
### Now Matching With Replacement
	matchout <- matching.wr.x3(treat, qx)
#
### and treatment effect estimation with robust s.e.'s
	wts <- rep(1, N)
	wts[treat == 0] <- matchout$cnts
	res <- wls.all2(cbind(rep(1, sum(wts > 0)), treat[wts > 0],covs[wts > 0,  ]), wts[wts > 0], outs[wts > 0], treat[wts > 0])
c(res[3],sqrt(res[2]))
}

wls.all2 <- function(X, w = wts, Y = y, treat = Trt)
{
	#
	# This produces coefficient estimates and both standard and robust variances 
	# estimates for regression with weights
	# the standard variance corresponds to a situation where an observation represents
	# the mean of w observations
	# the robust variance corresponds to a situation where weights represent 
	# probability or sampling weights
	#
	# first put together the necessary data inputs
	#
	nunits <-  sum(w > 0)
	k <-  ncol(X)
	## now the weights, properly normed
	wn <-  w * (nunits/sum(w))
	W <-  diag(wn * (nunits/sum(wn)))
	#
	# x prime x inverse (including weights)
	vhat <-   - sweep.inv((t(X) %*% W %*% X))
	#
	# estimated regression coefficients and variance for just the treatment coefficient
	b <-  vhat %*% t(X) %*% W %*% Y
	MSE <-  c(t(Y) %*% W %*% Y - t(b) %*% t(X) %*% W %*% Y)/(nunits - k)
	var.std <-  (vhat * MSE)[2, 2]
	#
	######  now for the robust variance calculations
	# now a matrix where each row represents the contribution to the score
	# for each observation
	U <-  c((Y - X %*% b) * wn) * X
	# finite sample adjustment
	qc <-  nunits/(nunits - 2)
	# the sum of outer products of each of the above score contributions for
	# each person is calculated here
	prodU <-  array(0, c(k, k, nunits))
	for(i in 1:nunits) {
		prodU[,  , i] <-  outer(U[i,  ], U[i,  ])
	}
	# putting it all together...
	Vrob <-  qc * vhat %*% apply(prodU, c(1, 2), sum) %*% vhat
	# and we pull off the variance just for the treatment effect 
	var.rob <-  Vrob[2, 2]
	###############
	results <-  c(var.std, var.rob, b[2])
	results
}

matching.wr.x3 <- function(z, score){
	nt <- sum(z)
	nc <- length(z) - nt
	cnts <- numeric(nc)
	scorec <- score[z == 0]
	scoret <- score[z == 1]
	indc <- c(99999)
	nearest <- numeric(nt)
	ind.mt <- matrix(0, nc, nt)
	ind.t <- (1:(nt + nc))[z == 1]
	for(j in 1:nt) {
		near <- (1:nc)[abs(scoret[j] - scorec) == min(abs(scoret[j] -
			scorec))]
		if(length(near) == 1) {
			nearest[j] <- near
			indc <- c(indc, near)
		}
		else {
			nearest[j] <- near[sample(1:length(near), 1, replace = 
				F)]
			indc <- c(indc, nearest[j])
		}
		cnts[nearest[j]] <- cnts[nearest[j]] + 1
		ind.mt[nearest[j], cnts[nearest[j]]] <- ind.t[j]
	}
	#
	ind.mt <- ind.mt[ind.mt[, 1] != 0, 1:max(cnts)]
	# now create list of indicators to pull off appropriate dataset
	ind <- numeric(nt + sum(cnts))
	# first get treat indicators
	ind[1:nt] <- (1:(nt + nc))[z == 1]
	# 
	#now the control indicators
	tmp <- (1:(nt + nc))[z == 0]
	ind[(nt + 1):length(ind)] <- tmp[indc[-1]]
	#
	out <- list(cnts = cnts, ind = ind, ind.mt = ind.mt)
	out
}

sweep.inv <- function(G){
	# sweeps a symmetric matrix on all positions
	# (so inverts the matrix)
	for(i in 1:nrow(G)) {
		G <- sweep.oper(G, i)
	}
	G
}

sweep.oper <- function(G = theta, k = 1.){
	# k is the sweep position
	p <- dim(G)[1.]
	H <- G
	#first do generic elements (those that don't involve k)
	H[] <- 0.
	tmp <- matrix(G[, k], p, 1.) %*% matrix(G[, k], 1., p)
	#now replace the row and col with index=k 
	H <- G - tmp/G[k, k]
	H[, k] <- G[, k]/G[k, k]
	#now replace the (k,k) diagonal element 
	H[k,  ] <- G[, k]/G[k, k]
	# and we're done
	H[k, k] <- -1./G[k, k]
	H
}

matching <- function(z, score){
	# argument z is the vector of indicators for treatment or control #
	# argument score is the vector of the propensity scores in the    #
	# same order as z                                                 #
	# the function returns a vector of indices that the corresponding #
	# unit is matched to. 0 means matched to nothing.                 #
	#                                                                 #
	# now also returns a number for each pair making it easier to     #
	# later recover those pairs
	n <- length(score)
	matched <- rep(0., n)
	pairs <- rep(0., n)
#	if(length(z) != n) print("Error: unequal lengths")
	b <- (sum(z) < n/2.) * 1
	tally <- 0
	for(i in (1:n)[z == b]) {
		available <- (1:n)[(z != b) & (matched == 0)]
		j <- available[order(abs(score[available] - score[i]))[1.]]
		matched[i] <- j
		matched[j] <- i
		tally <- tally + 1.
		pairs[c(i, j)] <- tally
	}
	out <- cbind.data.frame(matched = matched, pairs = pairs)
	out
}

balance.sum <- function(mat.out=NULL,matched=TRUE,num.cont=9,num.cov=41){
#
# this function runs the MatchIt balance function and then produces 
# summaries of the output for the matched group
# different summaries are calculated for the binary variables and the continuous
# (categorical are expected to be input as separate indicator variables)

# since I can't change what statistics are output by the summary function
# (i can add but not subtract variables) i will just have to be clever about
# how i extract them
#
n.cont=num.cont+1
k=num.cov+1
#
bal=summary(mat.out,interactions=FALSE,standardize=TRUE)
if(matched==TRUE){
BB = bal$sum.matched
}
else{
BB = bal$sum.all
}
#
## include the pscore as one of the continuous measures from the others
## 
# now stats for abs standardized diff. in means for the binary variables
# now n.cont doesn't include the "distance variable"
tmp = abs(BB[(n.cont+1):k,"Std. Mean Diff."])
cat("num. binary is",length(tmp),"\n")
tmp = fix(tmp)
bin.std.mn = mean(tmp)
bin.std.max = max(tmp)
bin.std.over.1 = sum(tmp>.1)
#
# now stats for continuous variables (no interactions)
# abs std. diff in means
tmp = abs(BB[1:n.cont,"Std. Mean Diff."])
cat("num. cont. is",length(tmp),"\n")
tmp = fix(tmp)
cont.std.mn = mean(tmp)
cont.std.max = max(tmp)
cont.std.over.1 = sum(tmp>.1)
# abs median diff. in QQ plots
tmp = abs(BB[1:n.cont,"eCDF Med"])
tmp = fix(tmp)
cont.medQQ.mn = mean(tmp)
cont.medQQ.max = max(tmp)
# abs max diff. in QQ plots
tmp = abs(BB[1:n.cont,"eCDF Max"])
tmp = fix(tmp)
cont.maxQQ.mn = mean(tmp)
cont.maxQQ.max = max(tmp)

# now stats for squared terms and interaction terms (in some sense this
# double counts the importance of the second moment but i'm not sure what
# to do about) for continuous variables with each other and continuous
# interacted with binary
stop.row = nrow(BB)
# abs std. diff in means
tmp = abs(BB[(k+1):stop.row,"Std. Mean Diff."])
cat("num. interacts is",length(tmp),"\n")
tmp = fix(tmp)
cont2.std.mn = mean(tmp)
cont2.std.max = max(tmp)
cont2.std.over.1 = sum(tmp>.1)
# abs median diff. in QQ plots
tmp = abs(BB[(k+1):stop.row,"eCDF Med"])
tmp = fix(tmp)
cont2.medQQ.mn = mean(tmp)
cont2.medQQ.max = max(tmp)
# abs max diff. in QQ plots
tmp = abs(BB[(k+1):stop.row,"eCDF Max"])
tmp = fix(tmp)
cont2.maxQQ.mn = mean(tmp)
cont2.maxQQ.max = max(tmp)

#
bal.bin = list(bin.std.mn=bin.std.mn,bin.std.max=bin.std.max,bin.std.over.1=bin.std.over.1)
#
bal.cont1 = list(cont.std.mn = cont.std.mn, cont.std.max = cont.std.max, cont.std.over.1 = cont.std.over.1,cont.medQQ.mn = cont.medQQ.mn, cont.medQQ.max = cont.medQQ.max, cont.maxQQ.mn = cont.maxQQ.mn,cont.maxQQ.max = cont.maxQQ.max)

bal.cont2 = list(cont2.std.mn = cont2.std.mn, cont2.std.max = cont2.std.max, cont2.std.over.1 = cont2.std.over.1,cont2.medQQ.mn = cont2.medQQ.mn, cont2.medQQ.max = cont2.medQQ.max, cont2.maxQQ.mn = cont2.maxQQ.mn,cont2.maxQQ.max = cont2.maxQQ.max)
#
bal.sum = list(bal.bin=bal.bin,bal.cont1=bal.cont1,bal.cont2=bal.cont2)
return(bal.sum)
}

balance.gm2.sum <- function(mat.out=NULL,matched=TRUE,dat=usek2,num.cont=5,num.cov=27,form.q=form.quad,
                           verbose=TRUE){
#
# this function runs the GenMatch balance function and then produces 
# summaries of the output for the matched group
# different summaries are calculated for the binary variables and the continuous
# (categorical are expected to be input as separate indicator variables)

# since I can't change what statistics are output by the summary function
# (i can add but not subtract variables) i will just have to be clever about
# how i extract them
#
k=num.cov+1
n.cont=num.cont+1
#
bal=MatchBalance(formul=form.q,data=dat,match.out=mat.out,ks=FALSE,print.level=0)
if(matched==TRUE){
BB = bal$AfterMatching
}
else{
BB = bal$BeforeMatching
}
kk=length(BB)
#
## include the pscore as one of the continuous measures from the others
## 
std.diffs = unlist(BB)[(0:(kk-1))*18+1]/100
mn.QQs = unlist(BB)[(0:(kk-1))*18+14]
max.QQs = unlist(BB)[(0:(kk-1))*18+15]

# now stats for abs standardized diff. in means for the binary variables
tmp = abs(std.diffs[(n.cont+1):k])
if(verbose)
  cat("num. binary is",length(tmp),"\n")
tmp = fix(tmp)
bin.std.mn = mean(tmp)
bin.std.max = max(tmp)
bin.std.over.1 = sum(tmp>.1)
#
# now stats for continuous variables (no interactions)
# abs std. diff in means
tmp = abs(std.diffs[1:n.cont])
if(verbose)
  cat("num. cont. is",length(tmp),"\n")
tmp = fix(tmp)
cont.std.mn = mean(tmp)
cont.std.max = max(tmp)
cont.std.over.1 = sum(tmp>.1)
# abs median diff. in QQ plots
tmp = abs(mn.QQs[1:n.cont])
tmp = fix(tmp)
cont.medQQ.mn = mean(tmp)
cont.medQQ.max = max(tmp)
# abs max diff. in QQ plots
tmp = abs(max.QQs[1:n.cont])
tmp = fix(tmp)
cont.maxQQ.mn = mean(tmp)
cont.maxQQ.max = max(tmp)

# now stats for squared terms and interaction terms (in some sense this
# double counts the importance of the second moment but i'm not sure what
# to do about) for continuous variables with each other and continuous
# interacted with binary
# abs std. diff in means
tmp = abs(std.diffs[(k+1):kk])
if(verbose)
  cat("num. interacts is",length(tmp),"\n")
tmp = fix(tmp)
cont2.std.mn = mean(tmp)
cont2.std.max = max(tmp)
cont2.std.over.1 = sum(tmp>.1)
# abs median diff. in QQ plots
tmp = abs(mn.QQs[(k+1):kk])
tmp = fix(tmp)
cont2.medQQ.mn = mean(tmp)
cont2.medQQ.max = max(tmp)
# abs max diff. in QQ plots
tmp = abs(max.QQs[(k+1):kk])
tmp = fix(tmp)
cont2.maxQQ.mn = mean(tmp)
cont2.maxQQ.max = max(tmp)

#
bal.bin = list(bin.std.mn=bin.std.mn,bin.std.max=bin.std.max,bin.std.over.1=bin.std.over.1)
#
bal.cont1 = list(cont.std.mn = cont.std.mn, cont.std.max = cont.std.max, cont.std.over.1 = cont.std.over.1,cont.medQQ.mn = cont.medQQ.mn, cont.medQQ.max = cont.medQQ.max, cont.maxQQ.mn = cont.maxQQ.mn,cont.maxQQ.max = cont.maxQQ.max)

bal.cont2 = list(cont2.std.mn = cont2.std.mn, cont2.std.max = cont2.std.max, cont2.std.over.1 = cont2.std.over.1,cont2.medQQ.mn = cont2.medQQ.mn, cont2.medQQ.max = cont2.medQQ.max, cont2.maxQQ.mn = cont2.maxQQ.mn,cont2.maxQQ.max = cont2.maxQQ.max)
#
bal.sum = list(bal.bin=bal.bin,bal.cont1=bal.cont1,bal.cont2=bal.cont2)
return(bal.sum)
}

balance.gm3.sum <- function(mat.out=NULL,matched=TRUE,dat=usek2,num.cont=5,num.cov=27,form.q=form.quad,
                           verbose=TRUE, xdata=NULL, Tr=NULL){
#
# this function runs the GenMatch balance function and then produces 
# summaries of the output for the matched group
# different summaries are calculated for the binary variables and the continuous
# (categorical are expected to be input as separate indicator variables)

# since I can't change what statistics are output by the summary function
# (i can add but not subtract variables) i will just have to be clever about
# how i extract them
#
k=num.cov+1
n.cont=num.cont+1
#

if(verbose==TRUE)
  {
    bal=MatchBalance(formul=form.q,data=dat,match.out=mat.out,ks=FALSE,print.level=2)
  } else{
    bal <- MatchBalance1(formul=form.q,data=dat,match.out=mat.out,ks=FALSE,print.level=0,
                         xdata=xdata, Tr=Tr)    
  }

if(matched==TRUE){
BB = bal$AfterMatching
}
else{
BB = bal$BeforeMatching
}
kk=length(BB)
#
## include the pscore as one of the continuous measures from the others
## 
std.diffs = unlist(BB)[(0:(kk-1))*18+1]/100
mn.QQs = unlist(BB)[(0:(kk-1))*18+14]
max.QQs = unlist(BB)[(0:(kk-1))*18+15]

# now stats for abs standardized diff. in means for the binary variables
tmp = abs(std.diffs[(n.cont+1):k])
if(verbose)
  cat("num. binary is",length(tmp),"\n")
tmp = fix(tmp)
bin.std.mn = mean(tmp)
bin.std.max = max(tmp)
bin.std.over.1 = sum(tmp>.1)
#
# now stats for continuous variables (no interactions)
# abs std. diff in means
tmp = abs(std.diffs[1:n.cont])
if(verbose)
  cat("num. cont. is",length(tmp),"\n")
tmp = fix(tmp)
cont.std.mn = mean(tmp)
cont.std.max = max(tmp)
cont.std.over.1 = sum(tmp>.1)
# abs median diff. in QQ plots
tmp = abs(mn.QQs[1:n.cont])
tmp = fix(tmp)
cont.medQQ.mn = mean(tmp)
cont.medQQ.max = max(tmp)
# abs max diff. in QQ plots
tmp = abs(max.QQs[1:n.cont])
tmp = fix(tmp)
cont.maxQQ.mn = mean(tmp)
cont.maxQQ.max = max(tmp)

# now stats for squared terms and interaction terms (in some sense this
# double counts the importance of the second moment but i'm not sure what
# to do about) for continuous variables with each other and continuous
# interacted with binary
# abs std. diff in means
tmp = abs(std.diffs[(k+1):kk])
if(verbose)
  cat("num. interacts is",length(tmp),"\n")
tmp = fix(tmp)
cont2.std.mn = mean(tmp)
cont2.std.max = max(tmp)
cont2.std.over.1 = sum(tmp>.1)
# abs median diff. in QQ plots
tmp = abs(mn.QQs[(k+1):kk])
tmp = fix(tmp)
cont2.medQQ.mn = mean(tmp)
cont2.medQQ.max = max(tmp)
# abs max diff. in QQ plots
tmp = abs(max.QQs[(k+1):kk])
tmp = fix(tmp)
cont2.maxQQ.mn = mean(tmp)
cont2.maxQQ.max = max(tmp)

#
bal.bin = list(bin.std.mn=bin.std.mn,bin.std.max=bin.std.max,bin.std.over.1=bin.std.over.1)
#
bal.cont1 = list(cont.std.mn = cont.std.mn, cont.std.max = cont.std.max, cont.std.over.1 = cont.std.over.1,cont.medQQ.mn = cont.medQQ.mn, cont.medQQ.max = cont.medQQ.max, cont.maxQQ.mn = cont.maxQQ.mn,cont.maxQQ.max = cont.maxQQ.max)

bal.cont2 = list(cont2.std.mn = cont2.std.mn, cont2.std.max = cont2.std.max, cont2.std.over.1 = cont2.std.over.1,cont2.medQQ.mn = cont2.medQQ.mn, cont2.medQQ.max = cont2.medQQ.max, cont2.maxQQ.mn = cont2.maxQQ.mn,cont2.maxQQ.max = cont2.maxQQ.max)
#
bal.sum = list(bal.bin=bal.bin,bal.cont1=bal.cont1,bal.cont2=bal.cont2)
return(bal.sum)
}


fix <- function(a){
a[a==Inf] = NA
a = na.omit(a)
}

MatchBalance1 <- function(formul, data=NULL, match.out=NULL, ks=FALSE, 
                          nboots=0, weights=NULL,
                          digits=5, paired=TRUE, print.level=1,
                          xdata=NULL, Tr=NULL)
  {
    nboots <- 0
    if(!is.list(match.out) & !is.null(match.out)) {
      warning("'Match' object contains no valid matches")
      match.out  <- NULL
    }

    if ( (class(match.out) != "Match") & (class(match.out) != "Matchby") & (!is.null(match.out)) ) {
      warning("Object not of class 'Match'")
      match.out  <- NULL
    }


    orig.na.action <- as.character(options("na.action"))
    options("na.action"=na.pass)

    if (is.null(data))
      {
        if(is.null(xdata)) {
          xdata <- as.data.frame(get.xdata(formul,datafr=environment(formul)))
        }

        if(is.null(Tr)) {
          Tr <- as.real(get.ydata(formul,datafr=environment(formul)))
        }

      } else {

        if(is.null(xdata) | is.null(Tr)) {
          data  <- as.data.frame(data)
        }
        
        if(is.null(xdata)) {
          xdata  <- as.data.frame(get.xdata(formul, data))
        }
        if(is.null(Tr)) {        
          Tr  <- as.real(get.ydata(formul, data))
        }
      }
    options("na.action"=orig.na.action)

    if (is.null(weights))
      weights <- rep(1,length(Tr))

    if(!is.numeric(weights))
      stop("'weights' must be a numeric vector")

    if( sum(is.na(xdata))!=0 | sum(is.na(Tr))!=0)
      {

        if(orig.na.action!="na.omit" & orig.na.action!="na.exclude" & orig.na.action!="na.fail")
          warning("'na.action' should be set to 'na.omit', 'na.exclude' or 'na.fail' see 'help(na.action)'")

        if (orig.na.action=="na.fail")
          {
            stop("NA's found in data input.")            
          } else {
            warning("NA's found in data input.  IT IS HIGHLY RECOMMENDED THAT YOU TEST IF THE MISSING VALUES ARE BALANCED INSTEAD OF JUST DELETING THEM.")
          }

        indx1 <- apply(is.na(xdata),1,sum)==0 & is.na(Tr)==0
        Tr <- Tr[indx1]
        xdata = xdata[indx1,]
        weights <- weights[indx1]
      } #end of na

    if (sum(Tr !=1 & Tr !=0) > 0) {
      stop("Treatment indicator must be a logical variable---i.e., TRUE (1) or FALSE (0)")
    }

    nvars  <- ncol(xdata)
    names.xdata  <- names(xdata)

    findx  <- 1
    if (sum(xdata[,1]==rep(1,nrow(xdata)))==nrow(xdata))
      {
        findx  <- 2
      }

    if(nboots < 10 & nboots > 0)
      nboots <- 10
    
    if (ks)
      {
        ks.bm <- KSbootBalanceSummary(index.treated=(Tr==0),
                                      index.control=(Tr==1),
                                      X=xdata[,findx:nvars],
                                      nboots=nboots)

        if (!is.null(match.out))
          {
            ks.am <- KSbootBalanceSummary(index.treated=match.out$index.treated,
                                          index.control=match.out$index.control,
                                          X=xdata[,findx:nvars],
                                          nboots=nboots)
          }
      } 

    BeforeMatchingBalance <- list()
    AfterMatchingBalance <- list()

    BMsmallest.p.value <- 1
    BMsmallest.number <- 1
    BMsmallest.name <- names.xdata[findx]

    AMsmallest.p.value <- NULL
    AMsmallest.number <- NULL
    AMsmallest.name <- NULL
    
    if (!is.null(match.out))
      {
        AMsmallest.p.value <- 1
        AMsmallest.number <- 1
        AMsmallest.name <- names.xdata[findx]        
      }

    for( i in findx:nvars)
      {
        count <- i-findx+1
        if(print.level > 0)
          cat("\n***** (V",count,") ", names.xdata[i]," *****\n",sep="")
        
        ks.do  <- FALSE
        is.dummy  <- length(unique( xdata[,i] )) < 3
        if (ks & !is.dummy)
          ks.do  <- TRUE

        BeforeMatchingBalance[[count]]  <-  balanceUV(xdata[,i][Tr==1], xdata[,i][Tr==0], nboots=0,
                                                      weights.Tr=weights[Tr==1], weights.Co=weights[Tr==0],
                                                      match=FALSE)
        
        if (BeforeMatchingBalance[[count]]$tt$p.value < BMsmallest.p.value)
          {
            BMsmallest.p.value <- BeforeMatchingBalance[[count]]$tt$p.value
            BMsmallest.number <- count
            BMsmallest.name <- names.xdata[i]            
          } else if (BeforeMatchingBalance[[count]]$tt$p.value == BMsmallest.p.value)
            {
              BMsmallest.number <- c(BMsmallest.number,count)
              BMsmallest.name <- c(BMsmallest.name,names.xdata[i])
            }
        
        if (ks.do)
          {
            BeforeMatchingBalance[[count]]$ks <- list()
            BeforeMatchingBalance[[count]]$ks$ks <- list()
            BeforeMatchingBalance[[count]]$ks$ks$p.value <- ks.bm$ks.naive.pval[count]
            BeforeMatchingBalance[[count]]$ks$ks$statistic <- ks.bm$ks.stat[count]              
            if (nboots > 0)
              {
                BeforeMatchingBalance[[count]]$ks$ks.boot.pvalue <- ks.bm$ks.boot.pval[count]

                if (ks.bm$ks.boot.pval[count] < BMsmallest.p.value)
                  {
                    BMsmallest.p.value <- ks.bm$ks.boot.pval[count]
                    BMsmallest.number <- count
                    BMsmallest.name <- names.xdata[i]            
                  } else if ( (ks.bm$ks.boot.pval[count] == BMsmallest.p.value) & (sum(BMsmallest.number==count)==0) )
                    {
                      BMsmallest.number <- c(BMsmallest.number,count)
                      BMsmallest.name <- c(BMsmallest.name,names.xdata[i])
                    }
              } else {
                BeforeMatchingBalance[[count]]$ks$ks.boot.pvalue <- NA

                if (ks.bm$ks.naive.pval[count] < BMsmallest.p.value)
                  {
                    BMsmallest.p.value <- ks.bm$ks.naive.pval[count]
                    BMsmallest.number <- count
                    BMsmallest.name <- names.xdata[i]            
                  } else if ( (ks.bm$ks.naive.pval[count] == BMsmallest.p.value) & (sum(BMsmallest.number==count)==0) )
                    {
                      BMsmallest.number <- c(BMsmallest.number,count)
                      BMsmallest.name <- c(BMsmallest.name,names.xdata[i])
                    }              
              }
            
          } else {
            BeforeMatchingBalance[[count]]$ks <- NULL
          }
        
        if (!is.null(match.out))
          {
            AfterMatchingBalance[[count]]  <- balanceUV(xdata[,i][match.out$index.treated],
                                                        xdata[,i][match.out$index.control],
                                                        weights=match.out$weights, nboots=0,
                                                        paired=paired, match=TRUE)
            
            if (AfterMatchingBalance[[count]]$tt$p.value < AMsmallest.p.value)
              {
                AMsmallest.p.value <- AfterMatchingBalance[[count]]$tt$p.value
                AMsmallest.number <- count
                AMsmallest.name <- names.xdata[i]            
              } else if ( (AfterMatchingBalance[[count]]$tt$p.value == AMsmallest.p.value) & (sum(AMsmallest.number==count)==0) )
                    {
                      AMsmallest.number <- c(AMsmallest.number,count)
                      AMsmallest.name <- c(AMsmallest.name,names.xdata[i])
                    }
            
            if (ks.do)
              {                
                AfterMatchingBalance[[count]]$ks <- list()
                AfterMatchingBalance[[count]]$ks$ks <- list()
                AfterMatchingBalance[[count]]$ks$ks$p.value <- ks.am$ks.naive.pval[count]
                AfterMatchingBalance[[count]]$ks$ks$statistic <- ks.am$ks.stat[count]
                if (nboots > 0)
                  {
                    AfterMatchingBalance[[count]]$ks$ks.boot.pvalue <- ks.am$ks.boot.pval[count]

                    if (ks.am$ks.boot.pval[count] < AMsmallest.p.value)
                      {
                        AMsmallest.p.value <- ks.am$ks.boot.pval[count]
                        AMsmallest.number <- count
                        AMsmallest.name <- names.xdata[i]            
                      } else if ( (ks.am$ks.boot.pval[count] == AMsmallest.p.value) & (sum(AMsmallest.number==count)==0) )
                        {
                          AMsmallest.number <- c(AMsmallest.number,count)
                          AMsmallest.name <- c(AMsmallest.name,names.xdata[i])
                        }
                  } else {
                    AfterMatchingBalance[[count]]$ks$ks.boot.pvalue <- NA

                    if (ks.am$ks.naive.pval[count] < AMsmallest.p.value)
                      {
                        AMsmallest.p.value <- ks.am$ks.naive.pval[count]
                        AMsmallest.number <- count
                        AMsmallest.name <- names.xdata[i]            
                      } else if ( (ks.am$ks.naive.pval[count] == AMsmallest.p.value) & (sum(AMsmallest.number==count)==0) )
                        {
                          AMsmallest.number <- c(AMsmallest.number,count)
                          AMsmallest.name <- c(AMsmallest.name,names.xdata[i])
                        }              
                  }                    
              } else {
                AfterMatchingBalance[[count]]$ks <- NULL
              }
            if(print.level > 0)
              PrintBalanceUV(BeforeMatchingBalance[[count]], AfterMatchingBalance[[count]], digits=digits)
          } else {
            if(print.level > 0)
              {
                cat("before matching:\n")
                summary(BeforeMatchingBalance[[count]], digits=digits)
              }
          } #end of if match.out
      } #end of for loop
    
    if(print.level & ( (nvars-findx+1) > 1))
      {
        cat("\n")

        if (BMsmallest.p.value < 1)
          {
            cat("Before Matching Minimum p.value:", format.pval(BMsmallest.p.value, digits=digits),"\n")
            cat("Variable Name(s):",BMsmallest.name, " Number(s):",BMsmallest.number,"\n\n")
          } else {
            cat("Before Matching Minimum p.value: 1\n\n")
          }

        if (!is.null(match.out))
          {
            if(AMsmallest.p.value < 1)
              {
                cat("After Matching Minimum p.value:", format.pval(AMsmallest.p.value, digits=digits),"\n")
                cat("Variable Name(s):",AMsmallest.name, " Number(s):",AMsmallest.number,"\n\n")
              } else {
                cat("After Matching Minimum p.value: 1\n\n")
              }
          } #end of !is.null(match.out)
      }#end of print.level & (nvars > 1)
  
    return(invisible(list(BeforeMatching=BeforeMatchingBalance,
                          AfterMatching=AfterMatchingBalance,
                          BMsmallest.p.value=BMsmallest.p.value,
                          BMsmallestVarName=BMsmallest.name,
                          BMsmallestVarNumber=BMsmallest.number,
                          AMsmallest.p.value=AMsmallest.p.value,
                          AMsmallestVarName=AMsmallest.name,
                          AMsmallestVarNumber=AMsmallest.number)))
  } #end of MatchBalance

#these are never called if xdata and Tr are provided
get.xdata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "term.labels"))==0 & attr(t1, "intercept")==0) {
    m <- NULL;  # no regressors specified for the model matrix
  } else {
    m <- model.matrix(formul, data=datafr, drop.unused.levels = TRUE)
  }
  return(m);
}


# get.ydata:
# Return response vector corresponding to the formula in formul
# 
get.ydata <- function(formul, datafr) {
  t1 <- terms(formul, data=datafr);
  if (length(attr(t1, "response"))==0) {
    m <- NULL;  # no response variable specified
  }  else {
    m <- model.response(model.frame(formul, data=datafr))
  }
  return(m);
}

