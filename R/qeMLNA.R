
# allows missing-value handling for most predictive functions in the 
# qeML package

# arguments:

# data,yName,holdout: as in qeML package
# mvFtn: function for handling NAs, e.g. 'mice'
# qeMLftn: qeML predictive function, e.g. 'qeRF' (random forests)
# qeMLopts: optional argments for 'qeMLftn'
# mvPredOpts: optional orguments for 'mvFtn', e.g. 'm' in 'mice'
# retainMVFtnOut: if TRUE, will be a component in the output
# seed: for random numbers in choosing holdout set
# holdout: size of holdout set (can be 0)

# examples:

# (load 'english', from 'data/'; after make this a package, just use 
# 'data(english)

# eng <-
#    english[,c("age","birth_order","ethnicity","sex","mom_ed","vocab")] 

# z <- qeMLna(eng,'vocab','qeLin','compCases',retainMVFtnOut=FALSE)
# z <- qeMLna(eng,'vocab','qeLin','mice',retainMVFtnOut=FALSE)

# just type 'z' to see what's in there, e.g. the estimated beta
# coefficients based on the non-NA data

# Not done!!  Need to define a 'predict' method for the class
# 'MVqeMLout'

qeMLna <- function(data,yName,qeMLftn,
   mvFtn,qeMLopts=NULL,mvPredOpts=NULL,retainMVFtnOut=TRUE,
   seed=9999,holdout=1000)
{
   require(qeML)
   # cases in which MV method involves pre-processing the data
   if (mvFtn == 'compCases') {
      tmp <- compCases(data)
      intactRows <<- tmp[[2]]
      nonintactRows <<- tmp[[3]]
      qeMLopts <- addListElement(qeMLopts,'data',tmp$intactData)
   } else if (mvFtn == 'mice') {
      warning('m will be set to 1')
      mvPredOpts[['data']] <- data
      mvPredOpts[['m']] <- 1
      MVFtnOut <- do.call('mice',mvPredOpts)
      imp <- complete(MVFtnOut)
      qeMLopts <- addListElement(qeMLopts,'data',imp)
   }

   qeMLopts <- addListElement(qeMLopts,'yName',yName)

   qeMLout <- do.call(qeMLftn,qeMLopts)

   qeMLout$intactRows <- intactRows
   qeMLout$nonintactRows <- nonintactRows
   if (retainMVFtnOut) qeMLout$MVFtnOut <- MVFtnOut

   class(qeMLout) <- 'MVqeMLout'

   qeMLout
}

# adds/replaces element of the given name in the list 'l'
addListElement <- function(l,newEltName,newElt) 
{
   if (is.null(l)) l <- list()
   l[[newEltName]] <- newElt
   l
}
