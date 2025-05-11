
# allows missing-value handling for most predictive functions in the 
# qeML package

# arguments:

# data,yName,holdout: as in qeML package
# mvFtn: function for handling NAs, e.g. 'mice'
# qeMLftn: qeML predictive function, e.g. 'qeRF' (random forests)
# qeMLopts: optional argments for 'qeMLftn'
# mvPredOpts: optional orguments for 'mvFtn', e.g. 'm' in 'mice'

qeMLna <- function(data,yName,qeMLftn,
   mvFtn,qeMLopts=NULL,mvPredOpts=NULL,seed=9999,holdout=1000)
{
   require(qeML)
   # cases in which MV method involves pre-processing the data
   if (mvFtn == 'compCases') {
      tmp <- compCases(data)
      data <<- tmp$data
      intactRows <<- tmp[[2]]
      nonintactRows <<- tmp[[3]]
   } else if (mvFtn == 'mice') {
   }

   qeMLout <- do.call(qeMLftn,qeMLopts)
   qeMLout$intactRows <- intactRows
   qeMLout$nonintactRows <- nonintactRows

   qeMLout
}

