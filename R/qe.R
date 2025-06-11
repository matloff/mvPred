

# allows missing-value handling for most predictive functions in the 
# qeML package

# arguments:

# data,yName,holdout: as in qeML package
# mvFtn: function for handling NAs, e.g. 'mice'
# qeMLftn: qeML predictive function, e.g. 'qeRF' (random forests)
# qeMLopts: optional argments for 'qeMLftn'
# mvPredOpts: optional orguments for 'mvFtn', e.g. 'm' in 'mice'

qeMLna <- function(data,yName,
    qeMLftn,mvFtn,qeMLopts,mvPredOpts,holdout=1000)
{

   # cases in which MV method involves pre-processing the data
   if (mvFtn == 'mice') print()

}
