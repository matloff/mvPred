
# Complete Cases method; throw out any row that has at least 1 NA

# returns the intact rows, and indices of intact and nonintact rows

compCases <- function(data) 
{
   intactRows <- complete.cases(data)
   list(intactData=data[intactRows,],intactRows=intactRows,
      nonintactRows=setdiff(1:nrow(data),intactRows))
}

