# Z'-factor
# http://en.wikipedia.org/wiki/Z-factor
# Zhang, J.-H., Chung, T., & Oldenburg, K. (1999). A Simple Statistical 
# Parameter for Use in Evaluation and Validation of High Throughput Screening 
# Assays. Journal of Biomolecular Screening, 4(2), 67–73.

f_zprime <- function(x,y)  1 - 3 * (sd(x) + sd(y)) / abs(mean(x)-mean(y))

# One tailed Z'-factor
# http://www.cellprofiler.org/CPmanual/CalculateStatistics.html
# OneTailedZfactor: This measure is an attempt to overcome a limitation of the 
# original Z'-factor formulation (it assumes a Gaussian distribution) and is 
# informative for populations with moderate or high amounts of skewness. In 
# these cases, long tails opposite to the mid-range point lead to a high 
# standard deviation for either population, which results in a low Z' factor #
# even though the population means and samples between the means may be 
# well-separated. 
# Therefore, the one-tailed Z' factor is calculated with the same formula but 
# using only those samples that lie between the positive/negative population 
# means.

mirror_mean <- function(x) {
  m <- mean(x)
  x <- x - m
  x <- x[x>0]
  x <- c(m-x,m+x)
  x
}

f_zprime_os <- function(x,y) {
  if (mean(x) > mean(y)) {
    x <- -x
    y <- -y
  }
  f_zprime(mirror_mean(x), -mirror_mean(-y))  
}


