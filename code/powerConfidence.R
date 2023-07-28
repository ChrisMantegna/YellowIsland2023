# testing statistical power and confidence intervals for data analysis

# Load the 'pwr' package (install it if you haven't already)
install.packages("pwr")
library(pwr)

# Enter the effect size (Cohen's d), sample sizes, and alpha level
effect_size <- 0.2  # Replace this with your desired effect size
n1 <- 166 # Sample size of group 1
n2 <- 110 # Sample size of group 2
n3 <- 110 # Sample size of group 3
alpha <- 0.05

# Calculate the statistical power

power <- pwr.anova.test(k = 3, n = c(n1, n2, n3), f = effect_size, sig.level = alpha, power = NULL)

# Print the statistical power
print(power)
