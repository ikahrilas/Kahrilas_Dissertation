# sensitivity power analyses
library(pwr)
library(WebPower)

## relation between depression and psychological variables (n = 132)
pwr.r.test(n = 132,
           sig.level = .05,
           power = .80,
           r = NULL,
           alternative = "greater")
r
plot(r)
## relation between depression and emotion reactivity (N = 107) with 3 repeated
## measures for each participant (n = 321)
wp.rmanova(n = 321, ng = 107, nm = 3, f = .48, nscor = 1,
           alpha = 0.05, power = NULL, type = 2)

## relation between headspace use and emotion reactivity (N = 55) with 6 repeated measures
## for each participant (n = 330)
wp.rmanova(n = 330, ng = 55, nm = 6, f = .46, nscor = 1,
           alpha = 0.05, power = NULL, type = 2)
