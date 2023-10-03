library(tidyverse)

# calculate NNS
set.seed(428)
stage.weights <- c(0.547, 0.075, 0.218, 0.16)
stage.weights <- stage.weights / sum(stage.weights)

detected <- c(65, 31, 51, 60)
total <- c(92, 35, 58, 61)

fp <- 71
nc <- 134

S <- 100000
results <- matrix(0, 5, S)
  
results[1, ] <- rbinom(S, nc, fp/nc)/nc
for (n in 2:5) {
    results[n, ] <- rbinom(S, total[n-1], detected[n-1]/total[n-1])/total[n-1]
}

specificity <- results[1, ]
sensitivity <- stage.weights %*% results[2:5, ]

npv <- function(se, sp, prev) {
    (sp * (1 - prev)) / ((sp * (1 - prev)) + ((1 - se) * prev))
}

ppv <- function(se, sp, prev) {
    (se * prev) / ((se * prev) + ((1 - sp) * (1 - prev)))
}

RR <- ppv(sensitivity, specificity, 0.007) / 
    (1-npv(sensitivity, specificity, 0.007))

quantile(RR, c(0.025, 0.5, 0.975))

1/ppv(stage.weights %*% (detected / total), 71 / 134, 0.007)
1/(1-npv(stage.weights %*% (detected / total), 71 / 134, 0.007))

ppv(stage.weights %*% (detected / total), 71 / 134, 0.007)/(1-npv(stage.weights %*% (detected / total), 71 / 134, 0.007))

NNS.positive <- 1/ppv(sensitivity, specificity, 0.007)
NNS.negative <- 1/(1-npv(sensitivity, specificity, 0.007))

mean(NNS.positive)
mean(NNS.negative)

NNS.positive.ul <- quantile(NNS.positive, 0.975)
NNS.negative.ul <- quantile(NNS.negative, 0.975)


# figure 4d
data <- tibble(result=c("Positive Test Result",
                        "Negative Test Result"),
               nns=c(84, 381),
               ll=c(81, 375),
               ul=c(87, 390))

pdf("fig4d.pdf", height=4)
ggplot(data, aes(x=result, y=nns)) +
    geom_bar(fill="grey", stat="identity", width=0.5) +
    #geom_segment(aes(x=c(1), y=381, xend=c(1), yend=NNS.negative.ul)) +
    #geom_segment(aes(x=c(0.875), y=NNS.negative.ul, xend=c(1.125), 
                     #yend=NNS.negative.ul)) +
    #geom_segment(aes(x=c(2), y=84, xend=c(2), yend=NNS.positive.ul)) +
    #geom_segment(aes(x=c(1.875), y=NNS.positive.ul, xend=c(2.125), 
                     #yend=NNS.positive.ul)) +
    geom_segment(aes(x=2, y=84, xend=2, yend=425), colour="lightblue") +
    geom_segment(aes(x=1, y=381, xend=1, yend=425), colour="lightblue") +
    geom_segment(aes(x=2, y=425, xend=1.7, yend=425), colour="lightblue") +
    geom_segment(aes(x=1, y=425, xend=1.3, yend=425), colour="lightblue") +
    geom_text(aes(label=nns),
              colour="white", hjust=1.5, size=6) +
    geom_text(label="RR: 4.5\n(95% CI: 2.8-7.5)",
              aes(x=1.5, y=425), size=5) +
    theme_classic(base_size=11) +
    theme(axis.text=element_text(size=10)) +
    scale_y_continuous(n.breaks=6) +
    scale_fill_manual(values=c("grey", "lightblue", "darkblue")) +
    expand_limits(y=c(0, 500)) +
    coord_flip() +
    labs(x="", y="Number Needed to Screen (NNS) by LDCT to detect one lung cancer")
dev.off()

