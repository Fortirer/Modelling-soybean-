# tukey test

setwd("C:/Users/Dell/Downloads")

data <- read.csv("NitroCarbTemp.csv", sep=";")
data$tratamento <- as.factor(data$tratamento)

data.aov <- aov(X.N ~ tratamento, data)
data.tukey <- glht (data.aov, linfct = mcp (tratamento = "Tukey"))
cld(data.tukey)


# outra opcao
mod <- lm(X.N ~ tratamento, data = data)
mod_means_contr <- emmeans::emmeans(object = mod, pairwise ~ "tratamento", adjust = "tukey")

mod_means <- multcomp::cld(object = mod_means_contr$emmeans,
                           Letters = letters)
library(ggplot2)

ggplot(data = mod_means,
       aes(x = tratamento, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, 
                    ymax = upper.CL), 
                width = 0.2) +
  geom_point() +
  geom_text(aes(label = gsub(" ", "", .group)),
            position = position_nudge(x = 0.2)) +
  labs(caption = "Means followed by a common letter are\nnot significantly different according to the Tukey-test")