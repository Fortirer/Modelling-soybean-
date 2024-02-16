setwd("C:/Users/Dell/Downloads")

# instale os pacotes se nao tiver
#install.packages(c("tidyverse", "agricolae"))
install.packages("car")
library(car)

# chamando bibliotecas
library(tidyverse)
library(agricolae)

dados <- read.csv("N_Temp.csv", sep=";")

ggplot(dados, aes(x = Treatment, y = N)) +
  geom_boxplot() +
  labs(title = "Boxplot dos Tratamentos")

# ANOVA
# importante verificar as premissas do teste ANOVA (normalidade dos residuos, homocedasticidade da variancia, independencia - residuos)
modelo_anova <- aov(N ~ Treatment, data = dados)
anova_resultado <- summary(modelo_anova)

print(anova_resultado)

# VERIFICANDO PREMISSAS TESTE ANOVA

# modelo_anova <- aov(Valor ~ Tratamento, data = dados)

# Normalidade dos Resíduos
residuos <- residuals(modelo_anova)

# Histograma dos Resíduos # PROBLEMA  N baixo não da para ver bem a normalidade
hist(residuos, main = "Histograma dos Resíduos", col = "lightblue")

# Gráfico Q-Q dos Resíduos
qqnorm(residuos)
qqline(residuos)

# Homogeneidade das Variâncias
# Gráfico dos Resíduos versus Valores Ajustados
plot(modelo_anova, 1)

# Independência dos Resíduos
# Gráfico dos Resíduos versus Observações
plot(modelo_anova, 2)

# para  N baixo fica dificil ver visualmente os teste das premissas do test ANOVA
# fazendo teste shapiro

# Teste de Shapiro-Wilk para normalidade
shapiro.test(residuos)   # hipotese nula eh que os dados tem normalidade valor abaixo 0.05 sugere que os dados nao tem normaidade

# Teste de Levene para homogeneidade de variâncias
leveneTest(modelo_anova)  #  hipotese nula eh de homocedasticidade dos residuos

# caso os dados não cumpra com as premissas, pode-se fazer a normalização dos dados / padronização

# feito isso, bora seguir  

# TUKEY TEST
tukey_resultado <- HSD.test(modelo_anova, "Treatment", group = TRUE)
print(tukey_resultado)

##### outras possibilidades

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
