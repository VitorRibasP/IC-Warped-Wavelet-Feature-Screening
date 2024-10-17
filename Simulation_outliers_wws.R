#A função numSelectFeat_wws foi criada no artigo 
#O link para o artigo é: https://www.tandfonline.com/doi/full/10.1080/10618600.2024.2342984?scroll=top&needAccess=true#abstract

#Para essa simulação, esse mesmo código foi rodado 2 vezes
#Na primeira, foi considerado a função original
#Na segunda, foi feita uma alteração na linha 60 do código do artigo para Y <- Y + rt(n, mean = 0, sd = error_std)
#Dessa forma, a primeira vez foi considerando uma distribuição normal
#E a segunda vez considerando uma distribuiçaõ t.

Medias <- c()

for (i in 1:15) {
  Medias <- c(Medias, mean(numSelectFeat_wws[((c(0:39)*15+i)[-(c(1:10)*4)])]))
}

Tabela <- data.frame(Blocks = Medias[c(1,2,3)],
                     Heavisine = Medias[c(4,5,6)],
                     Bumps = Medias[c(7,8,9)],
                     Doppler = Medias[c(10,11,12)],
                     Identity = Medias[c(13,14,15)])
Tabela <- round(Tabela, 2)
rownames(Tabela) <- c("64", "128", "256")
