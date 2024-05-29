# TCC_QCD2
Códigos para gerar os autovalores do operador de massa quadrática.

- Espectro.py:
Constrói a base de estados, calcula os elementos de matriz de M^2, encontra seus autovalores e imprime um gráfico do espectro. Também faz uma extrapolação da massa do férmion e bóson mais leves.
Parâmetros: K_min e K_max ditam o menor K e o maior K para calcular os autovalores de M^2. O parâmetro adimensional 'y' controla a intensidade da interação de gauge.

- Espectro_supersimetria.py:
Constrói a base de estados, calcula os elementos de matriz de M^2, encontra seus autovalores e imprime um gráfico do espectro. Também faz uma extrapolação da massa dos dois hádrons supersimétricos mais leves.
Parâmetros: K_min e K_max ditam o menor K e o maior K para calcular os autovalores de M^2. O parâmetro adimensional 'y' controla a intensidade da interação de gauge.

- Gerador_bases.py:
Constrói a base de estados para do subespaço com momento discreto K e imprime a dimensão deste.
Parâmetros: o momento discreto K
