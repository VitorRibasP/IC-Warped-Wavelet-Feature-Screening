import numpy as np #Pacote para Operações
import matplotlib.pyplot as plt #Pacote para Plotar Gráficos
import pywt #Pacote de Wavalets
import pandas as pd #Pacote de Bancos de Dados
import math #Outro pacote para operações
from skimage.restoration import denoise_wavelet #VisuShrink e BayesShrink
import threshold
import ops
import denoising as dsg

#Função que faz um Threshold a partir de um parâmetro
def Threshold(Dados, Wavelet, Modo, Parametro):
    Coeffs = pywt.wavedec(Dados, Wavelet, level = pywt.dwt_max_level(len(Dados), Wavelet))
    for i in range(len(Coeffs)-1):
        for j in range(len(Coeffs[i+1])):
            if Modo == "soft":
                Coeffs[i+1][j] = max(Coeffs[i+1][j]/np.abs(Coeffs[i+1][j]) * np.abs(Coeffs[i+1][j]) - Parametro, 0)
            if Modo == "hard":
                if np.abs(Coeffs[i+1][j]) <= Parametro:
                    Coeffs[i+1][j] = 0
    return pywt.waverec(Coeffs, Wavelet)

#Função que plota um gráfico de pontos de um Conjunto de Dados
def pontos_dados(Dados, title):
    l = list(range(1, len(Dados)+1))
    #plt.figure()
    plt.scatter(l, Dados)
    plt.title(title)
    plt.ylabel("Valor")

#Função que plota um gráfico de uma Função (conecta os pontos)
def funcao_dados(Funcao, title):
    #plt.figure()
    plt.plot(Funcao)
    plt.title(title)

#Função que plota o Gráfico com a Original e Reconstruída
def funcoes(Original, Reconstruida, Título):
    plt.plot(Original, color="orange", label="Original")
    plt.plot(Reconstruida, color="black", label="Reconstruída")
    plt.title(Título)
    plt.legend() 

#Função que plota o gráfico dos Coeficientes por Resolução em linhas
def coeficientes_linhas(Dados, Wavelet):
    L = pywt.dwt_max_level(len(Dados), Wavelet)
    coeffs = pywt.wavedec(Dados, Wavelet, level = L)
    plt.figure()
    for i in range(L+1):
        plt.subplot(L+1,1,L-i+1)
        l = list(range(1, len(coeffs[L-i])+1))
        plt.stem(l, coeffs[L-i], basefmt='w', markerfmt='', linefmt='')
        plt.axis('off')
        plt.title("{}".format(L-i-1), loc='left', y = 0.1)
    plt.suptitle("Nível de Resolução", x=0.05, y=0.75, fontsize=16, ha='left', rotation='vertical')

def criacao_even(Dados):
    Even = []
    for i in range(len(Dados)):
        if i%2 != 0:
            Even.append(Dados[i])
    return Even

def criacao_odd(Dados):
    Odd = []
    for i in range(len(Dados)):
        if i%2 == 0:
            Odd.append(Dados[i])
    return Odd

def cross_validation(Dados, Wavelet, Modo):
    #Separando o Banco em Dois
    Even = []
    Odd = []
    for i in range(len(Dados)):
        if i%2 == 0:
            Odd.append(Dados[i])
        else:
            Even.append(Dados[i])
    #Criando as funções y
    yOdd = []
    for i in range(int((len(Dados)/2))-1):
        yOdd.append((Odd[i]+Odd[i+1])/2)
    yOdd.append((Odd[i]+Odd[-1])/2)
    yEven = []
    for i in range(int((len(Dados)/2))-1):
        yEven.append((Even[i]+Even[-1])/2)
    yEven.append((Even[i]+Even[-1])/2)
    #Calculando todos os M(t)
    EvenCoeffs = pywt.wavedec(Even, Wavelet, level = pywt.dwt_max_level(len(Even), Wavelet))
    OddCoeffs = pywt.wavedec(Odd, Wavelet, level = pywt.dwt_max_level(len(Odd), Wavelet))
    Coeffs = []
    for Coeff in EvenCoeffs:
        for coeff in Coeff:
            Coeffs.append(coeff)
    for Coeff in OddCoeffs:
        for coeff in Coeff:
            Coeffs.append(coeff)
    t = 0
    Min = "Início"
    for coeff in Coeffs:
        M = 0
        fEven = Threshold(Even, Wavelet, Modo, np.abs(coeff))
        fOdd = Threshold(Odd, Wavelet, Modo, np.abs(coeff))
        for j in range(len(Odd)):
            M = M + (fEven[j]-yOdd[j])**2 + (fOdd[j]-yEven[j])**2   
        if Min == "Início" or M < Min:
            Min = M
            t = np.abs(coeff)
    return t

def hypotesis_testP5(Dados, Wavelet, Modo):
    Coeffs = pywt.wavedec(Dados, Wavelet, level = pywt.dwt_max_level(len(Dados), Wavelet))
    for coeff in Coeffs:
        t = 0
        for i in range(len(coeff)):
            if (coeff[i])**2 < 1.6449 and np.abs(coeff[i]) > t:
                t = np.abs(coeff[i])
        for i in range(len(coeff)-1):
            if Modo == "soft":
                coeff[i+1] = (coeff[i+1]/np.abs(coeff[i+1])) * np.maximum(np.abs(coeff[i+1]) - t, 0)
            if Modo == "hard":
                if np.abs(coeff[i+1]) <= t:
                    coeff[i+1] = 0
    return pywt.waverec(Coeffs, Wavelet)

#Dados
t = np.linspace(0,1,1024)
Z = np.sqrt(t*(1-t))*np.sin((2.1*np.pi)/(t+0.05))
Xa = Z + np.random.normal(0,.12,1024)

x = np.linspace(0, 1, 1024)
X = np.sin(2*np.pi*x)
X = X + np.random.normal(0,0.5, 1024)

#denoise_wavelet(Dados, method = "VisuShrink", mode = "hard", wavelet = "db2")
#dsg.tsd(Dados, method = "sureshrink", mode = "soft", wavelets_name = "db2")

#VisuShrink, SureShrink, Minimax, BayesShrink, CrossValidation, HypotesisTest
#plt.plot(uniform_design(X, "db6", "hard"))
#plt.show()

plt.subplot(1,2,1)
plt.plot(Z)
plt.title("Doppler Original")
plt.subplot(1,2,2)
plt.plot(Xa)
plt.title("Doppler com Ruído")
plt.show()

Modo = "soft"
Wavelet = "db13"

plt.subplot(2,3,1)
plt.plot(denoise_wavelet(Xa, method = "VisuShrink", mode = Modo, wavelet = Wavelet))
plt.title("Visu Shrink")
plt.subplot(2,3,2)
plt.plot(dsg.tsd(Xa, method = "sureshrink", mode = Modo, wavelets_name = Wavelet))
plt.title("Sure Shrink")
plt.subplot(2,3,3)
plt.plot(dsg.tsd(Xa, method = "minmax", mode = Modo, wavelets_name = Wavelet))
plt.title("Minimax")
plt.subplot(2,3,4)
plt.plot(denoise_wavelet(Xa, method = "BayesShrink", mode = Modo, wavelet = Wavelet))
plt.title("Bayes Shrink")
plt.subplot(2,3,5)
plt.plot(Threshold(Xa, Wavelet, Modo, cross_validation(Xa, Wavelet, Modo)))
plt.title("Cross Validation")
plt.subplot(2,3,6)
plt.plot(hypotesis_testP5(Xa, Wavelet, Modo))
plt.title("Hypotesis Test")
plt.show()
